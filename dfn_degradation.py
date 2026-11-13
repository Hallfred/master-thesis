"""Coupled O'Kane (2022) degradation and voltage-constrained cell dispatch.

The hourly LP proposes energy transfers. A nonlinear DFN prediction checks each
shorter constant-current action before it is committed. Rejected predictions
always restart from the same accepted electrochemical state.
"""

from dataclasses import asdict, replace
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
import pybamm
from scipy.integrate import trapezoid


SOURCE = "https://docs.pybamm.org/en/v25.12.2/source/examples/notebooks/models/coupled-degradation.html"
DEGRADATION_OPTIONS = {
    "SEI": "solvent-diffusion limited",
    "SEI porosity change": "true",
    "lithium plating": "partially reversible",
    "lithium plating porosity change": "true",
    "particle mechanics": ("swelling and cracking", "swelling only"),
    "SEI on cracks": "true",
    "loss of active material": "stress-driven",
}
VARIABLES = {
    "voltage_v": "Terminal voltage [V]",
    "temperature_k": "Volume-averaged cell temperature [K]",
    "negative_stoichiometry": "Negative electrode stoichiometry",
    "throughput_ah": "Throughput capacity [A.h]",
    "sei_loss_ah": "Loss of capacity to negative SEI [A.h]",
    "sei_cracks_loss_ah": "Loss of capacity to negative SEI on cracks [A.h]",
    "plating_loss_ah": "Loss of capacity to negative lithium plating [A.h]",
    "side_reaction_loss_ah": "Total capacity lost to side reactions [A.h]",
    "lli_pct": "Loss of lithium inventory [%]",
    "lam_negative_pct": "Loss of active material in negative electrode [%]",
    "lam_positive_pct": "Loss of active material in positive electrode [%]",
    "negative_porosity": "X-averaged negative electrode porosity",
}


def validate_config(c, cfg):
    if cfg.model not in ("DFN", "SPMe") or cfg.parameter_set != "OKane2022":
        raise ValueError("Use DFN or SPMe with the coupled OKane2022 parameter set")
    for name in ("hours", "series_cells", "output_period_s", "control_period_s",
                 "mesh_points", "particle_mesh_points", "current_search_iterations", "rest_prediction_s"):
        value = getattr(cfg, name)
        if not np.isfinite(value) or value < 1 or int(value) != value:
            raise ValueError(f"{name} must be a positive integer")
    if cfg.parallel_cells is not None and (cfg.parallel_cells < 1 or int(cfg.parallel_cells) != cfg.parallel_cells):
        raise ValueError("parallel_cells must be a positive integer")
    if not 0 < c.soc_min < c.initial_soc < min(c.soc_max, cfg.soc_ceiling) < 1:
        raise ValueError("Cell study requires SOC min < initial SOC < effective SOC max < 1")
    if c.max_c_rate <= 0 or not np.isfinite(c.max_c_rate):
        raise ValueError("max_c_rate must be positive and finite")
    if not (0 < c.charge_efficiency <= 1 and 0 < c.discharge_efficiency <= 1):
        raise ValueError("Invalid efficiency")
    if (cfg.thermal not in ("lumped", "isothermal") or cfg.cooling_w_m2_k <= 0
            or not np.isfinite(cfg.cooling_w_m2_k)):
        raise ValueError("Use lumped or isothermal thermal dynamics and positive cooling")
    if (not np.isfinite([cfg.voltage_margin_v, cfg.rtol, cfg.atol, cfg.cell_nominal_voltage_v]).all()
            or min(cfg.voltage_margin_v, cfg.rtol, cfg.atol, cfg.cell_nominal_voltage_v) <= 0):
        raise ValueError("Voltage margin, solver tolerances and nominal voltage must be positive")
    if not np.isfinite(cfg.voltage_headroom_v) or not 0 < cfg.voltage_headroom_v < cfg.voltage_margin_v:
        raise ValueError("Voltage headroom must be positive and smaller than the voltage margin")
    if not cfg.particle_radius_factors or any(not np.isfinite(f) or f <= 0 for f in cfg.particle_radius_factors):
        raise ValueError("Particle radius factors must be positive and finite")
    if len(set(cfg.particle_radius_factors)) != len(cfg.particle_radius_factors):
        raise ValueError("Particle radius factors must be distinct")


class DegradationPlant:
    """One representative cell in a balanced Ns x Np pack, with persistent state."""

    def __init__(self, c, cfg, pack_ah, ambient_c, radius_factor=1.0):
        validate_config(c, cfg)
        if not np.isfinite(pack_ah) or pack_ah <= 0 or not np.isfinite(radius_factor) or radius_factor <= 0:
            raise ValueError("Pack capacity and radius factor must be positive and finite")
        self.case, self.cfg = c, cfg
        self.parameters = pybamm.ParameterValues(cfg.parameter_set)
        self.cell_ah = float(self.parameters["Nominal cell capacity [A.h]"])
        self.np = cfg.parallel_cells or int(np.ceil(pack_ah / self.cell_ah))
        self.ns = cfg.series_cells
        self.nominal_kwh = self.ns * self.np * self.cell_ah * cfg.cell_nominal_voltage_v / 1000
        self.soc_min, self.soc_max = c.soc_min, min(c.soc_max, cfg.soc_ceiling)
        self.v_min = float(self.parameters["Lower voltage cut-off [V]"]) + cfg.voltage_margin_v
        self.v_max = float(self.parameters["Upper voltage cut-off [V]"]) - cfg.voltage_margin_v
        if self.v_min >= self.v_max:
            raise ValueError("Voltage margin leaves no operating voltage range")
        for electrode in ("Negative", "Positive"):
            key = f"{electrode} particle radius [m]"
            self.parameters.update({key: self.parameters[key] * radius_factor})
        # SOC uses electrode stoichiometry, rather than the EMS energy accumulator.
        self.x0, _ = pybamm.lithium_ion.get_initial_stoichiometries(0, self.parameters)
        self.x100, _ = pybamm.lithium_ion.get_initial_stoichiometries(1, self.parameters)
        self.parameters.update({
            "Current function [A]": "[input]",
            "Ambient temperature [K]": "[input]",
            "Initial temperature [K]": float(ambient_c) + 273.15,
            "Total heat transfer coefficient [W.m-2.K-1]": cfg.cooling_w_m2_k,
        })
        options = {**DEGRADATION_OPTIONS, "thermal": cfg.thermal}
        self.model = getattr(pybamm.lithium_ion, cfg.model)(options)
        voltage = self.model.variables["Terminal voltage [V]"]
        soc = (self.model.variables["Negative electrode stoichiometry"] - self.x0) / (self.x100 - self.x0)
        self.model.variables["Controller SOC"] = soc
        # Continuous-time root finding supplements sampled trajectory checks.
        # The original physical voltage/porosity events remain enabled.
        for name, expression in (
            ("Controller lower voltage margin", voltage - self.v_min),
            ("Controller upper voltage margin", self.v_max - voltage),
            ("Controller lower SOC", soc - self.soc_min),
            ("Controller upper SOC", self.soc_max - soc),
        ):
            self.model.events.append(pybamm.Event(name, expression))
        self.solver = pybamm.IDAKLUSolver(rtol=cfg.rtol, atol=cfg.atol)
        self.simulation = pybamm.Simulation(
            self.model, parameter_values=self.parameters, solver=self.solver,
            var_pts={"x_n": cfg.mesh_points, "x_s": cfg.mesh_points, "x_p": cfg.mesh_points,
                     "r_n": cfg.particle_mesh_points, "r_p": cfg.particle_mesh_points},
        )
        # A microsecond rest obtains a consistent initial DAE state. It is excluded
        # from reported elapsed time; no state is reset between dispatch steps.
        initial = self.simulation.solve(
            [0, 1e-6], initial_soc=c.initial_soc, inputs=self.inputs(0, ambient_c))
        self.state = initial.last_state
        self.time_origin = float(self.state.t[-1])
        self.predictions = self.rejected_predictions = 0
        self.last_rejection = ""

    @staticmethod
    def inputs(current_a, ambient_c):
        return {"Current function [A]": float(current_a), "Ambient temperature [K]": float(ambient_c) + 273.15}

    @property
    def soc(self):
        return float(self.state["Controller SOC"].entries[-1])

    @property
    def capacity_proxy_kwh(self):
        # A scheduling proxy only, not a measured discharge capacity.
        losses = [float(self.state[VARIABLES[key]].entries[-1])
                  for key in ("lli_pct", "lam_negative_pct", "lam_positive_pct")]
        return self.nominal_kwh * max(0.01, 1 - max(0, *losses) / 100)

    def _predict(self, current_a, ambient_c, duration_s, starting_state=None):
        self.predictions += 1
        return self.solver.step(
            self.state if starting_state is None else starting_state,
            self.simulation.built_model, duration_s, save=False,
            inputs=self.inputs(current_a, ambient_c),
            t_eval=np.array([0., duration_s]),
            t_interp=np.linspace(0., duration_s, max(2, int(np.ceil(duration_s / self.cfg.output_period_s)) + 1)),
        )

    def _feasible(self, trial, current_a, duration_s, power_limit_kw, starting_time=None):
        start = self.state.t[-1] if starting_time is None else starting_time
        if trial.termination != "final time" or abs(trial.t[-1] - start - duration_s) > 1e-4:
            self.last_rejection = f"termination={trial.termination}, completed_s={trial.t[-1]-start:.6g}"
            return False
        voltage = trial["Terminal voltage [V]"].entries
        soc = trial["Controller SOC"].entries
        power = np.abs(current_a * voltage * self.ns * self.np / 1000)
        # Small SOC reserve permits calendar side reactions during subsequent rest.
        lo = self.soc_min + (1e-4 if current_a > 0 else 0)
        hi = self.soc_max - (1e-4 if current_a < 0 else 0)
        v_lo = self.v_min + (self.cfg.voltage_headroom_v if current_a > 0 else 0)
        v_hi = self.v_max - (self.cfg.voltage_headroom_v if current_a < 0 else 0)
        feasible = bool(np.isfinite(voltage).all() and np.isfinite(soc).all()
                    and voltage.min() >= v_lo and voltage.max() <= v_hi
                    and soc.min() >= lo and soc.max() <= hi
                    and power.max() <= power_limit_kw + 1e-9)
        if not feasible:
            self.last_rejection = (f"V={voltage.min():.9g}..{voltage.max():.9g}, "
                                   f"SOC={soc.min():.9g}..{soc.max():.9g}, "
                                   f"power={power.max():.9g}/{power_limit_kw:.9g} kW")
        return feasible

    def advance(self, requested_dc_kw, ambient_c, duration_s):
        """Predict and commit a safe current, positive for discharge.

        The feasible prefix is found by bounded current search. Even predictions
        that hit an event are discarded in full: time and degradation advance only
        when a complete feasible interval is committed.
        """
        if (not np.isfinite([requested_dc_kw, ambient_c, duration_s]).all() or duration_s <= 0
                or ambient_c <= -273.15):
            raise ValueError("Invalid current request, ambient temperature or duration")
        request = float(requested_dc_kw)
        power_limit = abs(request)
        # Using the most adverse allowed voltage ensures the initial candidate
        # cannot consume more PV or supply more demand than the EMS permits.
        voltage_bound = self.v_max
        current = np.sign(request) * min(
            power_limit * 1000 / (self.ns * self.np * voltage_bound),
            self.case.max_c_rate * self.cell_ah)

        def attempt(value):
            try:
                trial = self._predict(value, ambient_c, duration_s)
                if self._feasible(trial, value, duration_s, power_limit):
                    if value != 0:
                        # Plated lithium can re-intercalate during rest and raise
                        # SOC/OCV even after charging stops. Reject actions that
                        # leave no feasible rest continuation at measured ambient.
                        rest = self._predict(0., ambient_c, self.cfg.rest_prediction_s,
                                             starting_state=trial.last_state)
                        if not self._feasible(rest, 0., self.cfg.rest_prediction_s, 0.,
                                              starting_time=trial.t[-1]):
                            self.rejected_predictions += 1
                            return None
                    return trial
            except pybamm.SolverError as exc:
                self.last_rejection = str(exc)
                self.rejected_predictions += 1
                return None
            self.rejected_predictions += 1
            return None

        accepted = attempt(current)
        accepted_current = current
        limited = False
        if accepted is None:
            limited = True
            accepted = attempt(0.)
            accepted_current = 0.
            if accepted is None:
                raise RuntimeError(
                    f"No feasible DFN rest step at {(self.state.t[-1]-self.time_origin)/3600:.3f} h; "
                    f"{self.last_rejection}")
            low, high = 0., 1.
            for _ in range(self.cfg.current_search_iterations):
                middle = (low + high) / 2
                trial = attempt(current * middle)
                if trial is None:
                    high = middle
                else:
                    low, accepted, accepted_current = middle, trial, current * middle
        frame = pd.DataFrame({key: np.asarray(accepted[name].entries).reshape(-1)
                              for key, name in VARIABLES.items()})
        frame["elapsed_s"] = accepted.t - self.time_origin
        frame["soc"] = accepted["Controller SOC"].entries
        frame["current_a"] = accepted_current
        frame["requested_dc_kw"] = request
        frame["pack_dc_kw"] = accepted_current * frame.voltage_v * self.ns * self.np / 1000
        frame["constraint_limited"] = limited
        self.state = accepted.last_state
        energy = float(trapezoid(frame.pack_dc_kw, frame.elapsed_s) / 3600)
        return frame, energy


def simulate_electrochemical(c, profiles, pv_kwp, pack_ah, cfg, mpc_config, radius_factor=1.0):
    """Closed-loop MPC with DFN degradation for any consecutive hourly window.

    Unlike a current-profile replay, each LP uses the accepted cell SOC and
    degradation state. Its first interval shrinks as the hour progresses.
    """
    from pv_battery_model import mpc_action

    validate_config(c, cfg)
    required = ["pv_kw_per_kwp", "load_kw", "ambient_c"]
    if (profiles.empty or not np.isfinite(profiles[required]).all().all()
            or (profiles[["pv_kw_per_kwp", "load_kw"]] < 0).any().any()
            or not np.isfinite(pv_kwp) or pv_kwp < 0):
        raise ValueError("Supply finite, nonnegative generation/load and a nonempty profile")
    if (not isinstance(profiles.index, pd.DatetimeIndex) or profiles.index.tz is None
            or (len(profiles) > 1 and not np.all(np.diff(profiles.index.asi8) == 3600_000_000_000))):
        raise ValueError("Cell profiles require consecutive timezone-aware hourly timestamps")
    if mpc_config.horizon_hours < 1 or int(mpc_config.horizon_hours) != mpc_config.horizon_hours:
        raise ValueError("horizon_hours must be a positive integer")
    started = perf_counter()
    plant = DegradationPlant(c, cfg, pack_ah, profiles.ambient_c.iloc[0], radius_factor)
    # Hold a small scheduling reserve for side reactions and solver tolerances.
    scheduling_case = replace(c, soc_min=plant.soc_min + 0.001, soc_max=plant.soc_max - 0.001)
    if scheduling_case.soc_min >= scheduling_case.soc_max:
        raise ValueError("SOC range is too narrow for the scheduling reserve")
    pv = profiles.pv_kw_per_kwp.to_numpy() * pv_kwp
    load = profiles.load_kw.to_numpy()
    temp = profiles.ambient_c.to_numpy()
    frames, hourly = [], []
    for i in range(len(profiles)):
        charge = discharge = 0.
        seconds = 0
        while seconds < 3600:
            remaining_hours = (3600 - seconds) / 3600
            capacity = plant.capacity_proxy_kwh
            # Physical SOC may sit just outside the tighter scheduling reserve.
            # Clipping the estimate only affects the LP; the DFN state is untouched.
            estimated_energy = np.clip(plant.soc, scheduling_case.soc_min,
                                       scheduling_case.soc_max) * capacity
            charge_plan, discharge_plan = mpc_action(
                scheduling_case, mpc_config, pv, load, i, estimated_energy,
                capacity, c.max_c_rate * plant.nominal_kwh,
                first_step_hours=remaining_hours)
            request = (discharge_plan - charge_plan) / remaining_hours
            duration = min(cfg.control_period_s, 3600 - seconds)
            frame, dc_energy = plant.advance(request, temp[i], duration)
            frame["time"] = profiles.index[0] + pd.to_timedelta(frame.elapsed_s, unit="s")
            frames.append(frame)
            charge += max(0., -dc_energy)
            discharge += max(0., dc_energy)
            seconds += duration
        direct = min(pv[i], load[i])
        hourly.append({
            "soc": plant.soc, "capacity_proxy_kwh": plant.capacity_proxy_kwh,
            "pv_kw": pv[i], "load_kw": load[i], "charge_dc_kwh": charge,
            "discharge_dc_kwh": discharge,
            "unserved_kwh": max(0., load[i] - direct - discharge * c.discharge_efficiency),
            "curtailed_kwh": max(0., pv[i] - direct - charge / c.charge_efficiency),
        })
    trace = pd.concat(frames, ignore_index=True).set_index("time")
    dispatch = pd.DataFrame(hourly, index=profiles.index)
    last = trace.iloc[-1]
    summary = {
        "model": cfg.model, "parameter_set": cfg.parameter_set, "forecast": mpc_config.forecast,
        "radius_factor": radius_factor, "hours": len(profiles), "runtime_s": perf_counter() - started,
        "series_cells": plant.ns, "parallel_cells": plant.np,
        "requested_pack_ah": pack_ah, "actual_pack_ah": plant.np * plant.cell_ah,
        "actual_pack_nominal_kwh": plant.nominal_kwh,
        "voltage_min_v": trace.voltage_v.min(), "voltage_max_v": trace.voltage_v.max(),
        "voltage_guard_min_v": plant.v_min, "voltage_guard_max_v": plant.v_max,
        "soc_min": trace.soc.min(), "soc_max": trace.soc.max(),
        "soc_limit_min": plant.soc_min, "soc_limit_max": plant.soc_max,
        "predictions": plant.predictions, "rejected_predictions": plant.rejected_predictions,
        "limited_intervals": int(sum(bool(f.constraint_limited.iloc[0]) for f in frames)),
        "completed_hours": (trace.elapsed_s.iloc[-1] - trace.elapsed_s.iloc[0]) / 3600,
        "termination": "final time",
        **{key: float(last[key]) for key in VARIABLES if key.endswith(("_ah", "_pct"))},
        **{key: float(dispatch[key].sum()) for key in (
            "charge_dc_kwh", "discharge_dc_kwh", "unserved_kwh", "curtailed_kwh")},
    }
    return summary, trace, dispatch


def plot_cell_results(c, cfg, trace, folder, label, show=False):
    import matplotlib.pyplot as plt
    from pv_battery_model import save_figure

    fig, axes = plt.subplots(4, 1, figsize=(12, 11), sharex=True)
    t = trace.elapsed_s / 3600
    axes[0].plot(t, trace.voltage_v)
    param = pybamm.ParameterValues(cfg.parameter_set)
    for v in (param["Lower voltage cut-off [V]"] + cfg.voltage_margin_v,
              param["Upper voltage cut-off [V]"] - cfg.voltage_margin_v):
        axes[0].axhline(v, color="tab:red", ls="--", alpha=.6)
    axes[0].set_ylabel("Cell voltage [V]")
    axes[1].plot(t, trace.soc * 100, label="Cell SOC")
    axes[1].set_ylabel("SOC [%]")
    temperature = axes[1].twinx()
    temperature.plot(t, trace.temperature_k - 273.15, color="tab:orange", alpha=.5)
    temperature.set_ylabel("Cell temperature [C]")
    axes[2].plot(t, trace.requested_dc_kw, label="Requested DC power", alpha=.6)
    axes[2].plot(t, trace.pack_dc_kw, label="Applied DC power")
    axes[2].set_ylabel("Pack power [kW]")
    for key, name in (("sei_loss_ah", "SEI"), ("sei_cracks_loss_ah", "SEI on cracks"),
                      ("plating_loss_ah", "Lithium plating")):
        axes[3].plot(t, trace[key], label=name)
    axes[3].set(ylabel="Cell capacity loss [Ah]", xlabel="Elapsed time [h]")
    for ax in axes:
        ax.grid(alpha=.2)
        if ax.get_legend_handles_labels()[0]:
            ax.legend()
    fig.suptitle(f"{c.city} | {cfg.model} coupled O'Kane 2022 degradation | {label}")
    save_figure(fig, Path(folder), "cell_voltage_soc_degradation", show)


def run_cell_study(c, profiles, selected, cfg, mpc_config, output_dir, show=False):
    """Fresh, independent seasonal windows on selected hardware, including radii."""
    validate_config(c, cfg)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    windows = []
    first_year = profiles.index[0].year
    for label, date in cfg.windows:
        # Match the requested study year instead of hard-coding 2019.
        start = pd.Timestamp(date, tz="UTC").replace(year=first_year)
        end = start + pd.Timedelta(hours=cfg.hours)
        window = profiles.loc[(profiles.index >= start) & (profiles.index < end)]
        if len(window) != cfg.hours:
            raise ValueError(f"Incomplete {label} cell window starting {start}: expected {cfg.hours} hours")
        windows.append((label, window))
    if not windows:
        raise ValueError("Configure at least one cell-study window")
    metadata = {
        "case": asdict(c), "cell_config": asdict(cfg), "mpc": asdict(mpc_config),
        "pybamm_version": pybamm.__version__, "model_options": DEGRADATION_OPTIONS,
        "source": SOURCE, "paper_doi": "10.1039/D2CP00417H",
        "scope": "Independent fresh-cell windows; coarse sizing still uses empirical ageing.",
        "soc_definition": "Negative electrode bulk stoichiometry relative to beginning-of-life SOC endpoints.",
        "capacity_proxy": "Nominal energy times (1 - max(LLI, negative LAM, positive LAM)/100); not a capacity test.",
        "pack_assumption": "Identical balanced cells; Ns*Np*Icell*Vcell; explicit DC/AC efficiencies.",
        "controller": "Hourly energy LP with live state feedback and shrinking first interval; nonlinear DFN feasibility prediction before every current action.",
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, default=str), encoding="utf-8")
    summaries = []
    for label, window in windows:
        for factor in cfg.particle_radius_factors:
            print(f"{cfg.model} coupled degradation: {label}, radius x{factor:g}, {len(window)} h...", flush=True)
            summary, trace, dispatch = simulate_electrochemical(
                c, window, float(selected.pv_kwp), float(selected.initial_capacity_ah), cfg, mpc_config, factor)
            folder = output_dir / f"{label}_radius_{factor:g}"
            folder.mkdir(parents=True, exist_ok=True)
            trace.to_csv(folder / "cell_trace.csv", index_label="time")
            dispatch.to_csv(folder / "hourly_dispatch.csv", index_label="time")
            summaries.append({"window": label, "start": window.index[0], **summary})
            # Persist completed runs even if a subsequent variant fails.
            pd.DataFrame(summaries).to_csv(output_dir / "cell_summary.csv", index=False)
            plot_cell_results(c, cfg, trace, folder, f"{label}, radius x{factor:g}", show)
            print(f"  Completed in {summary['runtime_s']:.1f} s; "
                  f"cell voltage {summary['voltage_min_v']:.3f}-{summary['voltage_max_v']:.3f} V", flush=True)
    return pd.DataFrame(summaries)
