from dataclasses import asdict, dataclass, replace
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from locations import city_slug, lookup_city, validate_location


@dataclass(frozen=True)
class Case:
    start: str = "2019-01-01"
    years: int = 5
    city: str = "Madrid"
    country_code: str | None = None  # optional ISO code, e.g. ES, NO, DE
    latitude_deg: float | None = None

    # snitt 2730 kWh i 2023 Spania
    annual_household_kwh: float = 2800.0 #kWh
    longitude_deg: float | None = None
    surface_tilt: float = 30.0
    surface_azimuth: float = 180.0  # pvlib: south=180
    pv_system_loss_pct: float = 14.0
    radiation_database: str = "auto"  # ERA5 north of SARAH3 coverage
    pvgis_url: str = "https://re.jrc.ec.europa.eu/api/v5_3/"
    local_timezone: str | None = None
    resolved_city: str = ""
    location_description: str = ""
    location_source: str = ""
    cache_dir: str = "pvgis_cache"
    refresh_pvgis: bool = False
    reference_pv_kwp: float = 6.0  # input plots only, not an optimized size

    # batteri
    battery_bus_v: float = 48.0
    initial_soc: float = 0.65
    soc_min: float = 0.15
    soc_max: float = 0.95
    max_c_rate: float = 0.5

    #effektivitet
    charge_efficiency: float = 0.95
    discharge_efficiency: float = 0.95

    # Degradation / tap
    # calendar fade = calendar_fade * sqrt(temperature-weighted years)
    # cycle fade = cycle_fade_per_efc * temperature-weighted full discharge cycles
    calendar_fade: float = 0.030 # fractional loss per sqrt(equivalent year)
    cycle_fade_per_efc: float = 0.000070 # 7% tap per 1000 discharge cycles
    reference_temp_c: float = 25.0
    temp_multiplier_per_10c: float = 2.0  # Q10-type acceleration

    # ønsket state of health ved lifecycle slutt
    min_soh: float = 0.80
    target_service_fraction: float = 0.995  # 99.5% av total last møtt
    # €900-1200 i Spania (2026) https://eu.jackery.com/blogs/buying-advice/how-much-do-solar-panels-cost-across-europe
    eur_per_kwp_pv: float = 1000.0
    # €400-500 (2025) https://en.u-energie.de/blogs/home-battery-storage-price-trends
    eur_per_kwh_battery: float = 400.0
    pv_options_kwp: tuple = (3.0, 4.0, 5.0, 6.0, 7.0, 8.0)
    battery_options_ah: tuple = (200, 250, 300, 350, 400, 500, 600, 750)
    random_seed: int = 67
    profile_csv: str | None = None  # columns: time, pv_kw_per_kwp, load_kw, ambient_c
    output_dir: str = "results"  # root; city and study years are appended

    @property
    def city_folder(self):
        suffix = f"_{city_slug(self.country_code)}" if self.country_code else ""
        return city_slug(self.city) + suffix

    @property
    def study_dir(self):
        hours = expected_hours(self)
        return Path(self.output_dir) / self.city_folder / f"{hours[0].year}_{hours[-1].year}"




@dataclass(frozen=True)
class MPCConfig:
    horizon_hours: int = 24
    forecast: str = "perfect"  # ideal benchmark; alternative: "persistence"
    unserved_weight: float = 1000.0
    throughput_weight: float = 0.001
    terminal_energy_value: float = 0.1



@dataclass(frozen=True)
class CellStudyConfig:
    enabled: bool = True
    model: str = "DFN"  # "SPMe" for faster exploratory runs
    parameter_set: str = "Chen2020"
    controller_label: str = "MPC perfect"
    # Fresh-cell, independent windows; no five-year electrochemical pre-ageing.
    windows: tuple = (("winter", "2019-01-15"), ("summer", "2019-07-15"))
    hours: int = 48
    particle_radius_factors: tuple = (1.0, 0.8, 1.2)
    series_cells: int = 13
    parallel_cells: int | None = None  # ceil(selected pack Ah / cell Ah)
    cell_nominal_voltage_v: float = 3.6  # pack nameplate estimate only
    thermal: str = "lumped"  # or "isothermal"
    cooling_w_m2_k: float = 10.0  # assumed, not calibrated pack cooling
    output_period_s: int = 60
    mesh_points: int = 15
    rtol: float = 1e-6
    atol: float = 1e-6
    output_dir: str = "pybamm_results"


# MAIN SETTINGS: change these before running subsequent cells.
# Case fields above control weather, load, sizing, SOC, efficiencies and ageing.
CASE = Case()
MPC = MPCConfig(horizon_hours=24, forecast="perfect")
CELL_STUDY = CellStudyConfig()
RUN_CONTROLLER_COMPARISON = False  # opt in: hourly MPC over five years is expensive


import hashlib
import json
from datetime import datetime, timezone
import pvlib


def resolve_case(c):
    """Resolve once at workflow entry; explicit coordinates allow offline studies."""
    if not c.city.strip():
        raise ValueError("City name cannot be empty")
    country = c.country_code.strip().upper() if c.country_code else None
    if country and (len(country) != 2 or not country.isascii() or not country.isalpha()):
        raise ValueError("country_code must be a two-letter ISO code, e.g. NO")
    c = replace(c, city=c.city.strip(), country_code=country)
    supplied = (c.latitude_deg, c.longitude_deg, c.local_timezone)
    if all(value is not None for value in supplied):
        validate_location(*supplied)
        return c
    if any(value is not None for value in supplied):
        raise ValueError("Provide latitude, longitude and timezone together, or let the city lookup set all three")
    return replace(c, **lookup_city(c.city, c.country_code, c.cache_dir))


def expected_hours(c):
    start = pd.Timestamp(c.start, tz="UTC")
    if (start != pd.Timestamp(year=start.year, month=1, day=1, tz="UTC")
            or c.years < 1 or int(c.years) != c.years):
        raise ValueError("Use a January 1 start and positive integer years")
    return pd.date_range(start, start + pd.DateOffset(years=c.years), freq="h", inclusive="left")


def validate_profiles(df, c):
    expected = expected_hours(c)
    required = ["pv_kw_per_kwp", "load_kw", "ambient_c"]
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Profiles are missing required columns: {', '.join(sorted(missing))}")
    if df.index.has_duplicates or not df.index.equals(expected):
        raise ValueError("Profiles must have every requested UTC hour exactly once, in order")
    values = df[required].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Profiles contain missing or non-finite values")
    if (df[["pv_kw_per_kwp", "load_kw"]] < 0).any().any() or df.load_kw.sum() <= 0:
        raise ValueError("PV/load must be nonnegative, with positive total demand")
    for col in ["poa_global", "wind_speed"]:
        if col in df and (not np.isfinite(df[col]).all() or (df[col] < 0).any()):
            raise ValueError(f"Invalid {col}")
    return df


def get_historical_pv(c):
    c = resolve_case(c)
    time = expected_hours(c)
    database = c.radiation_database
    if database == "auto":
        database = ("PVGIS-ERA5" if c.latitude_deg > 65 or abs(c.longitude_deg) > 65
                    else "PVGIS-SARAH3")
    request = dict(latitude=c.latitude_deg, longitude=c.longitude_deg,
                   start=time[0].year, end=time[-1].year,
                   raddatabase=database, surface_tilt=c.surface_tilt,
                   surface_azimuth=c.surface_azimuth, components=True,
                   pvcalculation=True, peakpower=1.0, pvtechchoice="crystSi",
                   mountingplace="free", loss=c.pv_system_loss_pct,
                   usehorizon=True, trackingtype=0, map_variables=True,
                   url=c.pvgis_url, outputformat="json")
    key = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()[:16]
    folder = Path(c.cache_dir) / c.city_folder / key
    raw_file, meta_file = folder / "hourly_raw.csv", folder / "metadata.json"
    if raw_file.exists() and meta_file.exists() and not c.refresh_pvgis:
        raw = pd.read_csv(raw_file, index_col=0)
        raw.index = pd.to_datetime(raw.index, utc=True)
        provenance = json.loads(meta_file.read_text())
        if provenance["request"] != request:
            raise ValueError("Cache metadata does not match request")
    else:
        try:
            raw, metadata = pvlib.iotools.get_pvgis_hourly(**request, timeout=60)
        except Exception as exc:
            raise RuntimeError(f"PVGIS download failed for {c.city} ({database}): {exc}. "
                               "Check internet/API availability, coverage and year range, "
                               "or use an existing matching cache. No synthetic fallback is used.") from exc
        provenance = {"city": c.city, "country_code": c.country_code,
                      "local_timezone": c.local_timezone,
                      "request": request, "pvlib_version": pvlib.__version__,
                      "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                      "pvgis_metadata": metadata}
        folder.mkdir(parents=True, exist_ok=True)
        raw.to_csv(raw_file, index_label="time")
        meta_file.write_text(json.dumps(provenance, indent=2, default=str))
    raw.index = pd.to_datetime(raw.index, utc=True)
    if not raw.index.is_monotonic_increasing or raw.index.has_duplicates:
        raise ValueError("PVGIS timestamps are duplicated or out of order")
    hours = raw.index.floor("h")
    if not hours.equals(time):
        raise ValueError("PVGIS response has missing/extra hours; refusing to fill gaps")
    fields = ["P", "temp_air", "poa_direct", "poa_sky_diffuse", "poa_ground_diffuse", "wind_speed"]
    if not np.isfinite(raw[fields].to_numpy(dtype=float)).all():
        raise ValueError("Non-finite PVGIS observations")
    df = pd.DataFrame(index=hours)
    df.index.name = "time"
    df["pv_kw_per_kwp"] = raw["P"].to_numpy() / 1000.0  # W -> kW, for 1 kWp
    df["ambient_c"] = raw["temp_air"].to_numpy()
    df["poa_global"] = raw[["poa_direct", "poa_sky_diffuse", "poa_ground_diffuse"]].sum(axis=1).to_numpy()
    df["wind_speed"] = raw["wind_speed"].to_numpy()
    df["source_time_utc"] = raw.index.astype(str)
    df.attrs["provenance"] = provenance
    return df


def synthetic_load(time, c):
    c = resolve_case(c)
    local = time.tz_convert(c.local_timezone)
    hour, doy = local.hour.to_numpy(), local.dayofyear.to_numpy()
    codes, dates = pd.factorize(local.date, sort=True)
    rng = np.random.default_rng(c.random_seed)
    variability = np.clip(rng.normal(1.0, 0.12, len(dates)), 0.70, 1.30)
    shape = (0.65 + 0.65 * np.exp(-0.5 * ((hour - 8) / 1.8)**2)
             + np.exp(-0.5 * ((hour - 20) / 2.4)**2))
    load = pd.Series(shape * (1 + 0.15 * np.cos(2*np.pi*(doy-15)/365.25))
                     * variability[codes], index=time)
    return load.groupby(time.year).transform(lambda x: x * c.annual_household_kwh / x.sum())


def make_profiles(c: Case) -> pd.DataFrame:
    c = resolve_case(c)
    if c.profile_csv is not None:
        df = pd.read_csv(c.profile_csv, index_col="time")
        df.index = pd.to_datetime(df.index, utc=True)  # naive CSV times explicitly treated as UTC
        df.attrs["provenance"] = {"source": str(c.profile_csv), "note": "User-supplied profiles"}
    else:
        df = get_historical_pv(c)
        df["load_kw"] = synthetic_load(df.index, c)
    return validate_profiles(df, c)


import matplotlib.dates as mdates


def save_figure(fig, output_dir, name, show=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    for extension in ("png", "pdf"):
        fig.savefig(output_dir / f"{name}.{extension}", dpi=180, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    return output_dir / f"{name}.png"


def plot_inputs(c, profiles, output_dir, show=False):
    c = resolve_case(c)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    load_label = "Supplied load" if c.profile_csv else "Synthetic load"
    energy = pd.DataFrame({"PV": profiles.pv_kw_per_kwp * c.reference_pv_kwp,
                           "Load": profiles.load_kw})  # hourly power * 1 h = kWh
    annual = energy.resample("YS").sum()
    annual["specific_pv_yield_kwh_per_kwp"] = profiles.pv_kw_per_kwp.resample("YS").sum()
    annual["hours"] = profiles.load_kw.resample("YS").count()
    annual.to_csv(output_dir / "annual_inputs.csv")
    monthly = energy.resample("MS").sum()
    monthly.to_csv(output_dir / "monthly_inputs.csv")
    def save(fig, name):
        return save_figure(fig, output_dir, name, show)
    fig, axes = plt.subplots(2, 1, figsize=(12, 7))
    positions = np.arange(len(annual))
    axes[0].bar(positions - .2, annual.PV, width=.4, color="#d79419", label="PV")
    axes[0].bar(positions + .2, annual.Load, width=.4, color="#285c87", label="Load")
    axes[0].set_xticks(positions, annual.index.year)
    axes[0].legend()
    axes[0].set(title="Annual energy input", ylabel="Energy [kWh/year]", xlabel="Year (UTC)")
    axes[1].plot(monthly.index, monthly.PV, color="#d79419", label="PV")
    axes[1].plot(monthly.index, monthly.Load, color="#285c87", label="Load")
    axes[1].legend()
    axes[1].set(title="Monthly energy — seasonal surplus and deficit", ylabel="Energy [kWh/month]", xlabel="Date (UTC)")
    for ax in axes: ax.grid(axis="y", alpha=.25)
    fig.suptitle(f"{c.city} | {c.reference_pv_kwp:g} kWp reference PV (not optimized) | {load_label.lower()}")
    save(fig, "input_energy_overview")

    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    if "poa_global" in profiles:
        daily_solar = profiles.poa_global.resample("D").sum() / 1000
        axes[0].plot(daily_solar.index, daily_solar, color="#d79419", lw=.65)
        axes[0].set_ylabel("Daily POA irradiation\n[kWh/m²/day]")
    else:
        axes[0].text(.5, .5, "Irradiance not supplied in CSV", transform=axes[0].transAxes, ha="center")
    temp = profiles.ambient_c.resample("D").agg(["min", "mean", "max"])
    axes[1].fill_between(temp.index, temp["min"], temp["max"], alpha=.2, color="#285c87", label="Daily range")
    axes[1].plot(temp.index, temp["mean"], lw=.7, label="Daily mean", color="#285c87")
    axes[1].set(ylabel="Ambient temperature [°C]", xlabel="Date (UTC)")
    axes[1].legend()
    for ax in axes: ax.grid(alpha=.2)
    fig.suptitle(f"{c.city} | Weather inputs | irradiation uses hourly sample × 1 h")
    save(fig, "input_weather_overview")

    fig, axes = plt.subplots(2, 2, figsize=(13, 6), sharex="col")
    local = profiles.tz_convert(c.local_timezone)
    year = pd.Timestamp(c.start).year
    for col, (month, season) in enumerate([(1, "Winter"), (7, "Summer")]):
        start = pd.Timestamp(year=year, month=month, day=15, tz=c.local_timezone)
        end_ns = start.value + 7 * 24 * 3600 * 1_000_000_000
        week = local.loc[(local.index.asi8 >= start.value) & (local.index.asi8 < end_ns)]
        axes[0,col].plot(week.index, week.pv_kw_per_kwp*c.reference_pv_kwp, color="#d79419", label="Modelled PV")
        axes[0,col].plot(week.index, week.load_kw, color="#285c87", label=load_label)
        axes[0,col].set(title=f"{season}: {start:%d %b %Y} — seven illustrative days", ylabel="Power [kW]")
        axes[0,col].legend()
        axes[1,col].plot(week.index, week.ambient_c, color="#a84646")
        axes[1,col].set(ylabel="Ambient temperature [°C]", xlabel=f"Local date ({c.local_timezone})")
        axes[1,col].xaxis.set_major_formatter(mdates.DateFormatter("%d %b", tz=local.index.tz))
        for ax in axes[:,col]: ax.grid(alpha=.2)
    fig.suptitle(f"{c.city} | Hourly inputs | {c.reference_pv_kwp:g} kWp reference PV | example weeks")
    save(fig, "input_winter_summer_weeks")
    return annual


from scipy.optimize import linprog
from functools import lru_cache


@lru_cache(maxsize=32)
def energy_matrix(h):
    # x = [DC charge, DC discharge, end-of-hour stored energy], in kWh.
    a = np.zeros((h, 3*h))
    for k in range(h):
        a[k, k], a[k, h+k], a[k, 2*h+k] = -1, 1, 1
        if k: a[k, 2*h+k-1] = -1
    return a


def mpc_action(c, cfg, pv, load, i, energy, capacity, power_limit):
    h = min(cfg.horizon_hours, len(load)-i)
    if cfg.forecast == "perfect":
        generation, demand = pv[i:i+h], load[i:i+h]
    elif cfg.forecast == "persistence":
        # Repeat the latest available UTC day, with measured current-hour inputs.
        # First day: zero future PV and current demand (no future data leakage).
        inds = i + np.arange(h)
        past = inds - 24 * (np.floor_divide(inds-i, 24)+1)
        generation = np.where(past >= 0, pv[np.maximum(0, past)], 0.0)
        demand = np.where(past >= 0, load[np.maximum(0, past)], load[i])
        generation[0], demand[0] = pv[i], load[i]
    else:
        raise ValueError("forecast must be perfect or persistence")
    surplus = np.maximum(generation-demand, 0)
    deficit = np.maximum(demand-generation, 0)
    # Direct PV serves load first. These disjoint bounds prohibit simultaneous
    # charge/discharge exactly without binary variables or artificial energy dumping.
    bounds = ([(0, min(s*c.charge_efficiency, power_limit)) for s in surplus]
              + [(0, min(d/c.discharge_efficiency, power_limit)) for d in deficit]
              + [(c.soc_min*capacity, c.soc_max*capacity)]*h)
    objective = np.r_[np.full(h, cfg.throughput_weight),
                      np.full(h, cfg.throughput_weight-cfg.unserved_weight*c.discharge_efficiency),
                      np.zeros(h)]
    objective[-1] = -cfg.terminal_energy_value
    rhs = np.zeros(h); rhs[0] = energy
    result = linprog(objective, A_eq=energy_matrix(h), b_eq=rhs,
                     bounds=bounds, method="highs")
    if not result.success:
        raise RuntimeError(f"MPC failed at hour {i}: {result.message}")
    return max(0., result.x[0]), max(0., result.x[h])





def simulate(c: Case, profiles: pd.DataFrame, pv_kwp: float,
             initial_capacity_ah: float, keep_trace: bool = False,
             controller: str = "rule_based", mpc_config: MPCConfig = MPC):
    if controller not in ("rule_based", "mpc"):
        raise ValueError("Unknown controller")
    if mpc_config.horizon_hours < 1 or int(mpc_config.horizon_hours) != mpc_config.horizon_hours:
        raise ValueError("horizon_hours must be a positive integer")
    if not (0 < c.soc_min <= c.initial_soc <= c.soc_max <= 1):
        raise ValueError("Invalid SOC limits/initial state")
    if not (0 < c.charge_efficiency <= 1 and 0 < c.discharge_efficiency <= 1):
        raise ValueError("Invalid efficiency")
    if initial_capacity_ah <= 0 or pv_kwp < 0 or c.max_c_rate <= 0:
        raise ValueError("Invalid sizing/power limit")
    if len(profiles) == 0 or not np.isfinite(profiles[["pv_kw_per_kwp", "load_kw", "ambient_c"]]).all().all():
        raise ValueError("Empty/nonfinite profiles")
    if (profiles[["pv_kw_per_kwp", "load_kw"]] < 0).any().any() or profiles.load_kw.sum() <= 0:
        raise ValueError("Invalid generation/demand")
    if len(profiles) > 1 and not np.all(np.diff(profiles.index.asi8) == 3600000000000):
        raise ValueError("Expected consecutive hourly records")

    pv_unit = profiles["pv_kw_per_kwp"].to_numpy()
    load = profiles["load_kw"].to_numpy()
    temp = profiles["ambient_c"].to_numpy()
    n = len(profiles)
    pv = pv_kwp * pv_unit
    details = []
    fade_spill_total = 0.0
    nominal_kwh = initial_capacity_ah * c.battery_bus_v / 1000.0
    battery_kwh = c.initial_soc * nominal_kwh
    max_dc_transfer_kwh = c.max_c_rate * nominal_kwh  # hourly time step

    temp_years = 0.0
    temp_efc = 0.0
    discharged_kwh = charged_kwh = unserved_kwh = curtailed_kwh = 0.0
    desired_capacity_ah_seen = initial_capacity_ah
    if keep_trace:
        soc_trace = np.empty(n)
        cap_trace = np.empty(n)

    for i in range(n):
        # nåværende kapasitet med tap
        soh = max(0.0, 1.0 - c.calendar_fade * np.sqrt(temp_years)
                  - c.cycle_fade_per_efc * temp_efc)
        actual_capacity_kwh = nominal_kwh * soh
        desired_capacity_ah_seen = min(desired_capacity_ah_seen, initial_capacity_ah * soh)
        energy_before = battery_kwh
        fade_spill = max(0., battery_kwh-c.soc_max*actual_capacity_kwh)
        fade_spill_total += fade_spill
        battery_kwh -= fade_spill
        if actual_capacity_kwh <= 0:
            raise RuntimeError("Battery exhausted: replacement model required")

        pv_kwh = pv_kwp * pv_unit[i]  # power kW * one-hour step
        load_kwh = load[i]
        direct_kwh = min(pv_kwh, load_kwh)
        surplus_kwh = pv_kwh - direct_kwh
        deficit_kwh = load_kwh - direct_kwh

        if controller == "mpc":
            charge_kwh, discharge_kwh = mpc_action(
                c, mpc_config, pv, load, i, battery_kwh,
                actual_capacity_kwh, max_dc_transfer_kwh)
        else:
            charge_kwh = min(surplus_kwh*c.charge_efficiency,
                             max(0., c.soc_max*actual_capacity_kwh-battery_kwh), max_dc_transfer_kwh)
            discharge_kwh = min(deficit_kwh/c.discharge_efficiency,
                                max(0., battery_kwh-c.soc_min*actual_capacity_kwh), max_dc_transfer_kwh)
        # Apply first action only; conservative clipping handles solver roundoff.
        charge_kwh = min(charge_kwh, surplus_kwh*c.charge_efficiency,
                         max(0., c.soc_max*actual_capacity_kwh-battery_kwh), max_dc_transfer_kwh)
        discharge_kwh = min(discharge_kwh, deficit_kwh/c.discharge_efficiency,
                            max(0., battery_kwh-c.soc_min*actual_capacity_kwh), max_dc_transfer_kwh)
        battery_kwh += charge_kwh-discharge_kwh
        unserved = max(0., deficit_kwh-discharge_kwh*c.discharge_efficiency)
        curtailed = max(0., surplus_kwh-charge_kwh/c.charge_efficiency)
        charged_kwh += charge_kwh
        discharged_kwh += discharge_kwh
        unserved_kwh += unserved
        curtailed_kwh += curtailed
        if keep_trace:
            details.append((energy_before, battery_kwh, actual_capacity_kwh,
                            charge_kwh, discharge_kwh, unserved, curtailed, fade_spill))

        # temperature-weighted calendar age and equivalent full discharge cycles
        multiplier = c.temp_multiplier_per_10c ** (
            (temp[i] - c.reference_temp_c) / 10.0
        )
        temp_years += multiplier / (365.25 * 24.0)
        temp_efc += multiplier * discharge_kwh / nominal_kwh

        if keep_trace:
            soc_trace[i] = battery_kwh / actual_capacity_kwh
            cap_trace[i] = initial_capacity_ah * soh

    soh_end = max(0.0, 1.0 - c.calendar_fade * np.sqrt(temp_years)
                  - c.cycle_fade_per_efc * temp_efc)
    final_capacity_ah = initial_capacity_ah * soh_end
    desired_capacity_ah_seen = min(desired_capacity_ah_seen, final_capacity_ah)
    total_load = float(load.sum())
    service_fraction = 1.0 - unserved_kwh / total_load
    cost_eur = (pv_kwp * c.eur_per_kwp_pv
                + nominal_kwh * c.eur_per_kwh_battery)
    row = {
        "controller": controller,
        "forecast": mpc_config.forecast if controller == "mpc" else "none",
        "fade_spill_kwh": fade_spill_total,
        "final_energy_kwh": battery_kwh,
        "pv_kwp": float(pv_kwp),
        "initial_capacity_ah": float(initial_capacity_ah),
        "initial_energy_kwh": nominal_kwh,
        "end_capacity_ah": final_capacity_ah,
        "desired_capacity_ah": desired_capacity_ah_seen,
        "end_soh_pct": 100.0 * soh_end,
        "served_pct": 100.0 * service_fraction,
        "unserved_kwh": unserved_kwh,
        "curtailed_kwh": curtailed_kwh,
        "charge_throughput_kwh": charged_kwh,
        "discharge_throughput_kwh": discharged_kwh,
        "example_cost_eur": cost_eur,
        "feasible": (
            desired_capacity_ah_seen >= c.min_soh * initial_capacity_ah
            and service_fraction >= c.target_service_fraction
        ),
    }
    if not keep_trace:
        return row
    trace = pd.DataFrame({"soc": soc_trace,
                          "capacity_ah": cap_trace}, index=profiles.index)
    trace[["energy_before_kwh", "energy_kwh", "available_capacity_kwh",
           "charge_dc_kwh", "discharge_dc_kwh", "unserved_kwh",
           "curtailed_kwh", "fade_spill_kwh"]] = np.asarray(details)
    trace["pv_kw"] = pv
    trace["load_kw"] = load
    # SOC and capacity both refer to the dispatch hour; exact final SOH is in row.
    return row, trace


def size_system(c: Case, profiles: pd.DataFrame, controller="rule_based", mpc_config=MPC,
                require_feasible=True):
    for options, minimum, name in [(c.pv_options_kwp, 0, "PV"), (c.battery_options_ah, 0, "battery")]:
        if (not options or len(set(options)) != len(options)
                or not np.isfinite(options).all() or min(options) < minimum
                or (name == "battery" and min(options) == 0)):
            raise ValueError(f"Invalid or duplicate {name} sizing options")
    trials = [simulate(c, profiles, p, q, controller=controller, mpc_config=mpc_config) for p in c.pv_options_kwp
              for q in c.battery_options_ah]
    table = pd.DataFrame(trials).sort_values("example_cost_eur")
    feasible = table.loc[table["feasible"]]
    # om ingen av de definerte oppsettene oppfyller krav
    if feasible.empty:
        if require_feasible:
            raise RuntimeError("No configuration satisfies both requirements; expand the PV/battery sizing grid")
        # A diagnostic candidate, explicitly not a feasible/optimal design.
        best = table.sort_values(["served_pct", "end_soh_pct", "example_cost_eur"],
                                 ascending=[False, False, True]).iloc[0]
        return table, best
    best = feasible.iloc[0]
    return table, best


def plot_results(c: Case, trace: pd.DataFrame, best: pd.Series, output_dir: Path,
                 show=False, label=None):
    daily = trace.resample("D").agg({
        "capacity_ah": "last", "soc": ["min", "mean", "max"]
    })
    cap = daily["capacity_ah", "last"]
    soc_min = daily["soc", "min"] * 100
    soc_mean = daily["soc", "mean"] * 100
    soc_max = daily["soc", "max"] * 100

    required_capacity_ah = (
        c.min_soh * best["initial_capacity_ah"]
    )

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    ax1.plot(cap.index, cap, label="Remaining battery capacity")
    ax1.axhline(
        required_capacity_ah,
        ls="--",
        color="tab:red",
        label=f"{c.min_soh:.0%} SOH ({required_capacity_ah:.0f} Ah)"
    )
    ax1.set_ylabel("Battery capacity [Ah]")
    ax1.grid(alpha=0.25)
    ax1.legend()

    ax2.fill_between(soc_mean.index, soc_min.to_numpy(), soc_max.to_numpy(),
                     alpha=0.3, label="Daily SOC range")
    ax2.plot(soc_mean.index, soc_mean, lw=0.85, label="Daily mean SOC")
    ax2.axhline(100 * c.soc_min, ls="--", color="tab:red", alpha=0.7)
    ax2.axhline(100 * c.soc_max, ls="--", color="tab:red", alpha=0.7)
    ax2.set_ylabel("State of charge [%]")
    ax2.set_xlabel("Date (UTC)")
    ax2.set_ylim(0, 100)
    ax2.grid(alpha=0.25)
    ax2.legend()
    status = "Feasible design" if best["feasible"] else "Diagnostic design: requirements not met"
    label = label or str(best["controller"])
    fig.suptitle(
        f"{c.city} | {c.years}-year off-grid case | {label}\n"
        f"{status} | PV {best['pv_kwp']:.1f} kWp | "
        f"battery {best['initial_capacity_ah']:.0f} Ah / "
        f"{best['initial_energy_kwh']:.1f} kWh"
    )
    return save_figure(fig, output_dir, "capacity_and_soc", show)


def plot_sizing(c, table, selected, output_dir, show=False):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ax, (metric, title, cmap) in zip(axes, [
        ("served_pct", "Demand served [%]", "viridis"),
        ("end_soh_pct", "Final state of health [%]", "viridis"),
        ("example_cost_eur", "Illustrative capital cost [EUR]", "cividis"),
    ]):
        grid = table.pivot(index="initial_capacity_ah", columns="pv_kwp", values=metric).sort_index()
        image = ax.imshow(grid.to_numpy(), origin="lower", aspect="auto", cmap=cmap)
        ax.set_xticks(range(len(grid.columns)), [f"{p:g}" for p in grid.columns])
        ax.set_yticks(range(len(grid.index)), [f"{q:g}" for q in grid.index])
        ax.set(xlabel="PV size [kWp]", ylabel="Battery capacity [Ah]", title=title)
        for _, row in table.loc[table.feasible].iterrows():
            ax.plot(grid.columns.get_loc(row.pv_kwp), grid.index.get_loc(row.initial_capacity_ah),
                    "o", ms=9, mfc="none", mec="white", label="Meets both targets")
        ax.plot(grid.columns.get_loc(selected.pv_kwp), grid.index.get_loc(selected.initial_capacity_ah),
                "*", ms=14, color="tomato", markeredgecolor="black")
        fig.colorbar(image, ax=ax, shrink=.8)
    status = "Selected feasible design" if selected.feasible else "No feasible design; star is a diagnostic candidate"
    fig.suptitle(f"{c.city} | {selected.controller} sizing | {status}\n"
                 f"Targets: demand served ≥ {c.target_service_fraction:.1%}, SOH ≥ {c.min_soh:.0%} "
                 "| white circles meet both targets")
    return save_figure(fig, output_dir, "sizing_overview", show)


def plot_dispatch(c, trace, output_dir, label, show=False):
    monthly = trace[["pv_kw", "load_kw", "unserved_kwh", "curtailed_kwh",
                     "charge_dc_kwh", "discharge_dc_kwh"]].resample("MS").sum()
    monthly = monthly.rename(columns={"pv_kw": "pv_kwh", "load_kw": "load_kwh"})
    monthly.to_csv(Path(output_dir) / "monthly_energy.csv", index_label="time")
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    axes[0].plot(monthly.index, monthly.pv_kwh, label="PV generation", color="#d79419")
    axes[0].plot(monthly.index, monthly.load_kwh, label="Demand", color="#285c87")
    axes[0].set_ylabel("Energy [kWh/month]")
    axes[1].plot(monthly.index, monthly.unserved_kwh, label="Unserved demand", color="#a84646")
    axes[1].plot(monthly.index, monthly.curtailed_kwh, label="Curtailed PV", color="#d79419")
    axes[1].set(ylabel="Energy [kWh/month]", xlabel="Date (UTC)")
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend()
    fig.suptitle(f"{c.city} | {label} | Monthly generation, demand and energy shortfall")
    save_figure(fig, output_dir, "monthly_energy_balance", show)

    local = trace.tz_convert(c.local_timezone)
    fig, axes = plt.subplots(3, 2, figsize=(14, 9), sharex="col")
    for col, (month, season) in enumerate([(1, "Winter"), (7, "Summer")]):
        start = pd.Timestamp(year=pd.Timestamp(c.start).year, month=month, day=15, tz=c.local_timezone)
        end_ns = start.value + 7 * 24 * 3600 * 1_000_000_000
        week = local.loc[(local.index.asi8 >= start.value) & (local.index.asi8 < end_ns)]
        axes[0, col].plot(week.index, week.pv_kw, label="PV", color="#d79419")
        axes[0, col].plot(week.index, week.load_kw, label="Demand", color="#285c87")
        axes[0, col].set(title=f"{season}: seven days from {start:%d %b %Y}", ylabel="Power [kW]")
        axes[1, col].plot(week.index, week.charge_dc_kwh, label="Charge (DC)")
        axes[1, col].plot(week.index, -week.discharge_dc_kwh, label="Discharge (DC)")
        axes[1, col].plot(week.index, week.unserved_kwh, label="Unserved", color="#a84646")
        axes[1, col].set_ylabel("Energy [kWh/hour]")
        axes[2, col].plot(week.index, 100 * week.soc, label="SOC", color="#285c87")
        axes[2, col].axhline(100 * c.soc_min, ls="--", color="tab:red", alpha=.6)
        axes[2, col].axhline(100 * c.soc_max, ls="--", color="tab:red", alpha=.6)
        axes[2, col].set(ylabel="SOC [%]", ylim=(0, 100), xlabel=f"Local date ({c.local_timezone})")
        axes[2, col].xaxis.set_major_formatter(mdates.DateFormatter("%d %b", tz=local.index.tz))
        for ax in axes[:, col]:
            ax.grid(alpha=.2)
            ax.legend(fontsize="small")
    fig.suptitle(f"{c.city} | {label} | Hourly dispatch in illustrative winter and summer weeks")
    save_figure(fig, output_dir, "dispatch_winter_summer_weeks", show)


def plot_controller_comparison(c, table, output_dir, show=False):
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    for ax, metric, title, target in zip(
        axes, ["served_pct", "end_soh_pct", "curtailed_kwh"],
        ["Demand served [%]", "Final state of health [%]", "Total curtailed PV [kWh]"],
        [100 * c.target_service_fraction, 100 * c.min_soh, None],
    ):
        ax.bar(table.label, table[metric], color=["#285c87", "#d79419", "#a84646"][:len(table)])
        ax.set_title(title)
        ax.tick_params(axis="x", labelrotation=20)
        ax.grid(axis="y", alpha=.2)
        if target is not None:
            ax.axhline(target, ls="--", color="tab:red", label="Target")
            ax.set_ylim(0, 105)
            ax.legend()
    fig.suptitle(f"{c.city} | Controller comparison using identical hardware (no MPC resizing)")
    return save_figure(fig, output_dir, "controller_comparison", show)


def run_study(c=CASE, controller="rule_based", mpc_config=MPC,
              compare_controllers=RUN_CONTROLLER_COMPARISON, show=False):
    """Save reproducible inputs, sizing trials and dispatch results for one city."""
    expected_hours(c)
    if not np.isfinite(c.annual_household_kwh) or c.annual_household_kwh <= 0:
        raise ValueError("Annual load must be positive and finite")
    c = resolve_case(c)
    root = c.study_dir
    inputs_dir = root / "inputs"
    result_name = controller if controller == "rule_based" else f"mpc_{mpc_config.forecast}"
    results_dir = root / "results" / result_name
    for folder in (inputs_dir, results_dir):
        folder.mkdir(parents=True, exist_ok=True)
    print(f"Location: {c.location_description or c.city} | "
          f"{c.latitude_deg:.5f}, {c.longitude_deg:.5f} | {c.local_timezone}", flush=True)
    profiles = make_profiles(c)
    profiles.to_csv(inputs_dir / "hourly_profiles.csv", index_label="time")
    metadata = {"city": c.city, "case": asdict(c), "mpc": asdict(mpc_config),
                "sizing_controller": controller, "compare_controllers": compare_controllers,
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "provenance": profiles.attrs.get("provenance", {}),
                "assumptions": ["Hourly off-grid energy model; no grid trading.",
                                "Load and economic assumptions are user inputs, not city statistics.",
                                "Battery ageing and costs are illustrative."]}
    # Each controller run keeps its own configuration; shared inputs live above it.
    (results_dir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    plot_inputs(c, profiles, inputs_dir, show)
    print(f"Sizing {len(c.pv_options_kwp) * len(c.battery_options_ah)} designs using {controller}...", flush=True)
    table, selected = size_system(c, profiles, controller, mpc_config, require_feasible=False)
    table.to_csv(results_dir / "sizing_trials.csv", index=False)
    pd.DataFrame([selected]).to_csv(results_dir / "selected_design.csv", index=False)
    plot_sizing(c, table, selected, results_dir, show)

    jobs = [(result_name, controller, mpc_config)]
    if compare_controllers:
        jobs = [("rule_based", "rule_based", mpc_config),
                ("mpc_perfect", "mpc", replace(mpc_config, forecast="perfect")),
                ("mpc_persistence", "mpc", replace(mpc_config, forecast="persistence"))]
    summaries = []
    for label, strategy, cfg in jobs:
        print(f"Simulating {label} for {c.city}...", flush=True)
        summary, trace = simulate(c, profiles, float(selected.pv_kwp),
                                  float(selected.initial_capacity_ah), True, strategy, cfg)
        folder = results_dir / "dispatch" / label
        folder.mkdir(parents=True, exist_ok=True)
        trace.to_csv(folder / "hourly_dispatch.csv", index_label="time")
        summaries.append({"label": label, **summary})
        plot_results(c, trace, pd.Series(summary), folder, show, label)
        plot_dispatch(c, trace, folder, label, show)
    summary_table = pd.DataFrame(summaries)
    summary_table.to_csv(results_dir / "controller_summary.csv", index=False)
    if compare_controllers:
        plot_controller_comparison(c, summary_table, results_dir, show)
    status = "feasible" if selected.feasible else "infeasible: expand the sizing grid or revise the assumptions"
    print(f"{c.city}: {selected.pv_kwp:g} kWp PV, {selected.initial_capacity_ah:g} Ah battery "
          f"({status}).\nOutputs: {root.resolve()}", flush=True)
    return {"case": c, "output_dir": root, "sizing": table,
            "selected": selected, "controllers": summary_table}


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description="Size and visualize a PV-battery system for a European city.")
    parser.add_argument("--city", default=CASE.city, help="City name, optionally 'City, Region'")
    parser.add_argument("--country", help="Two-letter ISO country code, e.g. NO, ES or DE")
    parser.add_argument("--start", default=CASE.start, help="January 1 of the first historical year")
    parser.add_argument("--years", type=int, default=CASE.years)
    parser.add_argument("--annual-load-kwh", type=float, default=CASE.annual_household_kwh)
    parser.add_argument("--pv-options", nargs="+", type=float, default=CASE.pv_options_kwp, metavar="KWP")
    parser.add_argument("--battery-options", nargs="+", type=float, default=CASE.battery_options_ah, metavar="AH")
    parser.add_argument("--output-dir", default=CASE.output_dir, help="Root for city/year output folders")
    parser.add_argument("--cache-dir", default=CASE.cache_dir)
    parser.add_argument("--profile-csv", help="Use supplied UTC hourly profiles instead of downloading PVGIS")
    parser.add_argument("--latitude", type=float, help="Offline location override; requires longitude and timezone")
    parser.add_argument("--longitude", type=float)
    parser.add_argument("--timezone", help="Local IANA timezone, e.g. Europe/Oslo")
    parser.add_argument("--radiation-database", default=CASE.radiation_database,
                        choices=["auto", "PVGIS-SARAH3", "PVGIS-ERA5"])
    parser.add_argument("--refresh-pvgis", action="store_true")
    parser.add_argument("--controller", choices=["rule_based", "mpc"], default="rule_based")
    parser.add_argument("--forecast", choices=["perfect", "persistence"], default=MPC.forecast)
    parser.add_argument("--horizon-hours", type=int, default=MPC.horizon_hours)
    parser.add_argument("--compare-controllers", action="store_true",
                        help="Run rule-based, perfect MPC and persistence MPC on the selected hardware (slow)")
    parser.add_argument("--show", action="store_true", help="Also open figures interactively")
    args = parser.parse_args(argv)
    c = replace(CASE, city=args.city.strip(), country_code=args.country,
                latitude_deg=args.latitude, longitude_deg=args.longitude, local_timezone=args.timezone,
                start=args.start, years=args.years, annual_household_kwh=args.annual_load_kwh,
                pv_options_kwp=tuple(args.pv_options), battery_options_ah=tuple(args.battery_options),
                output_dir=args.output_dir, cache_dir=args.cache_dir, profile_csv=args.profile_csv,
                radiation_database=args.radiation_database, refresh_pvgis=args.refresh_pvgis)
    cfg = replace(MPC, horizon_hours=args.horizon_hours, forecast=args.forecast)
    try:
        expected_hours(c)
        if c.annual_household_kwh <= 0:
            raise ValueError("Annual load must be positive")
        run_study(c, args.controller, cfg, args.compare_controllers, args.show)
    except (ValueError, RuntimeError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
