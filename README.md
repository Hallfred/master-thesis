# Master's thesis preparatory work

A PV-battery microgrid model for Smarter Solar: Adaptive Batteries.

Exploring how lithium-ion cells and battery packs can be designed and optimized specifically for use in solar energy systems, with particular attention to how real world environmental conditions shape their performance and long term behavior. 
Develop and use computational models that capture the influence of temperature variations, irradiation patterns, charging and discharging behavior, and geographic differences, to evaluate how these factors affect cell efficiency, stability, and aging. 
By combining detailed solar input data with battery simulations, investigate trends, limitations, and opportunities for more robust and adaptive storage solutions.

## Run a European city study

Install the dependencies and run the script (Python 3.10 or newer):

```sh
python -m pip install -r requirements.txt
python pv_battery_model.py --city Oslo --country NO
python pv_battery_model.py --city "Tromsø" --country NO --years 1
python pv_battery_model.py --city Berlin --country DE --annual-load-kwh 3500
```

`--city` accepts a city name rather than a fixed list. Open-Meteo/GeoNames supplies
coordinates and the local timezone. A country code is optional; for ambiguous names,
use `--country` or a qualifier such as `--city "Bergen, Vestland"`. The most populous
exact European name match is preferred, with the most populous search match used
for aliases. The resolved place, coordinates and timezone are printed before the
study and recorded in its metadata. Check that they describe your intended city.
The configured city name is retained in every figure title and its output path,
including names with accents. European IANA timezones, European Atlantic islands
and Cyprus are accepted; non-European search matches are excluded.

The defaults retain the original **2019–2023** historical period, 2,800 kWh/year
synthetic household load and illustrative sizing/economic assumptions. Changing
city changes the weather and load timezone; it does **not** supply local household
consumption statistics, electricity prices or calibrated battery parameters.
`--years 1` is useful for a quicker exploratory run. Use `--help` for all options.

PVGIS downloads and city lookups are cached in `pvgis_cache/`. Automatic radiation
database selection uses SARAH3 within its conservative coverage bounds and ERA5
above 65°N or outside ±65° longitude. `--radiation-database PVGIS-ERA5` overrides
this choice. Availability still depends on PVGIS coverage and historical years;
unsupported requests fail explicitly without fabricated weather. Use supplied
profiles for locations/years the service cannot provide.

For an offline study, supply a CSV and all three location overrides:

```sh
python pv_battery_model.py --city Oslo --country NO --years 1 --profile-csv my_profiles.csv --latitude 59.91273 --longitude 10.74609 --timezone Europe/Oslo
```

The CSV must have `time,pv_kw_per_kwp,load_kw,ambient_c`, with every requested UTC
hour exactly once. Optional `poa_global` [W/m²] and `wind_speed` [m/s] preserve
additional weather inputs. Naive timestamps are explicitly interpreted as UTC.
An offline run with these overrides makes no city or weather API calls.

## Saved inputs and results

For Oslo's default five-year study:

```text
results/
  oslo_no/
    2019_2023/
      inputs/
        hourly_profiles.csv
        annual_inputs.csv
        monthly_inputs.csv
        input_energy_overview.{png,pdf}
        input_weather_overview.{png,pdf}
        input_winter_summer_weeks.{png,pdf}
      results/
        rule_based/                 # or mpc_perfect / mpc_persistence sizing
          run_metadata.json         # resolved location, settings, data provenance
          sizing_trials.csv         # every tested PV/battery design
          selected_design.csv
          sizing_overview.{png,pdf}
          controller_summary.csv
          dispatch/
            rule_based/             # additional controllers when requested
              hourly_dispatch.csv
              monthly_energy.csv
              capacity_and_soc.{png,pdf}
              monthly_energy_balance.{png,pdf}
              dispatch_winter_summer_weeks.{png,pdf}
```

`--output-dir` sets the root; the city, country (when known), and study years are
always appended. Repeat runs with the same city, years and sizing controller
replace matching files. Use a different root to retain distinct parameter studies.
Figures are saved without opening windows; add `--show` for interactive display.
Energy aggregates use UTC; illustrative winter/summer weeks use the city's local
timezone, including daylight saving time.

The sizing plot shows service, final SOH, illustrative capital cost, feasible
designs and the selected candidate. If no candidate meets both targets, all trials
and diagnostic plots are still saved. The candidate is explicitly marked
**infeasible**, rather than reported as an optimum. Northern off-grid sites may
need a much larger grid of sizes or revised assumptions. For example:

```sh
python pv_battery_model.py --city Oslo --country NO --pv-options 6 10 15 20 --battery-options 500 1000 2000
```

Python use follows the same workflow:

```python
from pv_battery_model import Case, run_study

study = run_study(Case(city="Vienna", country_code="AT", years=1))
print(study["output_dir"])
print(study["selected"])
```

For lower-level calls, use `c = resolve_case(Case(city="Vienna", country_code="AT"))`
before calling `make_profiles`, `size_system` and the plotting functions. Importing
the module does not download data or start simulations.

## Hourly MPC

Read `MPC_METHOD.md` for equations, source attribution, and limitations.
`pv_battery_model.py` exposes the same definitions without notebook execution.
Run `python -m unittest discover -v` for offline coupled-DFN integration tests,
voltage/SOC constraints, degradation, and controller energy-balance checks.

MPC tuning: `MPCConfig(horizon_hours=24, forecast="perfect")`.
Rule-based sizing is the default. For a fresh MPC sizing sweep:
`python pv_battery_model.py --city Madrid --country ES --controller mpc`.
This is slower than the rule-based sweep and is not performed automatically.
Add `--forecast persistence` or change `--horizon-hours` to tune MPC.

`--compare-controllers` runs rule-based, perfect-forecast MPC and persistence MPC
on the same selected hardware, saves each dispatch trace, and adds
`controller_comparison.{png,pdf}`. This can take a long time over five years.
A comparison on fixed hardware is not evidence of an MPC-optimal system size.
The sizing sweep uses illustrative empirical ageing. Selected hardware also runs
the detailed electrochemical studies below. No converter switching or reactive
power controller from Hu et al. is implemented.

## Coupled DFN degradation and voltage protection

The default workflow now runs PyBaMM's **DFN with OKane2022 parameters** after
sizing. `dfn_degradation.py` enables solvent-diffusion-limited SEI growth, SEI on
particle cracks, partially reversible lithium plating, porosity changes,
negative-electrode swelling/cracking, positive-electrode swelling, and
stress-driven loss of active material. See the
[PyBaMM coupled-degradation example](https://docs.pybamm.org/en/v25.12.2/source/examples/notebooks/models/coupled-degradation.html)
and [O'Kane et al. (2022)](https://doi.org/10.1039/D2CP00417H).

By default, independent fresh-cell winter and summer windows each cover 48 hours
in the first study year. Each runs particle-radius factors 1.0, 0.8 and 1.2,
applied to both electrodes. These are mechanism studies on selected hardware;
they do not replace the multi-year empirical sizing calculation or reproduce
five years of DFN ageing. The low-level `simulate_electrochemical` function also
accepts longer consecutive hourly profiles without resetting cell state.

The cell controller runs the energy MPC using live electrode SOC and a capacity
proxy based on LLI/LAM. Every 600 seconds it updates the estimate, solves the LP
with the remaining current-hour duration, and predicts a full DFN trajectory for
the proposed constant current. Predictions must stay within **2.55-4.15 V**,
50 mV inside the unchanged physical 2.5/4.2 V cutoffs. Continuous voltage and SOC
events supplement 60-second output checks. A bounded current search reduces
infeasible actions before committing them. Nonzero actions also need a feasible
one-hour rest prediction at the measured ambient temperature: voltage can rise
after charging as plated lithium re-intercalates. An additional 5 mV current-action
reserve leaves room for relaxation. Rejected trajectories consume no
simulated time and add no degradation to the accepted state. Excess PV is
curtailed and unmet demand is recorded, preserving the AC energy balance using
actual integrated cell voltage/current. Physical cutoff events remain enabled;
an infeasible rest or unusable cell raises an explicit error.

The default SOC upper limit is **99.5%** (`--soc-max 0.995`); the voltage guard
often limits charging earlier. MPC holds a small additional SOC reserve for
side reactions. There are no long constant-voltage charging holds. SOC is the
negative-electrode bulk stoichiometry relative to beginning-of-life endpoints,
and the scheduling capacity proxy is **not** a measured capacity/SOH test.

The pack uses 13 series cells and `ceil(selected pack Ah / 5 Ah)` parallel cells,
with power calculated as `Ns * Np * Vcell * Icell`. Actual rounded pack capacity
and nameplate energy are reported, including their difference from the sizing
model's 48 V assumption. Lumped thermal dynamics use measured hourly ambient
temperature, an assumed cooling coefficient, and persistent cell temperature.
The cracking mesh uses 30 radial points and 5 points in each electrode/separator;
mesh and tolerance convergence remain necessary before thesis conclusions.

```sh
# Full default seasonal/radius study after sizing
python pv_battery_model.py --city Madrid --country ES --years 1
# Faster exploratory run with the same coupled mechanisms
python pv_battery_model.py --city Madrid --country ES --years 1 --cell-model SPMe --cell-hours 24 --particle-radius-factors 1.0
# Empirical sizing/dispatch only
python pv_battery_model.py --city Madrid --country ES --years 1 --skip-cell-study
```

`CellStudyConfig` exposes control/output periods, rest prediction duration, thermal
assumptions, mesh, solver tolerances, voltage reserves and pack layout. For Python-only coarse runs,
pass `cell_config=CellStudyConfig(enabled=False)` to `run_study`.

Results are saved under the sizing controller's `pybamm_results/` directory:
`run_metadata.json`, `cell_summary.csv`, and one folder per season/radius with
`cell_trace.csv`, `hourly_dispatch.csv` and voltage/SOC/degradation PNG/PDF figures.
The trace includes SEI, crack-SEI and plating capacity losses, LLI, electrode LAM,
porosity, temperature, current, requested/applied pack power and limiting flags.
The summary records runtime, completed duration, voltage extrema and rejected
predictions. Short-window losses should not be extrapolated into lifetime SOH.

Sources:

- [Open-Meteo geocoding / GeoNames](https://open-meteo.com/en/docs/geocoding-api)
- [pvlib hourly PVGIS interface](https://pvlib-python.readthedocs.io/en/stable/reference/generated/pvlib.iotools.get_pvgis_hourly.html)
- [PVGIS API](https://joint-research-centre.ec.europa.eu/photovoltaic-geographical-information-system-pvgis/using-pvgis-5/api-non-interactive-service_en)
- [PVGIS 5.3 datasets and coverage](https://publications.jrc.ec.europa.eu/repository/handle/JRC139355)
