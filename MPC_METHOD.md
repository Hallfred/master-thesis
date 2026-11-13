## Hourly MPC: scope and source adaptation

Sandeep Machani (2024), *Modelling battery operations for household energy management*, §§2.3, 3.1–3.4 supplies the multiperiod energy balance, hourly rolling horizon, state feedback structure. His §3.4 also discusses forecast-error settlement through the grid, which is not reproduced in this off-grid model. His grid-price objective, LSTM and DoD ageing formulation are **not reproduced** here.

Jiefeng **Hu** et al. (2018), *A model predictive control strategy of PV-Battery microgrid under variable power generations and load conditions*, Applied Energy 221, 195–203, DOI: [10.1016/j.apenergy.2018.03.085](https://doi.org/10.1016/j.apenergy.2018.03.085), §§3–5: SOC limits, balancing, load shedding and PV curtailment inform this EMS. Its §§4 MPVC/MPPC are converter-level switching control, **not this hourly scheduling algorithm**.

For a 24-hour horizon, choose battery-side charge/discharge energies $c_k,d_k$ and stored energy $E_{k+1}$:

$$E_{k+1}=E_k+c_k-d_k,$$
$$PV_k+d_k\eta_d+u_k=L_k+c_k/\eta_c+w_k.$$

Here $u$ is unserved energy, $w$ is curtailed PV. Direct PV serves demand first; charge only from net surplus, discharge only into net deficit. This restriction prevents simultaneous charging/discharging and intentionally excludes grid trading and flexible demand. All powers are integrated over **one hour**. SOC limits use current faded capacity; the C-rate power rating uses initial capacity, preserving the baseline assumption.

Minimize $1000\sum u_k+0.001\sum(c_k+d_k)-0.1E_H$. These are **chosen numerical priorities, not sourced prices or calibrated ageing costs**. Curtailment has no penalty, so wasting energy through cycling is not rewarded. Terminal stored-energy value limits end-of-horizon depletion; it is a heuristic, not a long-term feasibility guarantee. Test horizon/weight sensitivity before thesis conclusions. The controller freezes capacity over each horizon, updates actual ageing hourly, and re-solves. The existing hypothetical ageing model is unchanged; this is not degradation-optimized MPC.

Only the first action is applied. Current-hour PV/load are assumed measured at dispatch; future values use either `perfect` (ideal-forecast benchmark) or `persistence` (repeat the latest available UTC day; first-day zero future PV). This does not test subhourly forecast errors. The horizon shrinks at the dataset end; terminal reward remains active. Solver failure stops with an explicit error, never a silent fallback.

The five-year comparison uses the **same baseline-selected hardware** for all controllers. The existing sizing sweep remains rule-based by default; `size_system(..., controller="mpc")` enables a new MPC sizing study, which is computationally more expensive. No MPC sizing optimum is claimed from a single-design comparison.

Solver: SciPy `linprog(method="highs")`, https://docs.scipy.org/doc/scipy/reference/optimize.linprog-highs.html .

Sandeep uses a price-availability-dependent prediction horizon (35 to 11 hours after the initial phase). The fixed 24-hour horizon here is a new choice for this off-grid case, not a replication of his schedule. Current-hour inputs are treated as known, so the persistence comparison tests future forecast quality, not ex-post correction of the current-hour dispatch.
