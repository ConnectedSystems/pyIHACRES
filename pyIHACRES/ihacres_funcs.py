from typing import Dict, Tuple, Optional
from math import atan, exp, log, pi, tan
import numpy as np


def calc_cmd(prev_cmd: float, rainfall: float, et: float, 
             effective_rainfall: float, recharge: float) -> float:
    """Calculate Catchment Moisture Deficit.

    Min value of CMD is 0.0 and is in represented in mm depth.
    A value of 0 indicates that the catchment is fully saturated.
    A value greater than 0 means that there is a moisture deficit.

    Returns
    -------
    cmd (in mm)
    """
    cmd = prev_cmd + et + effective_rainfall + recharge - rainfall
    return max(0.0, cmd)
# End calc_cmd()


def calc_interim_cmd(cmd: float, param_d: float, rainfall: float) -> float:
    """Calculate interim CMD (M_{f}) in its linear form.

    Parameters
    ----------
    cmd: current Catchment Moisture Deficit (M_{k})
    param_d: model parameter factor `d`
    rainfall: rainfall for current time step in mm

    Returns
    -------
    interim CMD (M_{f})
    """
    if cmd < param_d:
        Mf = cmd * exp(-rainfall / param_d)
    elif cmd < (param_d + rainfall):
        Mf = param_d * exp((-rainfall + cmd - param_d) / param_d)
    else:
        Mf = cmd - rainfall
    # End if

    return Mf
# End calc_interim_cmd()


def calc_trig_interim_cmd(cmd: float, param_d: float, rainfall: float) -> float:
    """Calculate interim CMD (M_{f}) in its trigonometric form.

    Parameters
    ----------
    cmd: current Catchment Moisture Deficit (M_{k})
    param_d: model parameter factor `d`
    rainfall: rainfall for current time step in mm

    Returns
    -------
    interim CMD (M_{f})
    """
    if cmd < param_d:
        Mf = 1.0 / tan((cmd / param_d) * (pi / 2.0))
        Mf = (2.0 * param_d / pi) * atan(1.0 / (pi * rainfall / (2.0 * param_d) + Mf))
    elif (rainfall < (param_d + rainfall)):
        Mf = (2.0 * param_d / pi) * atan(2.0 * param_d / (pi * (param_d - cmd + rainfall)))
    else:
        Mf = cmd - rainfall
    # End if

    return Mf
# End calc_trig_interim_cmd()


def calc_ft_interim(cmd: float, rain: float, 
                    d: float, d2: float, alpha: float) -> Tuple[float, float, float]:
    """Direct port of original Fortran implementation to calculate interim CMD.

    Calculates effective rainfall and recharge as a by-product.

    Parameters
    ----------
    cmd: Catchment Moisture Deficit
    rain: rainfall for time step in mm
    d: flow threshold value
    d2: scaling factor applied to `d`
    alpha: took value from IHACRESparams.csv file for 406219 for dev purposes only

    Returns
    -------
    (interim CMD value, effective rainfall, recharge)
    """
    d2 = d * d2

    tmp_cmd = cmd
    e_rain = 0.0
    recharge = 0.0
    if rain == 0.0:
        return cmd, e_rain, recharge
    # End if

    if tmp_cmd > (d2 + rain):
        # CMD never reaches d2, so all rain is effective
        Mf = tmp_cmd - rain
    else:
        if tmp_cmd > d2:
            tmp_rain = rain - (tmp_cmd - d2)  # leftover rain after reaching d2 threshold
            tmp_cmd = d2
        else:
            tmp_rain = rain
        # End if

        d1a = d * (2.0 - exp(-(rain / 50.0)**2))
        if tmp_cmd > d1a:
            eps = d2 / (1.0 - alpha)

            # original comment: now get rainfall to reach cmd = d1a
            # amount of rain necessary to get to threshold `d`
            depth_to_d = eps * log((alpha + tmp_cmd / eps) / (alpha + d1a / eps))
            if depth_to_d > tmp_rain:
                lam = exp(tmp_rain * (1.0 - alpha) / d2)
                epsilon = alpha * eps

                Mf = tmp_cmd / lam - epsilon * (1.0 - 1.0 / lam)
                e_rain = 0.0
            else:
                if (tmp_cmd > d1a):
                    tmp_rain = tmp_rain - depth_to_d

                tmp_cmd = d1a
                gamma = (alpha * d2 + (1.0 - alpha) * d1a) / (d1a * d2)
                Mf = tmp_cmd * exp(-tmp_rain * gamma)
                e_rain = alpha * (tmp_rain + 1.0 / d1a / gamma * (cmd - tmp_cmd))
            # End if
        else:
            gamma = (alpha * d2 + (1.0 - alpha) * d1a) / (d1a * d2)
            Mf = tmp_cmd * exp(-tmp_rain * gamma)
            e_rain = alpha * (tmp_rain + 1.0 / d1a / gamma * (cmd - tmp_cmd))
        # End if

        recharge = rain - (cmd - Mf) - e_rain
    # End if

    return cmd, e_rain, recharge
# End calc_ft_interim()


def calc_effective_rainfall(rainfall: float, cmd: float, 
                            d: float, d2: float, n:Optional[float] = 0.1) -> float:
    """
    Estimate effective rainfall.

    Parameters
    ---------
    rainfall: rainfall for time step
    cmd: previous CMD value
    Mf, interim CMD value
    d: threshold value
    d2: scaling factor applied to `d`
    n: scaling factor (default = 0.1)
       Default value is suitable for most cases (Croke & Jakeman, 2004)

    Returns
    -------
    estimate of effective rainfall

    References
    ----------
    .. [1] Croke, B.F.W., Jakeman, A.J. 2004
           A catchment moisture deficit module for the IHACRES rainfall-runoff model, 
           Environmental Modelling & Software, 19(1), pp. 1–5. 
           doi: 10.1016/j.envsoft.2003.09.001.
    """
    d2 = d * d2
    if cmd > d2:
        e_rainfall = rainfall
    else:
        f1 = min(1.0, cmd / d)
        f2 = min(1.0, cmd / d2)
        e_rainfall = rainfall * ((1.0 - n) * (1.0 - f1) + (n * (1.0 - f2)))
    # End if

    return e_rainfall
# End calc_effective_rainfall()


def calc_ET_from_temperature(e: float, T: float, interim_cmd: float, 
                             f: float, d: float) -> float:
    """Calculate evapotranspiration based on temperature data.

    Parameters `f` and `d` are used to calculate `g`, the value of the CMD
    which the ET rate will begin to decline due to insufficient
    water availability for plant transpiration.

    Parameters
    ----------
    e: temperature to PET conversion factor (a stress threshold)
    T: float or None, temperature in degrees C
    interim_cmd: Catchment Moisture Deficit prior to accounting for ET losses (`M_{f}`)
    f: multiplication factor on `d`
    d: flow threshold factor

    Returns
    -------
    estimate of evapotranspiration
    """
    param_g = f * d
    et = e * T * exp(2.0 * (1.0 - (interim_cmd / param_g)))

    # temperature can be negative, so we have a min cap of 0.0
    return max(0.0, et)
# End calc_ET_from_temperature()


def calc_ET(e: float, evap: float, interim_cmd: float, f: float, d: float) -> float:
    """Calculate evapotranspiration.

    Parameters
    ----------
    e: temperature to PET conversion factor (a stress threshold)
    evap: evaporation for given time step.
    interim_cmd: Catchment Moisture Deficit prior to accounting for ET losses (`M_{f}`)
    f: calibrated parameter that acts as a multiplication factor on `d`
    d: flow threshold factor

    Returns
    -------
    float, evapotranspiration
    """
    et = e * evap
    param_g = f * d
    et = et * min(1.0, exp(2.0 * (1.0 - (interim_cmd / param_g))))

    return et
# End calc_ET()


def calc_flow(tau: float, flow: float, v: float, e_rainfall: float) -> float:
    """Estimate quick and slow flow.

    Parameters
    ----------
    tau: `tau` value for flow
    flow: proportional quick or slow flow in ML/day
    v: `v` parameter
    e_rainfall: effective rainfall in mm (`E` parameter in literature)

    Returns
    -------
    float, adjusted quick or slow flow in ML/day
    """
    # convert from 'tau' and 'v' to 'alpha' and 'beta'
    alpha = exp(-1.0 / tau)
    beta = v * (1.0 - alpha)

    flow = (-alpha * flow) + (1.0 - alpha) * (beta * e_rainfall)

    return flow
# End calc_flow()


def calc_flows(prev_flows: Tuple[float], v_s: float, 
               e_rainfall: float, taus: Dict[str, float]) -> Tuple[float, float, float]:
    """
    Calculate quick/slow flow, and outflow.

    Estimates flows for current time step based on previous flows and current effective rainfall.

    prev_flows: previous quick and slow flow in ML/day
    v_s: proportional amount that goes to slow flow. v_s <= 1.0
    e_rainfall: current and previous effective rainfall
    taus: time constant, quick and slow flow tau variables ("q", "s").
            Represent the time required for the quickflow and slowflow
            responses to fall to $1/e$ of their initial values 
            after an impulse of rainfall

    Returns
    -------
    (quick, slow, outflow) in ML/day
    """
    prev_quick, prev_slow = prev_flows
    v_q = 1.0 - v_s  # proportional quick flow
    quick = calc_flow(taus['q'], prev_quick, v_q, e_rainfall)
    slow = calc_flow(taus['s'], prev_slow, v_s, e_rainfall)
    outflow = (quick + slow)

    return quick, slow, outflow
# End calc_flows()


def routing(volume: float, storage_coef: float, inflow: float, flow: float, 
            irrig_ext: float, gw_exchange: Optional[float]=0.0) -> Tuple[float, float]:
    """Linear routing used to convert effective rainfall into streamflow for a given time step.

    Parameters
    ----------
    volume: catchment moisture deficit
    storage_coef: storage factor
    inflow: incoming streamflow (flow from previous node)
    flow: flow for the node (local flow)
    irrig_ext: volume of irrigation extraction in ML
    gw_exchange: groundwater flux (default = 0.0)

    Returns
    -------
    (cmd in mm, and streamflow in ML/day)
    """
    threshold = volume + (inflow + flow + gw_exchange) - irrig_ext
    if threshold > 0.0 and not np.isclose(threshold, 0.0):
        volume = 1.0 / (1.0 + storage_coef) * threshold
        outflow = storage_coef * volume
    else:
        volume = threshold
        outflow = 0.0
    # End if

    return volume, outflow
# End routing()


def calc_outflow(flow: float, extractions: float) -> float:
    """Calculate streamflow of node taking into account extractions.

    Parameters
    ----------
    flow: unmodified sum of quickflow and slowflow in ML/day
    extractions: water extractions that occurred in ML/day

    Returns
    -------
    Outflow
    """
    return max(0.0, flow - extractions)
# End calc_outflow()


def calc_ft_flows(prev_quick: float, prev_slow: float, e_rain: float, 
                  recharge: float, area: float, 
                  a: float, b: float, loss=0.0) -> Tuple[float, float, float]:
    """Fortran port of flow calculation.

    prev_quick: previous quickflow storage
    prev_slow: previous slowflow storage
    e_rain: effective rainfall in mm
    recharge: recharge amount in mm
    area: catchment area in km^2
    a: `a` factor controlling quickflow rate
    b: `b` factor controlling slowflow rate
    loss: losses in mm depth

    Returns
    -------
    (quick store, slow store, outflow)
    """
    a2 = 0.5
    tmp_calc = prev_quick + (e_rain * area)
    if (tmp_calc - 0.5 * loss) > 0.0 \
            and not np.isclose(0.0, (tmp_calc - 0.5 * loss)):

        quick_store = 1.0 / (1.0 + a) * (tmp_calc - loss / 2.0)
        outflow = a * quick_store
    else:
        a2 = 0.0 if loss == 0.0 else max(0.0, min(1.0, (tmp_calc / loss)))
        quick_store = tmp_calc - a2 * loss
        outflow = 0.0
    # End if

    assert outflow >= 0.0, "Calculating quick store: Outflow cannot be negative"

    b2 = 1.0 - a2
    tmp_calc = prev_slow + (recharge * area)
    if ((tmp_calc - b2 * loss) > 0.0) \
            and not np.isclose(0.0, (tmp_calc - b2 * loss)):
        slow_store = 1.0 / (1.0 + b) * (tmp_calc - loss * b2)
        outflow = outflow + b * slow_store
    else:
        slow_store = tmp_calc - b2 * loss
    # End if

    assert outflow >= 0.0, "Calculating slow store: Outflow cannot be negative"

    return quick_store, slow_store, outflow
# End calc_ft_flows()
