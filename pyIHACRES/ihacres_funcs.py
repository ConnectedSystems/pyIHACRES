from __future__ import division

from math import atan, exp, log, pi, tan
import numpy as np

def calc_cmd(prev_cmd, rainfall, et, effective_rainfall, recharge):
    """Calculate Catchment Moisture Deficit.

    Min value of CMD is 0.0 and is in represented in mm depth.
    A value of 0 indicates that the catchment is fully saturated.
    A value greater than 0 means that there is a moisture deficit.
    """
    cmd = prev_cmd + et + effective_rainfall + recharge - rainfall  # units in mm
    # cmd = interim_cmd - rainfall + et + effective_rainfall
    return max(0.0, cmd)
# End calc_cmd()


def calc_interim_cmd(cmd, param_d, rainfall):
    """Calculate interim CMD (M_{f}) in its linear form.

    Based on HydroMad implementation.

    :param cmd: float, current Catchment Moisture Deficit (M_{k})
    :param param_d: float, model parameter factor `d`
    :param rainfall: float, rainfall for current time step in mm

    :returns: float, interim CMD (M_{f})
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


def calc_trig_interim_cmd(cmd, param_d, rainfall):
    """Calculate interim CMD (M_{f}) in its trigonometric form.

    Based on HydroMad implementation.

    :param cmd: float, current Catchment Moisture Deficit (M_{k})
    :param param_d: float, model parameter factor `d`
    :param rainfall: float, rainfall for current time step in mm

    :returns: float, interim CMD (M_{f})
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


def calc_ft_interim(cmd, rain, d, d2, alpha):
    """Direct port of original Fortran implementation to calculate interim CMD.

    :param cmd: float, Catchment Moisture Deficit
    :param rain: float, rainfall for time step in mm
    :param d: float, flow threshold value
    :param d2: float, scaling factor applied to `d`
    :param alpha: float, took value from IHACRESparams.csv file for 406219 for dev purposes only

    :returns: tuple[float], interim CMD value, effective rainfall, recharge (all in mm)
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
        cmd = tmp_cmd - rain
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

                cmd = tmp_cmd / lam - epsilon * (1.0 - 1.0 / lam)
                e_rain = 0.0
            else:
                if (tmp_cmd > d1a):
                    tmp_rain = tmp_rain - depth_to_d

                tmp_cmd = d1a
                gamma = (alpha * d2 + (1.0 - alpha) * d1a) / (d1a * d2)
                cmd = tmp_cmd * exp(-tmp_rain * gamma)
                e_rain = alpha * (tmp_rain + 1.0 / d1a / gamma * (cmd - tmp_cmd))
            # End if
        else:
            gamma = (alpha * d2 + (1.0 - alpha) * d1a) / (d1a * d2)
            cmd = tmp_cmd * exp(-tmp_rain * gamma)
            e_rain = alpha * (tmp_rain + 1.0 / d1a / gamma * (cmd - tmp_cmd))
        # End if

        recharge = rain - (tmp_cmd - cmd) - e_rain
    # End if

    return cmd, e_rain, recharge
# End calc_ft_interim()


def calc_effective_rainfall(rainfall, cmd, Mf, d, d2, n=0.1):
    """
    :param rainfall: float, rainfall for time step
    :param cmd: float, previous CMD value
    :param Mf, float, interim CMD value
    :param d: float, threshold value
    :param d2: float, scaling factor applied to `d`
    :param n: float, scaling factor taken from Croke & Jakeman (2004).
              `n` = 0.1 in suitable for most cases
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


def calc_effective_rainfall_orig(rainfall, cmd, Mf, d, d2):
    d2 = d * d2
    if cmd > d2:
        e_rainfall = rainfall
    else:
        # effective rainfall = rainfall - prev_cmd + (cmd before ET loss is accounted for)
        e_rainfall = max(0.0, rainfall - cmd + Mf)
    # End if

    return e_rainfall
# End calc_effective_rainfall2()


def calc_ET_from_temp(e, T, interim_cmd, f, d):
    """Calculate evapotranspiration based on temperature data.

    Parameters `f` and `d` are used to calculate `g`, the value of the CMD
    which the ET rate will begin to decline due to insufficient
    water availability for plant transpiration.

    :param e: float, temperature to PET conversion factor (a stress threshold)
    :param T: float or None, temperature in degrees C
    :param interim_cmd: float, Catchment Moisture Deficit prior to accounting for ET losses (`M_{f}`)
    :param f: float, multiplication factor on `d`
    :param d: float, flow threshold factor
    """
    param_g = f * d
    et = e * T * exp(2.0 * (1.0 - (interim_cmd / param_g)))

    # temperature can be negative, so we have a min cap of 0.0
    return max(0.0, et)
# End calc_ET_from_temp()


def calc_ET(e, evap, interim_cmd, f, d):
    """Calculate evapotranspiration

    :param e: float, temperature to PET conversion factor (a stress threshold)
    :param evap: float, evaporation for given time step.
    :param interim_cmd: float, Catchment Moisture Deficit prior to accounting for ET losses (`M_{f}`)
    :param f: float, calibrated parameter that acts as a multiplication factor on `d`
    :param d: float, flow threshold factor
    """
    et = e * evap
    param_g = f * d
    et = et * min(1.0, exp(2.0 * (1.0 - (interim_cmd / param_g))))

    return et
# End calc_ET()


def calc_flow(tau, flow, v, e_rainfall):
    """Common function to calculate quick and slow flows.

    :param tau: float, `tau` value for flow
    :param flow: float, proportional quick or slow flow in ML/day
    :param v: float, `v` parameter
    :param e_rainfall: float, effective rainfall in mm (`E` parameter in literature)

    :returns: float, adjusted quick or slow flow in ML/day
    """
    # convert from 'tau' and 'v' to 'alpha' and 'beta'
    alpha = exp(-1.0 / tau)
    beta = v * (1.0 - alpha)

    flow = (-alpha * flow) + (1.0 - alpha) * (beta * e_rainfall)

    return flow
# End calc_flow()


def calc_flows(prev_flows, v_s, e_rainfall, taus):
    """
    Calculate quick and slow flow, and outflow.

    Calculates flows for current time step based on previous flows and current effective rainfall.

    :param prev_flows: tuple[float], previous quick and slow flow in ML/day
    :param v_s: float, proportional amount that goes to slow flow. v_s <= 1.0
    :param e_rainfall: float, current and previous effective rainfall
    :param taus: dict[q, s], time constant, quick and slow flow tau variables.
                               Represent the time required for the quickflow and slowflow
                               responses to fall to :math:`1/e` of their initial values after an impulse of rainfall

    :returns: tuple[float], quick, slow, outflow in ML/day
    """
    prev_quick, prev_slow = prev_flows
    v_q = 1.0 - v_s  # proportional quick flow
    quick = calc_flow(taus['q'], prev_quick, v_q, e_rainfall)
    slow = calc_flow(taus['s'], prev_slow, v_s, e_rainfall)
    outflow = (quick + slow)

    return quick, slow, outflow
# End calc_flows()


def routing(volume, storage_coef, inflow, flow, irrig_ext, gamma=0.0):
    """Linear routing used to convert effective rainfall into streamflow for a given time step.

    :param volume: float, catchment moisture deficit
    :param storage_coef: float, unknown parameter
    :param inflow: float, incoming streamflow (flow from previous node)
    :param flow: float, flow for the node (local flow)
    :param irrig_ext: float, volume of irrigation extraction in ML.
    :param gamma: float, unknown parameter which is always set to 0.0 in the original Fortran implementation.

    :returns: tuple[float], (cmd in mm, and streamflow in ML/day)
    """
    # print("Vol, inflow, l flow, gamma, irrig")
    # print(volume, inflow, flow, gamma, irrig_ext)
    threshold = volume + (inflow + flow + gamma) - irrig_ext
    if threshold > 0.0 and not np.isclose(threshold, 0.0):
        volume = 1.0 / (1.0 + storage_coef) * threshold
        outflow = storage_coef * volume
    else:
        volume = threshold
        outflow = 0.0
    # End if

    return volume, outflow
# End routing()


def calc_outflow(flow, extractions):
    """Calculate streamflow of node taking into account extractions

    :param flow: float, unmodified sum of quickflow and slowflow in ML/day
    :param extractions: float, water extractions that occurred in ML/day
    """
    return max(0.0, flow - extractions)
# End calc_outflow()


def calc_ft_flows(prev_quick, prev_slow, e_rain, recharge, area, a, b, loss=0.0):
    """Fortran port of flow calculation.

    :param prev_quick: float, previous quickflow storage
    :param prev_slow: float, previous slowflow storage
    :param e_rain: float, effective rainfall in mm
    :param recharge: float, recharge amount in mm
    :param area: float, catchment area in km^2
    :param a: float, `a` factor controlling quickflow rate
    :param b: float, `b` factor controlling slowflow rate
    :param loss: float, losses in mm depth

    :returns: tuple[float], quick store, slow store, outflow
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
