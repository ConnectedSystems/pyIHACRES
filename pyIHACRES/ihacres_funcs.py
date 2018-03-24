from __future__ import division

import numpy as np


def calc_cmd(prev_cmd, rainfall, et, effective_rainfall):
    """Calculate Catchment Moisture Deficit.

    Min value of CMD is 0.0 and is in represented in mm depth.
    A value of 0 indicates that the catchment is fully saturated.
    A value greater than 0 means that there is a moisture deficit.
    """
    cmd = prev_cmd - rainfall + et + effective_rainfall  # units in mm
    return max(0.0, cmd)
# End calc_cmd()


def calc_effective_rainfall(cmd, d):
    e_rain = 1 - min(1, cmd / d)
    return e_rain
# End calc_effective_rainfall()


def calc_ET_from_temp(e, T, p_cmd, f, d):
    """Calculate evapotranspiration based on temperature data.

    Parameters `f` and `d` are used to calculate `g`, the value of the CMD
    which the ET rate will begin to decline due to insufficient
    water availability for plant transpiration.

    :param e: float, temperature to PET conversion factor (a stress threshold)
    :param T: float or None, temperature in degrees C
    :param p_cmd: float, Catchment Moisture Deficit prior to accounting for ET losses
    :param f: float, multiplication factor on `d`
    :param d: float, flow threshold factor
    """
    g = f * d
    et = e * T * np.exp(2.0 * (1 - (p_cmd / g)))
    return et
# End calc_ET()


def calc_ET(e, et, p_cmd, f, d):
    """Stub function. Available data provides ET, but currently unsure if some scaling should be applied.
    """
    return et


def partition_flow(alpha, flow, v, e_rainfall):
    # q_k = (a[0] * q[k - 1] + a[1] * q[k - 2]) + b[0] * u[k] + b[1] * u[k - 1]
    return (alpha * flow) + (1 - alpha) * (v * e_rainfall)


def calc_flows(timestep, flows, e_rainfall, ab):
    """Calculate total streamflow.

    Total streamflow is the sum of quickflow and slowflow.

    :param flows: list[float], time series of outflow
    :param e_rainfall: float, effective rainfall in mm
    :param ab: tuple[float] of length 2, calibrated constants for quick and slow flows,
                   values between 0 and 1.
    # :param p_factors: tuple[float] of length 2, proportion of effective rainfall diverted into quickflow and slowflow,
    #                   Constitutes parameters `v_q` and `v_s` (must sum to 1.0)

    :returns: tuple[float] of length 2, quickflow and slowflow
    """
    # a, b = tau_v_to_ab(taus, p_factors)
    a, b = ab

    k = timestep
    try:
        q_k1 = flows[k - 1]
        q_k2 = flows[k - 2]
    except IndexError:
        q_k1 = 0.0  # flows[-1]
        q_k2 = 0.0
    # End if

    try:
        e_k = e_rainfall[k]
        e_k1 = e_rainfall[k - 1]
    except IndexError:
        e_k = e_rainfall[-1]
        e_k1 = 0.0
    # End if

    quickflow = max(0.0, a[0] * q_k1)
    slowflow = max(0.0, a[1] * q_k2)
    # quickflow = partition_flow(quick_alpha, prev_quickflow, quick_v, e_rainfall)
    # slowflow = partition_flow(slow_alpha, prev_slowflow, slow_v, e_rainfall)

    outflow = (quickflow + slowflow) + b[0] * e_k + b[1] * e_k1

    return quickflow, slowflow, outflow
# End calc_flows()


def calc_outflow(flow, extractions):
    """Calculate streamflow of node taking into account extractions

    :param flow: float, unmodified sum of quickflow and slowflow in ML/day
    :param extractions: float, water extractions that occurred in ML/day
    """
    # outflow = flow + b[0] * e_rainfall[k] + b[1] * e_rainfall[k - 1]
    return flow - extractions
# End calc_outflow()


def tau_v_to_ab(tau, v):
    """Calculate tau and v to `a` and `b` parameters.

    Subject to:
    * tau_q < tau_s
    * v_q + v_s = 1.0

    :param tau: tuple[float]
    :param v: tuple[float], are the proportions of effective rainfall diverted to quickflow and slowflow
    """
    assert tau[0] < tau[1], "`tau` quickflow parameter must be less than slowflow `tau` parameter"
    assert np.sum(v) == 1.0, "`v` parameter must sum to 1.0"

    alpha = [np.exp(-1.0 / tau_i) for tau_i in tau]
    beta = [v[i] * (1.0 - alpha[i]) for i in range(len(v))]

    if len(alpha) == 1:
        a = [alpha[0]]
        b = [beta[0]]
    elif len(alpha) == 2:
        a = [(alpha[0] + alpha[1]),
             -1.0 * (alpha[0] * alpha[1])]

        b = [beta[0] + beta[1],
             -1.0 * (beta[0] * (alpha[1]) +
                     beta[1] * (alpha[0]))]
    else:
        raise ValueError("Length of `alpha` must be 1 or 2")
    # End if

    return a, b
# End tau_v_to_ab()
