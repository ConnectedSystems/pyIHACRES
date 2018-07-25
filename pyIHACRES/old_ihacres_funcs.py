"""Collection of functions adapted from mjasher
https://github.com/mjasher/aus-hydro-data
"""

from __future__ import division

import numpy as np


def alphabeta_to_vtau(alpha, beta):

    v = [beta[i] / (1.0 + alpha[i]) for i in range(len(alpha))]
    tau = [-1.0 / np.log(-1.0 * a) for a in alpha]

    return v, tau
# End alphabeta_to_vtau()


def ab_to_alphabeta(a, b):

    # TODO check for complex or no poles
    beta = np.roots(a)

    # TODO generalize for 3 or more stores
    size = 2
    M = np.empty((size, size))
    M[0, :] = 1.0
    M[1, :] = beta
    alpha = np.linalg.solve(M, b)

    return alpha, beta[::-1]
# End ab_to_alphabeta()


def arma(a, b, q, u, k):
    """arma calculation.
    Returns flow for current time step.

    :param a:
    :param b:
    :param q:
    :param u:
    :param k: int, current time step
    """
    q_k = (a[0] * q[k - 1] + a[1] * q[k - 2]) + b[0] * u[k] + b[1] * u[k - 1]
    return q_k
# End arma()


def ihacres_cmd_mod(prev_cmd, rainfall, evap, area, scaled_ext, d1, d, e, f, a, b):
    """Catchment Moisture Deficit based on initial work by mjasher.

    Dependent Publications:
        A catchment moisture deficit module for the IHACRES rainfall-runoff model
        and
        Corrigendum to the above

    :param prev_cmd: float, the catchment moisture deficit (CMD) for the previous time step
    :param rainfall: float, rainfall depth (mm)
    :param evap: float, evaporation in mm
    :param area: float, subcatchment area in unknown unit
    :param scaled_ext: float, scaled water extractions
    :param d: float, flow threshold
    :param e: float, calibrated constant to calculate Potential Evapotranspiration (PET)
    :param f: float, multiplication factor on `d`
    :param a: float, node parameter a, currently unknown
    :param b: float, node parameter b, currently unknown

    :returns: tuple[float], effective rainfall, storage (CMD), PET
    """

    if prev_cmd < d:
        cmd = prev_cmd * np.exp(-1.0 * rainfall / d)
    elif cmd < d + rainfall:
        cmd = d * np.exp(-1.0 * (rainfall - prev_cmd + d) / d)
    else:
        cmd = prev_cmd - rainfall
    # End if

    f = f * d

    evap = e * evap
    if cmd > f:
        evap = evap * np.exp(2.0 * (1.0 - cmd / f))
    # End if

    cmd = cmd + evap
    effective_rainfall = rainfall + cmd - prev_cmd

    return effective_rainfall, cmd

    d1_a = d1 * (2.0 - np.exp(-(rainfall / 50.0)**2))

    # TODO: REPLACE WITH INPUT PARAMETER
    d2 = 0.0
    alpha = 0.0
    prev_S1flow = 0.0  # Previous quick flow?
    prev_S2flow = 0.0  # Previous slow flow?
    loss = 0.0  # Not sure what this is, is always 0.0 in original Fortran code

    cmd = prev_cmd
    if rainfall > 0.0:
        if prev_cmd > (d2 + rainfall):
            cmd = prev_cmd - rainfall
        else:
            if prev_cmd > d2:
                mod_rain = rainfall - (prev_cmd - d2)
                tmp_cmd = d2
            else:
                tmp_cmd = prev_cmd
                mod_rain = rainfall
            # End if

            if tmp_cmd > d1_a:
                epsilon = d2 / (1.0 - alpha)  # epislonn in original fortran
                tmp_rain = epsilon * np.log((alpha + (tmp_cmd / epsilon)) / (alpha + (d1_a / epsilon)))

                if tmp_rain >= mod_rain:
                    lmbda = np.exp(mod_rain * (1.0 - alpha) / d2)
                    ep = alpha * epsilon
                    cmd = tmp_cmd / lmbda - ep * (1.0 - (1.0 / lmbda))
                else:
                    if tmp_cmd > d1_a:
                        mod_rain = mod_rain - tmp_rain
                    # End if

                    mod_cmd = d1_a
                    gamma = (alpha * d2 + (1 - alpha) * d1_a) / (d1_a * d2)
                    cmd = mod_cmd * np.exp(-mod_rain * gamma)
                    U_k = alpha * (mod_rain + 1.0 / d1_a / gamma * (cmd - mod_cmd))
                # End if
            else:
                gamma = (alpha * d2 + (1 - alpha) * d1_a) / (d1_a * d2)
                cmd = mod_cmd * np.exp(-mod_rain * gamma)
                U_k = alpha * (mod_rain + 1.0 / d1_a / gamma * (cmd - mod_cmd))
            # End if

            r = rainfall - (mod_cmd - cmd) - U_k
        # End if
    # End if

    et = e * evap
    if cmd > f:
        et = et * np.exp((1.0 - cmd / f) * 2.0)
    # End if

    tmp_cmd = prev_cmd + et + U_k + r - rainfall
    a2 = 0.5
    if (prev_S2flow + U_k * area - 0.5 * loss) > 0.0:
        tmp_flow1 = 1 / (1 + a) * (prev_S2flow + U_k * area - loss / 2.0)
        flow = a * tmp_flow1
    else:
        if loss < 0.0:
            a2 = 0.0
        else:
            a2 = (prev_S2flow + U_k * area) / loss  # division by zero?
            a2 = min(1.0, max(0.0, a2))
        # End if

        tmp_flow1 = prev_S1flow + U_k * area - a2 * loss
        flow = 0.0
    # End if

    b2 = 1.0 - a2
    if ((prev_S2flow + r * area - b2 * loss) < 0):
        tmp_flow2 = 1 / (1 + b) * (prev_S2flow + r * area - loss * b2)
        flow = flow + b * tmp_flow2
    else:
        tmp_flow2 = prev_S2flow + r * area - b2 * loss
    # End if

    return quickflow, slowflow, flow, mod_cmd, U_k, r, crd

# End ihacres_cmd_mod()


def ihacres_cmd(M_k_m_1, P_k, E_k, d, e, f):
    """Catchment Moisture Deficit based on initial work by mjasher.

    Publications:
        A catchment moisture deficit module for the IHACRES rainfall-runoff model
        and
        Corrigendum to the above

    :param M_k_m_1: float, the catchment moisture deficit (CMD)
    :param P_k: float, rainfall depth
    :param E_k: float, evaporation
    :param d: float, flow threshold
    :param e: float, calibrated constant to calculate Potential Evapotranspiration (PET)
    :param f: float, multiplication factor on `d`

    :returns: tuple[float], effective rainfall, storage (CMD), PET
    """
    if M_k_m_1 < d:
        M_k = M_k_m_1 * np.exp(-1.0 * P_k / d)
    elif M_k_m_1 < d + P_k:
        M_k = d * np.exp(-1.0 * (P_k - M_k_m_1 + d) / d)
    else:
        M_k = M_k_m_1 - P_k
    # End if

    f = f * d

    if M_k > f:
        E_k = e * E_k * np.exp(2.0 * (1.0 - M_k / f))
    else:
        E_k = e * E_k
    # End if

    M_k = M_k + E_k
    U_k = P_k + M_k - M_k_m_1

    return U_k, M_k, E_k

# End ihacres_cmd()


def cmd(M_k_m_1, P_k, T_k, d, e, f):
    """
    Publications:
        A catchment moisture deficit module for the IHACRES rainfall-runoff model
        and
        Corrigendum to the above

    :param M_k: float, the catchment moisture deficit (CMD)
    :param P_k: float, rainfall depth
    :param T_k: float, temperature, not needed if we're passing in ET
    :param d: float, flow threshold
    :param e: float, calibrated constant to calculate Potential Evapotranspiration (PET)
    :param f: float, represents a function value between 0 and 1 with the following properties:

              Monotonically increasing function if M_k is between 0 and :math:`d`

              if :math:`M_k` is zero then the catchment is fully saturated and :math:`f` is 0.

              :math:`f` is 1 if :math:`M_k > 0`


    Original notes:
    ET[k] , M[k] and T[k] are the ET, CMD and temperature for timestep k,
    c_1 and c_2 are parameters
    u[k] effective rainfall

    P[k] is the rainfall depth for timestep k and
    c_3 and c_4 are parameters
    For M[k] < 0, u[k] is
    increased by M k and M k set to zero
            M[k-1] = 100 # 100
            P[k] = [1,2,3,4]
            T[k] = [1,2,3,4]
            d = 200 # 200
            e = 0.166 # 0.166
            f = 1. # ? 0.5-1.3

    (1) of "Evaluation of streamflow predictions by the IHACRES rainfall-runoff model in two South African catchments"
    """
    if M_k_m_1 < d:
        M_k = M_k_m_1 * np.exp(-1.0 * P_k / d)
    elif M_k_m_1 < d + P_k:
        M_k = d * np.exp(-1.0 * (P_k - M_k_m_1 + d) / d)
    else:
        M_k = M_k_m_1 - P_k
    # End if

    f = f * d

    if T_k < 0:
        E_k = 0
    elif M_k < f:
        E_k = e * T_k
    else:
        E_k = e * T_k * np.exp(2.0 * (1.0 - M_k / f))
    # End if

    M_k = M_k + E_k
    U_k = P_k + M_k - M_k_m_1

    return U_k, M_k, E_k

# End cmd()


def routing(volume, storage_coef, node_inflow, local_inflow, irrig_ext, gamma=0.0):
    """Linear routing used to convert effective rainfall into streamflow for a given time step.

    :param volume: float, local water volume
    :param storage_coef: float, unknown parameter
    :param node_inflow: float, incoming streamflow
    :param local_inflow: float, unknown
    :param irrig_ext: float, volume of irrigation extraction in ML.
    :param gamma: float, unknown parameter which is always set to 0.0 in the original Fortran implementation.

    :returns: tuple[float], (local volume, and streamflow in unknown unit)
    """
    tmp_vol = volume + (node_inflow + local_inflow + gamma) - irrig_ext

    if tmp_vol > 0.0:
        volume = 1 / (1 + storage_coef) * tmp_vol
        # volume = (1 / (1 + storage_coef)) * tmp_vol
        outflow = storage_coef * volume
    else:
        volume = tmp_vol
        outflow = 0.0
    # End if

    return volume, outflow
# End routing()


def tau_v_to_ab(tau, v):

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


def run(x):  # dates, rain, mean_temp
    M_k = 100.0  # Initial Catchment moisture deficit (CMD)
    U = []  # effective rainfall
    M = []  # CMD for each timestep
    E = []  # evapotranspiration (ET) for each time step
    for k in range(len(dates)):
        U_k, M_k, E_k = cmd(
            M_k_m_1=M_k,  # 100
            P_k=rain[k],
            T_k=mean_temp[k],
            d=x[0],  # 200, flow threshold
            e=x[1],  # 0.166, calibrated constant?
            f=x[2],  # ? 0.5-1.3
        )
        U.append(U_k)
        M.append(M_k)
        E.append(E_k)
    # End for

    # flow (using IHACRES)
    a, b = tau_v_to_ab(tau=x[3:5], v=[x[5], x[6] * (1.0 - x[2])])

    Q = [1.0, 1.0]
    for k in range(2, len(dates)):
        Q_k = arma(a, b, Q, U, k)
        Q.append(Q_k)

    return U, M, E, Q
# End run()


def objective(x):
    U, M, E, Q = run(x)
    return np.sum((flows - Q)**2) / obj_denominator
# End objective()
