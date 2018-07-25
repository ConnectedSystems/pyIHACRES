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

    # print("a, b", a, b)

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


if __name__ == '__main__':

    import os
    import json
    import datetime

    site_id = "405229"
    with open(os.path.join("sw_data", site_id + '.json')) as f:
        site_data = json.loads(f.read())

    site_details = site_data["site_details"]
    dates = np.array([datetime.datetime.strptime(d, "%d/%m/%Y %H:%M:%S") for d in site_data["dates"]])
    flows = np.array(site_data["flows"])

    print(dates)
    print(flows)
    raise TypeError("debug exit")

    catchments_area = 1e6 * float(site_details['Drainage Area (Sq.km)'])
    rain = np.array(site_data["closest_rain"]) * catchments_area
    mean_temp = (np.array(site_data["closest_max_temp"]) + np.array(site_data["closest_min_temp"])) / 2.
    max_temp = np.array(site_data["closest_max_temp"])
    obj_denominator = np.sum((flows - np.mean(flows))**2)

    # 1000 m3 = 1 Ml

    from scipy import optimize
    import matplotlib.pyplot as plt
    # x_0 = [-1.7836, 0.7859, 20.4317, 19.9046, 200, 0.166, 0.9]
    # x_0 = [-1.7836, 0.7859, 20.4317, 19.9046]

    # `x` values appears to be:
    # x_0 = [d, e, f, tau_q, tau_s, v_q, v_s]
    # where
    # * tau_q < tau_s
    # * v_q + v_s = 1.0

    # tau_v_to_ab(tau=x[3:5], v=[x[5], x[6] * (1.0 - x[2])])

    x_0 = [200, 0.166, 0.9, 1.1, 30., 0.5, 0.5]
    # x_0 = [1.1, 30., 0.5, 0.5]
    # x_0 = [1.1, 0.5]
    # TODO ensure sum(vs) = 0
    res = optimize.minimize(fun=objective,
                            x0=x_0,
                            # method = 'CG',
                            method='L-BFGS-B',
                            # bounds = [(2., 3.), (0., 1.)]
                            bounds=[(150., 350.), (0.1, 0.2), (0.5, 1.3),
                                    (0.1, 2.), (10., 200.), (0., 1.), (0., 1.)]
                            )
    print res

    U, M, E, Q = run(res.x)
    plt.plot(Q, label='mod Q')
    # plt.plot(rain, label='rain')
    # plt.plot(max_temp, label='max_temp')
    # plt.plot(U, label='mod U')
    plt.plot(flows, '--', label='obs flow')
    plt.legend()
    plt.show()

    print "objective", objective(res.x)

    # import statsmodels.tsa.api as tsa
    # arma =tsa.ARMA(flows, order=(2,2), exog=U)
    # results= arma.fit()
    # print results.predict(30, 40, exog=U)
