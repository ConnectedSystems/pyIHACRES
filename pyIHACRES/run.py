"""Generic run methods."""

import numpy as np

from . import climate
from . import flow
from . import cmd as ihacres_cmd


def run_ihacres(cmd, rainfall, evaporation, inflow, quickflow, slowflow,
                gw_state, gw_exchange, extraction):
    '''Example run function.

    Parameters
    ----------
    cmd : catchment moisture deficit, $M_{k}$ in the literature.
    rainfall : precipitation ($P$)
    evaporation : ($E$)
    inflow : flow from previous node
    quickflow : quickflow at previous time step, $(V_{q-1})$
    slowflow : slowflow at previous time step, $(V_{s-1})$
    gw_state : groundwater storage index at $t-1$
    gw_exchange : volume flux, interaction with groundwater
    extraction : volume extracted from stream
    '''
    d, d2, e, f = run_ihacres.d, run_ihacres.d2, run_ihacres.e, run_ihacres.f
    a, b, alpha, s = run_ihacres.a, run_ihacres.b, run_ihacres.alpha, run_ihacres.s
    catchment_area = run_ihacres.catchment_area

    mf, U, r = ihacres_cmd.calc_ft_interim_cmd(cmd, rainfall, d, d2, alpha)
    ET = climate.calc_ET(e, evaporation, mf, f, d)
    cmd = ihacres_cmd.calc_cmd(cmd, rainfall, ET, U, r)

    Vq, Vs, outflow = flow.calc_ft_flows(quickflow, slowflow, U, r,
                                         catchment_area, a, b)

    # if node routes to another node
    gw_state, outflow = flow.routing(gw_state, s, inflow, outflow, extraction, gw_exchange)

    return {
        "flow": (gw_state, outflow),
        "state": (Vq, Vs, cmd)
    }