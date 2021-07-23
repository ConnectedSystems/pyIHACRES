## Functions related to estimating ET and U (evapotranspiration and effective rainfall)

from math import exp
import numpy as np


def calc_effective_rainfall(rainfall: float, cmd: float, d: float, d2: float, n: float=0.1) -> float:
    """
    Estimate effective rainfall.
    
    :References:
        Croke, B.F.W., Jakeman, A.J. 2004
            A catchment moisture deficit module for the IHACRES rainfall-runoff model,
            Environmental Modelling & Software, 19(1), pp. 1–5.
            doi: 10.1016/j.envsoft.2003.09.001
    
        Croke, B.F.W., Jakeman, A.J. 2005
            Corrigendum to "A Catchment Moisture Deficit module for the IHACRES
            rainfall-runoff model [Environ. Model. Softw. 19 (1) (2004) 1–5]"
            Environmental Modelling & Software, 20(7), p. 997.
            doi: https://doi.org/10.1016/j.envsoft.2004.11.004
    
    :Parameters:
        - rainfall : rainfall for time step
        - cmd      : previous CMD value
        - Mf       : interim CMD value
        - d        : threshold value
        - d2       : scaling factor applied to `d`
        - n        : scaling factor (default = 0.1)
                     Default value is suitable for most cases (Croke & Jakeman, 2004)
    
    :Returns:
        effective rainfall
    """
    d2: float = d * d2
    e_rainfall: float

    if (cmd > d2):
        e_rainfall = rainfall
    else:
        f1: float = np.minimum(1.0, cmd / d)
        f2: float = np.minimum(1.0, cmd / d2)
        e_rainfall = rainfall * ((1.0 - n) * (1.0 - f1) + (n * (1.0 - f2)))

    return np.max(0.0, e_rainfall)


def calc_ET_from_E(e: float, evap: float, Mf: float, f: float, d: float) -> float:
    """
    Calculate evapotranspiration from evaporation.
    
    :Parameters:
        - e    : temperature to PET conversion factor (a stress threshold)
        - evap : evaporation for given time step.
        - Mf   : Catchment Moisture Deficit prior to accounting for ET losses (`M_{f}`)
        - f    : calibrated parameter that acts as a multiplication factor on `d`
        - d    : flow threshold factor
    
    :Results:
        estimate of ET
    """
    param_g: float = f * d
    et: float = e * evap

    if Mf > param_g:
        et = et * np.minimum(1.0, exp((1.0 - Mf/param_g)*2.0))

    return np.maximum(0.0, et)


def calc_ET(e: float, evap: float, Mf: float, f: float, d: float) -> float:
    """Deprecated function - call calc_ET_from_E instead."""
    return calc_ET_from_E(e, evap, Mf, f, d)


def calc_ET_from_T(e: float, T: float, Mf: float, f: float, d: float) -> float:
    """
    Calculate evapotranspiration based on temperature data.
    
    Parameters `f` and `d` are used to calculate `g`, the value of the CMD
    which the ET rate will begin to decline due to insufficient
    water availability for plant transpiration.
    
    :Parameters:
        - e  : temperature to PET conversion factor (a stress threshold)
        - T  : temperature in degrees C
        - Mf : Catchment Moisture Deficit prior to accounting for ET losses (`M_{f}`)
        - f  : multiplication factor on `d`
        - d  : flow threshold factor
    
    :Returns:
        estimate of ET from temperature (for catchment area)
    """
    # temperature can be negative, so we have a min cap of 0.0
    if T <= 0.0:
        return 0.0

    param_g: float = f * d
    et: float = e * T * np.minimum(1.0, exp(2.0 * (1.0 - (Mf / param_g))) )

    return np.maximum(0.0, et)