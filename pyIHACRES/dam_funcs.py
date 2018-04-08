"""Functions to represent dam storage and discharge"""


def calc_dam_level(volume):
    """Calculate dam level from its volume.

    Adapted from Fortran code by Barry Croke.

    :math:`l = 156.8 + (0.9463 * v**0.2922)`

    where
    * :math:`l` is dam level
    * :math:`v` is dam volume
    * 156.8 is a base ground level in mAHD
    * 0.9463 is a coefficient representing the maximum capacity of the dam
    * 0.2922 - unsure what this is. Area coefficient?

    :param volume: float, volume of water in the dam

    :returns: float, dam level in mAHD
    """
    level = 156.8 + 0.9463 * volume**0.2922

    return level
# End calc_dam_level()


def calc_dam_area(volume):
    """Calculate dam area from its volume.

    Adapted from Fortran code by Barry Croke.

    :math:`A = 0.0021 * v**0.762`

    where
    * :math:`A` is dam area
    * :math:`v` is dam volume
    * 0.0021 was not specified
    * 0.762 was not specified

    :param volume: float, volume of water in the dam

    :returns: float, dam area in m^2 (?)
    """
    area = 0.0021 * volume**0.762

    return area
# End calc_dam_area()


def calc_dam_discharge(volume, max_storage):
    """Calculate dam discharge based on its volume and maximum storage.

    Adapted from Fortran code by Barry Croke.

    .. math::
        D = 0.0 if V <= S else 0.001492 * (V - S)**1.5280`

    where
    * :math:`D` is dam discharge
    * :math:`V` is dam volume
    * :math:`S` is dam storage
    * 0.001492 was not specified
    * 1.5280 was not specified

    :param volume: float, volume of water in the dam
    :param max_storage: float, maximum storage

    :returns: float, dam discharge (unknown unit)
    """
    discharge = 0.0
    if volume > max_storage:
        discharge = 0.001492 * (volume - max_storage)**1.5280
    # End if

    return max(0.0, discharge)
# End calc_dam_discharge()


def dam_volume_update(volume, node_inflow, gamma, rain, evap, infiltration, area, extractions, discharge):
    """Update dam volume for timestep

    gamma is gw exchange
    """
    vol = volume + (node_inflow + gamma) + (rain - evap - infiltration) * area - extractions - discharge

    return max(0.0, vol)
# End dam_volume_update()


def calc_dam_outflow(discharge, irrigation_extraction):
    """Calculate outflow from dam
    """
    return discharge + irrigation_extraction
# End calc_dam_outflow()


if __name__ == '__main__':

    volume = 304655.0
    irrig_ext = 3.0
    other_extraction = 2.0

    level = calc_dam_level(volume)
    area = calc_dam_area(volume)
    discharge = calc_dam_discharge(volume, 304650.0)
    volume = dam_volume_update(volume, 10, 0.3, 20.0, 3.0, 0.0, area, irrig_ext + other_extraction, discharge)

    print("""
    Dam Level: {}
    Dam Area: {}
    Dam Discharge: {}
    Dam Volume: {}""".format(level, area, discharge, volume))

    outflow = calc_dam_outflow(discharge, irrig_ext)

    print("Adj Dam Vol: {}".format(volume - outflow))

    print("Dam outflow: {}".format(outflow))

    # output(2,timestep,node)=outflow
    # output(1,timestep,node)=level
    # if (timestep.lt.timesteps) states(5,timestep+1,node)=volume
