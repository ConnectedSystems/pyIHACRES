from StreamNode import StreamNode

from .. import dam_funcs


class DamNode(StreamNode):

    """Dam node representation"""

    # node_id, next_node, formula_type, def_col, initial_storage, storage_coef
    def __init__(self, node_id, prev_node, next_node, initial_storage, max_storage, formula_type, storage_coef):
        self.node_id = node_id
        self.prev_node = {k: None for k in prev_node} if prev_node else {}
        self.next_node = next_node
        self.formula_type = formula_type
        self.max_storage = max_storage

        self.storage_coef = storage_coef

        self._initial_storage = initial_storage
        self._storage = [initial_storage]  # storage volume or Catchment Moisture Deficit
        self._effective_rainfall = []
        self._et = []  # evapotranspiration
        self._inflow = []  # inflow from each parent node

        self._level = []
        self._area = []
        self._discharge = []
        self._outflow = []
    # End init()

    @property
    def level(self):
        return self._level[-1]

    @level.setter
    def level(self, value):
        self._level.append(value)

    @property
    def area(self):
        return self._area[-1]

    @area.setter
    def area(self, value):
        self._area.append(value)

    @property
    def discharge(self):
        return self._discharge[-1]

    @discharge.setter
    def discharge(self, value):
        self._discharge.append(value)

    def update_state(self, timestep, storage, rainfall, et, area, discharge, outflow):
        self.storage = storage
        self.effective_rainfall = rainfall
        self.et = et
        self.level = dam_funcs.calc_dam_level(storage)
        self.outflow = (timestep, outflow)
    # End update_state()

    def run(self, timestep, rain_et, extractions):

        try:
            outflow = self.get_outflow(timestep)
        except IndexError:
            rain = rain_et["{}_rain".format(self.node_id)][timestep]
            et = rain_et["{}_evap".format(self.node_id)][timestep]
            irrig_ext = extractions["{}_irrig".format(self.node_id)][timestep]
            other = extractions["{}_other".format(self.node_id)][timestep]
            ext = irrig_ext + other

            area = dam_funcs.calc_dam_area(self.storage)
            discharge = dam_funcs.calc_dam_discharge(self.storage, self.max_storage)

            inflow = 0.0
            for nid in self.prev_node:
                inflow += self.prev_node[nid].run(timestep, rain_et, extractions)
            # End for
            self.inflow = (timestep, inflow)

            # volume, node_inflow, gamma, rain, evap, infiltration, area, extractions, discharge
            storage = dam_funcs.dam_volume_update(self.storage, inflow, 0.0, rain, et, 0.0,
                                                  area, ext, discharge)

            outflow = dam_funcs.calc_dam_outflow(discharge, irrig_ext)

            # area, discharge, storage, outflow = self.calc_outflow(rain, et, irrig_ext, extractions)
            self.update_state(timestep, storage, rain, et, area, discharge, outflow)
        # End try

        return outflow
    # End run()

# End DamNode()
