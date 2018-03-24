import numpy as np

from NetworkNode import NetworkNode

from .. import ihacres_funcs


class StreamNode(NetworkNode):

    def __init__(self, node_id, prev_node, next_node, formula_type, def_col, initial_storage, storage_coef):
        """next_node is intended to be a reference to the next StreamNode object"""
        self.node_id = node_id
        self.prev_node = {k: None for k in prev_node} if prev_node else {}
        self.next_node = next_node
        self.formula_type = formula_type
        # self.a = a
        # self.b = b
        # self.h0 = h0
        self._def_col = def_col  # stored to enable resetting
        self.d = def_col[0]
        self.e = def_col[1]
        self.f = def_col[2]
        self.storage_coef = storage_coef
        self._outflow = []

        self._initial_storage = initial_storage
        self._storage = [initial_storage]  # storage volume or Catchment Moisture Deficit
        self._effective_rainfall = []
        self._et = []  # evapotranspiration
        self._inflow = []  # inflow from each parent node
    # End init()

    def update_state(self, timestep, storage, effective_rainfall, et, outflow):
        self.storage = storage
        self.effective_rainfall = effective_rainfall
        self.et = et
        self.outflow = (timestep, outflow)
    # End update_state()

    # def run(self, timestep, rain, et, irrig_ext, extractions):
    #     # extractions are ignored for stream nodes
    #     e_rain, storage = ihacres_funcs.ihacres_cmd(M_k=self.storage, P_k=rain, E_k=et, d=self.d, e=self.e, f=self.f)
    #     storage, outflow = ihacres_funcs.routing(storage, self.storage_coef, self.inflow, e_rain, irrig_ext, gamma=0.0)
    #
    #     # storage, effective_rainfall, outflow = self.calc_outflow(rain, et, irrig_ext + extractions)
    #     self.update_state(timestep, storage, e_rain, et, outflow)
    # # End run()

    def run(self, timestep, rain_et, extractions):
        """Run node to calculate outflow.

        :param timestep: int, time step
        :param rain_et: np.ndarray, rainfall and evapotranspiration data for all nodes
        :param extractions: np.ndarray, irrigation and other water extractions for all nodes

        :returns: float, outflow from node
        """
        try:
            outflow = self.get_outflow(timestep)
        except IndexError:
            rain = rain_et["{}_rain".format(self.node_id)][timestep]
            et = rain_et["{}_evap".format(self.node_id)][timestep]

            # other extractions are ignored for stream nodes, so only extract irrigation ext.
            irrig_ext = extractions["{}_irrig".format(self.node_id)]
            ext = irrig_ext[timestep]

            e_rain, storage, et = ihacres_funcs.ihacres_cmd(
                M_k_m_1=self.storage, P_k=rain, E_k=et, d=self.d, e=self.e, f=self.f)

            inflow = 0.0
            for nid in self.prev_node:
                inflow += self.prev_node[nid].run(timestep, rain_et, extractions)
            # End for
            self.inflow = (timestep, inflow)

            # TODO: Calc stream level
            # if self.formula_type == 1:
            #     waterlevel = 1.0 * np.exp

            #       if (formula.eq.1) then
            # c      write(*,*) 'i'
            #        waterlevel=1.0d0
            #      :  *exp(par(1))*(tmp_flow**par(2))
            #      :  *1.0d0/((1.0d0+(tmp_flow/par(3))**par(4))**(par(5)/par(4)))
            #      :  *exp(par(6)/(1+exp(-par(7)*par(8))*tmp_flow**par(7)))
            #      :  +CTF

            # routing(volume, storage_coef, node_inflow, local_inflow, irrig_ext, gamma=0.0)
            # not sure what local_inflow is...
            storage, outflow = ihacres_funcs.routing(storage, self.storage_coef, inflow, inflow, ext, gamma=0.0)

            self.update_state(timestep, storage, e_rain, et, outflow)
        # End try

        return outflow
    # End run()

# End StreamNode()
