from __future__ import division

from .NetworkNode import NetworkNode

from .. import ihacres_funcs


class StreamNode(NetworkNode):

    def __init__(self, node_id, prev_node, next_node, formula_type, def_col, initial_storage, storage_coef, alpha, a, b, area):
        """Stream node"""
        self.node_id = node_id
        self.prev_node = {k: None for k in prev_node} if prev_node else {}
        self.next_node = next_node
        self.formula_type = formula_type
        self.area = area

        self.set_calib_params(def_col, storage_coef, alpha, a, b, flow_mod=1.0)

        self._quickflow = [0.0]
        self._slowflow = [0.0]
        self._outflow = []

        self._initial_storage = initial_storage
        self._storage = [initial_storage]  # storage volume or Catchment Moisture Deficit
        self._effective_rainfall = []
        self._et = []  # evapotranspiration
        self._inflow = []  # inflow from each parent node
    # End init()

    @property
    def quickflow(self):
        return self._quickflow[-1]
    # End quickflow()

    @property
    def slowflow(self):
        return self._slowflow[-1]
    # End slowflow()

    def set_calib_params(self, def_col, s_coef, alpha, a, b, flow_mod=None):
        self.d = def_col[0]
        self.d2 = def_col[1]
        self.e = def_col[2]
        self.f = def_col[3]
        self.storage_coef = s_coef
        self.alpha = alpha
        self.a = a
        self.b = b

        self.flow_mod = flow_mod

    # End set_calib_params()

    def get_flows(self):
        return self._quickflow[-1], self._slowflow[-1], self._outflow[-1]

    def update_state(self, timestep, storage, effective_rainfall, et, qflow_store, sflow_store, outflow):
        self.storage = (timestep, storage)
        self.effective_rainfall = effective_rainfall
        self.et = et
        self.outflow = (timestep, outflow)

        self.append_timestep(self._quickflow, (timestep, qflow_store))
        self.append_timestep(self._slowflow, (timestep, sflow_store))
    # End update_state()

    def run(self, timestep, rain_evap, extractions):
        """Run node to calculate outflow and update state.

        :param timestep: int, time step
        :param rain_evap: np.ndarray, rainfall and evapotranspiration data for all nodes
        :param extractions: np.ndarray, irrigation and other water extractions for all nodes

        :returns: float, outflow from node
        """
        # try:
        #     outflow = self.get_outflow(timestep)
        # except IndexError:
        #     pass
        # # End try
        rainfall = rain_evap["{}_rain".format(self.node_id)]
        evap = rain_evap["{}_evap".format(self.node_id)][timestep]

        # other extractions are ignored for stream nodes, so only extract irrigation ext.
        irrig_ext = extractions["{}_irrig".format(self.node_id)]
        ext = irrig_ext[timestep]
        ts_rainfall = rainfall[timestep]

        Mf, e_rainfall, recharge = ihacres_funcs.calc_ft_interim(self.storage, ts_rainfall, self.d,
                                                                 self.d2, self.alpha)

        et = ihacres_funcs.calc_ET(self.e, evap, Mf, self.f, self.d)
        cmd = ihacres_funcs.calc_cmd(self.storage, ts_rainfall, et, e_rainfall, recharge)

        inflow = 0.0
        for nid in self.prev_node:
            inflow += self.prev_node[nid].run(timestep, rain_evap, extractions)
        # End for
        self.inflow = (timestep, inflow)

        quick_store, slow_store, outflow = ihacres_funcs.calc_ft_flows(self.quickflow, self.slowflow,
                                                                       e_rainfall, recharge, self.area,
                                                                       self.a, self.b, loss=0.0)

        # DEBUG: modifier
        if self.flow_mod:
            quick_store = quick_store / self.flow_mod
            slow_store = slow_store / self.flow_mod
            outflow = outflow / self.flow_mod
            # print("Modded Quick, slow, outflow")
            # print(quick_store, slow_store, outflow)

        if self.next_node and ('dam' not in type(self.next_node).__name__.lower()):
            cmd, outflow = ihacres_funcs.routing(cmd, self.storage_coef, inflow, outflow, ext, gamma=0.0)
        else:
            outflow = ihacres_funcs.calc_outflow(outflow, ext)
        # End if

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

        self.update_state(timestep, cmd, e_rainfall, et, quick_store, slow_store, outflow)

        return outflow
    # End run()

# End StreamNode()
