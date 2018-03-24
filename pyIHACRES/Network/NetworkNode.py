class NetworkNode(object):

    @property
    def storage(self):
        return self._storage[-1]
    # End storage()

    @storage.setter
    def storage(self, value):
        self._storage.append(value)
    # End storage()

    @property
    def effective_rainfall(self):
        return self._effective_rainfall[-1]
    # End effective_rainfall()

    @effective_rainfall.setter
    def effective_rainfall(self, value):
        self._effective_rainfall.append(value)
    # End storage()

    @property
    def et(self):
        return self._et[-1]
    # End et()

    @et.setter
    def et(self, value):
        self._et.append(value)
    # End et()

    @property
    def inflow(self):
        # TODO: Subtle bug exist here, wherein mismatch runs from different
        #       parent nodes can cause strange summed values
        return self._inflow[-1]
    # End inflow()

    @inflow.setter
    def inflow(self, ts_value):
        """Node inflow.

        :param ts_value: tuple, (int, float) indicating the time step and the inflow amount.
        """
        self.append_timestep(self._inflow, ts_value)
    # End inflow.setter()

    @property
    def outflow(self):
        return self._outflow[-1]
    # End outflow()

    @outflow.setter
    def outflow(self, ts_value):
        """Node outflow.

        :param ts_value: tuple, (int, float) indicating the time step and the outflow amount.
        """
        self.append_timestep(self._outflow, ts_value)
    # End outflow.setter()

    def get_outflow(self, timestep):
        """Get outflow for a given timestep

        :param timestep: int, 0-index value indicating timestep
        """
        return self._outflow[timestep]
    # End get_outflow()

    def append_timestep(self, attribute, ts_value):
        timestep, value = ts_value
        try:
            attribute[timestep] = value
        except IndexError:
            attribute.append(value)
            if len(attribute) != (timestep + 1):
                err_msg = "Timestep mismatch detected!\n"
                err_msg += "Attempted to assign value for time step '{}', but array length was '{}'"\
                           .format(timestep, len(attribute))
                raise IndexError(err_msg)
            # End if
        # End try
    # End append_timestep()

    def reset(self):
        # self.a = a
        # self.b = b
        # self.h0 = h0
        if hasattr(self, '_def_col'):
            self.d = self._def_col[0]
            self.e = self._def_col[1]
            self.f = self._def_col[2]

        self._outflow = []
        self._storage = [self._initial_storage]  # storage volume or Catchment Moisture Deficit
        self._effective_rainfall = []
        self._et = []  # evapotranspiration
        self._inflow = []  # inflow from each parent node
    # End reset()

# End NetworkNode()
