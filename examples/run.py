from pyIHACRES.run import run_ihacres

# Here we leverage the fact that everything in Python is an object
# including functions. We assign the model parameters to the function
# such that the function itself represents a node in a stream network.
# Of course, you could (should?) define a Class instead.
run_ihacres.d = 200.0  # flow threshold
run_ihacres.d2 = 2.0   # flow threshold, multiplier applied to `d`
run_ihacres.e = 0.1    # temperature to PET conversion factor
run_ihacres.f = 0.1    # plant stress threshold factor (applied to `d`). Determines effective rainfall.
run_ihacres.a = 54.35  # quickflow scaling factor
run_ihacres.b = 0.012  # slowflow scaling factor
run_ihacres.alpha = 0.727   #  effective rainfall scaling factor
run_ihacres.s = 2.5  # groundwater store factor
run_ihacres.catchment_area = 100.0  # area in km^2


# Initial conditions
cmd, quickflow, slowflow = 214.65, 0, 0

# Assume these are all 0.0
inflow = 0.0
gw_exchange = 0.0
extraction = 0.0

# Ideally, these would be read in via Numpy or Pandas
rainfall_ts = [70.0, 10.0, 0.0, 0.0, 200.0]
evaporation_ts = [2.0, 6.5, 7.0, 5.0, 1.0]

# Set up arrays to record state
outflow = [None] * len(rainfall_ts)
gw_state = [None] * (len(rainfall_ts)+1)
gw_state[0] = 0.0

# Run model (this would be in its own function)
for i in range(len(outflow)):
    progress = run_ihacres(cmd, rainfall_ts[i], evaporation_ts[i], inflow, quickflow, slowflow,
                           gw_state[i], gw_exchange, extraction)
    quickflow, slowflow, cmd = progress["state"]
    gw_state[i+1], outflow[i] = progress["flow"]

# Remove initial gw state to line up records
gw_state = gw_state[1:]

print(outflow)