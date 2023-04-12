import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

# Define the one-compartment pharmacokinetic model
def one_compartment_model(param, t, dose, weight, stop_time):
    # Extract the parameters
    Vd, Kel = param
    
    # Calculate the clearance
    Cl = Kel * weight
    
    # Initialize an array to store the ketamine concentration at each time step
    Cp = np.zeros(len(t))
    
    # Loop through each time step
    for i, t_i in enumerate(t):
        # If it's the first time step
        if i == 0:
            # Calculate the initial ketamine concentration using the initial dose
            Cp[i] = dose[i] / Vd
        else:
            # If the infusion is still ongoing
            if t_i < stop_time:
                # Calculate the current ketamine concentration using the previous concentration and the new dose
                Cp[i] = (dose[i] / Vd) + (Cp[i-1] * np.exp(-Kel * (t_i - t[i-1])))
            # If the infusion has stopped
            else:
                # Calculate the current ketamine concentration using only the previous concentration
                Cp[i] = Cp[i-1] * np.exp(-Kel * (t_i - t[i-1]))
    # Return the array of ketamine concentrations
    return Cp

# Define the time period for simulation
t = np.arange(0, 120, 10)

# Define the time the infusion stops
stop_time = 40 # minutes

# Define the weight of the patient
weight = 70 # kg

# Define the infusion amounts at each time step
dose = np.zeros(len(t))
for i, t_i in enumerate(t):
    if t_i < stop_time:
        dose[i] = 2 * weight # mg/kg/hr * weight (kg)
    else:
        dose[i] = 0

# Define the initial guesses for the parameters
Vd0 = weight
Kel0 = 1 / weight
param0 = [Vd0, Kel0]

# Compute the estimated ketamine concentration over time
Cp = one_compartment_model(param0, t, dose, weight, stop_time)

# Apply a moving average filter to the data
Cp_smooth = savgol_filter(Cp, 5, 4)

# Plot the estimated ketamine concentration over time
plt.plot(t, Cp_smooth)
plt.xlabel('Time (minutes)')
plt.ylabel('Ketamine Concentration (mg/L)')
plt.show()

