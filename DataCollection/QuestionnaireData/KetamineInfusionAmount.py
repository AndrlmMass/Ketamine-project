#Import relevant libraries
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
import numpy as np

#Define weight and doses
time_stamps = [0,5,10,15,20,25,30,35,40,45,50,55,60]

def ketamine_amount(doses, time_stamps):
    halflife = 150 #minutes
    ketamine_amount = []
    for i in range(len(time_stamps)):
        if i == 0:
            ketamine_amount.append(doses[i])
        else:
            ketamine_amount.append(ketamine_amount[i-1]*np.exp(-np.log(2)*time_stamps[i]/halflife) + doses[i])
    time_stamps = np.array(time_stamps)
    ketamine_amount = np.array(ketamine_amount)
    X_Y_Spline = make_interp_spline(time_stamps, ketamine_amount)
    # Returns evenly spaced numbers
    # over a specified interval.
    X_ = np.linspace(time_stamps.min(), time_stamps.max(), 500)
    Y_ = X_Y_Spline(X_)
    return X_, Y_

#SD5001
doses = [0,6.0,2.0,1.5,1.5,1.0,0,0,0,0,0,0,0]
SD01_x,SD01_y = ketamine_amount(doses, time_stamps)

#SD5002
doses = [0,21.0,1.5,1.0,1.0,0,0,0,0,0,0,0,0]
SD02_x,SD02_y = ketamine_amount(doses, time_stamps)

#SD5007
doses = [0,4.0,2.0,1.5,1.0,1.0,0,0,0,0,0,0,0]
SD07_x,SD07_y = ketamine_amount(doses, time_stamps)

#SD5008
doses = [0,5.0,2.0,1.0,1.0,1.0,0,0,0,0,0,0,0]
SD08_x,SD08_y = ketamine_amount(doses, time_stamps)

#SD5009
doses = [0,10.0,0,0,0,0,0,0,0,0,0,0,0]
SD09_x,SD09_y = ketamine_amount(doses, time_stamps)

#SD5010
doses = [0,7.0,1.0,1.5,1.0,0.5,0,0,0,0,0,0,0]
SD10_x,SD10_y = ketamine_amount(doses, time_stamps)

# Plotting the Graph
plt.plot(SD01_x, SD01_y, label = 'SD5001')
plt.plot(SD02_x, SD02_y, label = 'SD5002')
plt.plot(SD07_x, SD07_y, label = 'SD5007')
plt.plot(SD08_x, SD08_y, label = 'SD5008') 
plt.plot(SD09_x, SD09_y, label = 'SD5009') 
plt.plot(SD10_x, SD10_y, label = 'SD5010') 
plt.plot(SD01_x, (SD01_y+SD02_y+SD07_y+SD08_y+SD09_y+SD10_y)/6, label = 'Average Response', color = 'black', linestyle = '--', linewidth = 2)
plt.title('Ketamine Infusion Amount')
plt.xlabel('Time (minutes)')
plt.ylabel('Ketamine (mg/kg)')
plt.legend()
plt.show()




