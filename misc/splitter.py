# here's an example of how one could split up and read in 
# the csv files that h52csv makes

import pandas as pd
import io
import matplotlib.pyplot as plt
import numpy as np


file_name = 'test_Run_0_1555509928.h5_A_4.csv'

f = open(file_name)
data = f.read()
f.close()

chunks = data.split('\n\n')

header = chunks[0]

del(chunks[0])

for chunk in chunks:
    this_buf = io.StringIO(chunk)
    what_this_is = this_buf.readline().strip()
    dframe = pd.read_csv(this_buf)
    i = np.array(dframe.current)
    v = np.array(dframe.voltage)
    t = np.array(dframe.time)
    s = np.array(dframe.status)

    if what_this_is == 'I_sc dwell':
        plt.figure()
        plt.title(what_this_is)
        plt.plot(t, i, '.')
        plt.ylabel('Current [A]')
        plt.xlabel('Time [s]')
        plt.grid()
    elif what_this_is in ['Snaith', 'Sweep']:
        plt.figure()
        plt.title(what_this_is)
        plt.plot(v, i, '.')
        plt.xlabel('Voltage [V]')
        plt.ylabel('Current [I]')
        plt.grid()
    elif what_this_is == 'MPPT':
        plt.figure()
        plt.title(what_this_is)
        plt.plot(t, i*v*1000, '.')
        plt.xlabel('Time [s]')
        plt.ylabel('Power [mW]')
        plt.grid()
    elif what_this_is == 'V_oc dwell':
        plt.figure()
        plt.title(what_this_is)
        plt.plot(t, v, '.')
        plt.ylabel('Voltage [V]')
        plt.xlabel('Time [s]')
        plt.grid()
    pass
plt.show(True)
