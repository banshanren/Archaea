import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import rcParams

rcParams['savefig.dpi'] = 400
rcParams['font.family'] = 'Times New Roman'
rcParams['mathtext.default'] = 'regular'
# rcParams['text.usetex'] = True
font = {'family': 'Cambria',
        'weight': 'normal',
        'size': 12,
        }
result = pd.read_csv('./split_result.txt')

# Size = result['Size'].values
Ratio = result['Ratio'].values
ACC = result['ACC'].values
MCC = result['MCC'].values
AUC = result['AUC'].values
Sn = result['Sn'].values
Sp = result['Sp'].values
MAC = (Sn+Sp)/2
Sprod = Sn*Sp


fig, ax = plt.subplots() 
ax.plot(Ratio, ACC, color='orange', label='ACC')
ax.plot(Ratio, MCC, color='blue', label='MCC')
# ax.plot(Ratio, AUC, label='AUC')
ax.plot(Ratio, Sn, color='red', label='Sn')
ax.plot(Ratio, Sp, color='magenta', label='Sp')
ax.plot(Ratio, MAC, color='black', label='BACC')
# ax.plot(Ratio, Sprod, color='coral', label='Sprod')
ax.legend(loc=0,ncol=2)

# plt.ylim([0.3, 1.0])
plt.xlim([1, 2.1])
plt.tick_params(direction='in', top=True, right=True)
plt.xlabel('Ratio of negtive samples and positive samples',font)
plt.ylabel('Performance', font)

plt.show()

