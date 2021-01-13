# plot the Pc versus diff
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

fsize = 15

file_dir = os.getcwd() + '/'
master_dir = file_dir + '/../../Results' # gets master directory
img_dir = master_dir + '/'
Ethane_4Bar = pd.read_csv(img_dir + 'Ethane_298K_4Bar_CH3_1A/cm3overcm3/' + "RF_VarImp.csv")
# change "-inf" to -9
Ethane_4Bar['midpoints'][0] = -9
Ethane_4Bar['midpoints'][Ethane_4Bar['midpoints'] == 2] = 1
Ethane_4Bar = Ethane_4Bar.sort_values(by = ['midpoints'])

Ethane_20Bar = pd.read_csv(img_dir + 'Ethane_298K_20Bar_CH3_1A/cm3overcm3/' + "RF_VarImp.csv")
Ethane_20Bar['midpoints'][0] = -9
Ethane_20Bar['midpoints'][Ethane_20Bar['midpoints'] == 2] = 1
Ethane_20Bar = Ethane_20Bar.sort_values(by = ['midpoints'])

Ethane_40Bar = pd.read_csv(img_dir + 'Ethane_298K_40Bar_CH3_1A/cm3overcm3/' + "RF_VarImp.csv")
Ethane_40Bar['midpoints'][0] = -9
Ethane_40Bar['midpoints'][Ethane_40Bar['midpoints'] == 2] = 1
Ethane_40Bar = Ethane_40Bar.sort_values(by = ['midpoints'])

Propane_1Bar = pd.read_csv(img_dir + 'Propane_298K_1Bar_CH3_1A/cm3overcm3/' + "RF_VarImp.csv")
Propane_1Bar['midpoints'][0] = -9
Propane_1Bar['midpoints'][Propane_1Bar['midpoints'] == 2] = 1
Propane_1Bar = Propane_1Bar.sort_values(by = ['midpoints'])

Propane_5Bar = pd.read_csv(img_dir + 'Propane_298K_5Bar_CH3_1A/cm3overcm3/' + "RF_VarImp.csv")
Propane_5Bar['midpoints'][0] = -9
Propane_5Bar['midpoints'][Propane_5Bar['midpoints'] == 2] = 1
Propane_5Bar = Propane_5Bar.sort_values(by = ['midpoints'])

Propane_10Bar = pd.read_csv(img_dir + 'Propane_298K_10Bar_CH3_1A/cm3overcm3/' +"RF_VarImp.csv")
Propane_10Bar['midpoints'][0] = -9
Propane_10Bar['midpoints'][Propane_10Bar['midpoints'] == 2] = 1
Propane_10Bar = Propane_10Bar.sort_values(by = ['midpoints'])

linewidth = 5
highlight = 10
f, (ax,axes) = plt.subplots(1,2, sharey = True, figsize = (20, 6))
f.subplots_adjust(wspace = 0.05)

ax.plot(Ethane_4Bar['midpoints'], Ethane_4Bar['%IncMSE'], 
          label='Ethane at 4 Bar', color = 'blue', linewidth = linewidth)
ax.plot(Ethane_20Bar['midpoints'], Ethane_20Bar['%IncMSE'], 
          label='Ethane at 20 Bar', color = 'green', linewidth = linewidth)
ax.plot(Ethane_40Bar['midpoints'], Ethane_40Bar['%IncMSE'], 
          label='Ethane at 40 Bar', color = 'red', linewidth = linewidth)
# highlight circles, top 3
ax.scatter(Ethane_4Bar.sort_values(by = ['%IncMSE'], ascending=False)[0:highlight]['midpoints'], 
        Ethane_4Bar.sort_values(by = ['%IncMSE'], ascending=False)[0:highlight]['%IncMSE'], 
        edgecolors = 'blue', s = 100, facecolors='none', linewidth = linewidth/2)
ax.scatter(Ethane_20Bar.sort_values(by = ['%IncMSE'], ascending=False)[0:highlight]['midpoints'], 
        Ethane_20Bar.sort_values(by = ['%IncMSE'], ascending=False)[0:highlight]['%IncMSE'], 
        edgecolors = 'green', s = 100, facecolors='none', linewidth = linewidth/2)
ax.scatter(Ethane_40Bar.sort_values(by = ['%IncMSE'], ascending=False)[0:highlight]['midpoints'], 
        Ethane_40Bar.sort_values(by = ['%IncMSE'], ascending=False)[0:highlight]['%IncMSE'], 
        edgecolors = 'red', s = 100, facecolors='none', linewidth = linewidth/2)
ax.axvline(x = -8.5, linestyle = '--')
ax.axvline(x = 0.5, linestyle = '--')

x=list(range(5))
xticks=list(np.arange(-10, 2 + 2, 2))
yticks=list(np.arange(0, 35 + 5, 5))
xlabels=xticks.copy()
xlabels[-1] = r"$\mathbf{\infty}$"
xlabels[0] = r"-$\mathbf{\infty}$"
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xticklabels(xlabels, fontsize=fsize, fontweight = 'bold')
ax.set_yticklabels(yticks, fontsize=fsize, fontweight = 'bold')

#axes.ticklabel_format(axis = 'x', style = 'sci')
#axes.get_xticklabels()[-1].set_fontsize(fsize)
#axes.get_xticklabels()[0].set_fontsize(fsize)

ax.set_xlabel("Energies [kJ/mol]", fontsize = fsize, fontweight='bold')
ax.set_ylabel("Feature %IncMSE", fontsize = fsize, fontweight='bold')
ax.legend(loc = 'upper center', fontsize = fsize)
#plt.show()

#f.savefig(file_dir + 'Variable-Importance-Ethane.png', dpi=400, bbox_inches = 'tight')

axes.plot(Propane_1Bar['midpoints'], Propane_1Bar['%IncMSE'], 
          label='Propane at 1 Bar', color = 'blue', linewidth = linewidth)
axes.plot(Propane_5Bar['midpoints'], Propane_5Bar['%IncMSE'], 
          label='Propane at 5 Bar', color = 'green', linewidth = linewidth)
axes.plot(Propane_10Bar['midpoints'], Propane_10Bar['%IncMSE'], 
          label='Propane at 10 Bar', color = 'red', linewidth = linewidth)

# highlight circles
axes.scatter(Propane_1Bar.sort_values(by = ['%IncMSE'], ascending=False)[0:highlight]['midpoints'], 
        Propane_1Bar.sort_values(by = ['%IncMSE'], ascending=False)[0:highlight]['%IncMSE'], 
        edgecolors = 'blue', s = 100, facecolors='none', linewidth = linewidth/2)
axes.scatter(Propane_5Bar.sort_values(by = ['%IncMSE'], ascending=False)[0:highlight]['midpoints'], 
        Propane_5Bar.sort_values(by = ['%IncMSE'], ascending=False)[0:highlight]['%IncMSE'], 
        edgecolors = 'green', s = 100, facecolors='none', linewidth = linewidth/2)
axes.scatter(Propane_10Bar.sort_values(by = ['%IncMSE'], ascending=False)[0:highlight]['midpoints'], 
        Propane_10Bar.sort_values(by = ['%IncMSE'], ascending=False)[0:highlight]['%IncMSE'], 
        edgecolors = 'red', s = 100, facecolors='none', linewidth = linewidth/2)

axes.axvline(x = -8.5, linestyle = '--')
axes.axvline(x = 0.5, linestyle = '--')

x=list(range(5))
xticks=list(np.arange(-10, 2 + 2, 2))
xlabels=xticks.copy()
xlabels[-1] = r"$\mathbf{\infty}$"
xlabels[0] = r"-$\mathbf{\infty}$"
axes.set_xticks(xticks)
axes.set_xticklabels(xlabels, fontsize=fsize, fontweight = 'bold')

#axes.ticklabel_format(axis = 'x', style = 'sci')
#axes.get_xticklabels()[-1].set_fontsize(fsize)
#axes.get_xticklabels()[0].set_fontsize(fsize)

axes.set_xlabel("Energies [kJ/mol]", fontsize = fsize, fontweight='bold')
#plt.ylabel("Variable Importance", fontsize = fsize, fontweight='bold')
axes.legend(loc = 'upper center', fontsize = fsize)
#axes.show()

f.savefig(file_dir + 'Figure_8.png', dpi=400, bbox_inches = 'tight')
f.clf()