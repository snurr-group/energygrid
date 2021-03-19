# plot the Pc versus diff
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from re import search

import PIL
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw

fsize = 8

file_dir = os.getcwd()
master_dir = file_dir + '/../../../All_data/' # gets master directory
f_dir = master_dir + '/'

rangemax = 20.0

list_files = ['XeKr_Mix_273K_1Bar', 'XeKr_Mix_273K_10Bar']


list_names = ["1 Bar", "10 Bar"]
list_labels = ["a", "b"]

list_img = []

F1bar = pd.read_excel(master_dir + 'SI.xlsx', 'XeKr_Mix_273K_1Bar')
F10bar = pd.read_excel(master_dir + 'SI.xlsx', 'XeKr_Mix_273K_10Bar')

F1bar = F1bar.rename(columns = {F1bar.columns[0] : "ID", F1bar.columns[1] : "Temp", 
                              F1bar.columns[2] : "Pres", F1bar.columns[3] : "Kr", 
                              F1bar.columns[4] : "Xe"})

F10bar = F10bar.rename(columns = {F10bar.columns[0] : "ID", F10bar.columns[1] : "Temp", 
                              F10bar.columns[2] : "Pres", F10bar.columns[3] : "Kr", 
                              F10bar.columns[4] : "Xe"})

F1bar['Select'] = (F1bar['Xe']/0.2)/(F1bar['Kr']/0.8)
F10bar['Select'] = (F10bar['Xe']/0.2)/(F10bar['Kr']/0.8)
# change the y-axis
F1bar.loc[F1bar['Select'] > rangemax, 'Select'] = rangemax
F10bar.loc[F10bar['Select'] > rangemax, 'Select'] = rangemax
fig, ax = plt.subplots()
#right = ax.spines['right']
#right.set_visible(False)
#up = ax.spines['top']
#up.set_visible(False)
#ax.set_xscale('log')
#ax.set_yscale('log')
line = plt.scatter(F1bar['Xe'], F1bar['Select'], color = 'orange', label = "1 Bar")
line2 = plt.scatter(F10bar['Xe'], F10bar['Select'], color = 'blue', label = "10 Bar")
line.set_clip_on(False)
line2.set_clip_on(False)
#plt.xticks(np.arange(0, rangemax + int(rangemax/4), int(rangemax/4)), fontsize = fsize)
plt.yticks(np.arange(0, rangemax + int(rangemax/4), int(rangemax/4)), fontsize = fsize)
plt.xlabel(r"GCMC Mixture Xe Loading [cm$^{\rm 3}$/cm$^{\rm 3}$]", fontsize = fsize, fontweight='bold')
plt.ylabel("GCMC Selectivity", fontsize = fsize, fontweight='bold')


#plt.axis('square')
plt.ylim([0, rangemax])
#plt.xlim([0, np.max(GCMC['Xe'])])

plt.legend()
plt.savefig(file_dir + '/'  + 'pareto.png', dpi=900, bbox_inches = 'tight')

#plt.clf()
