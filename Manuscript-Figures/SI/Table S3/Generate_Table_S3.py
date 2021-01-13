# plot the Pc versus diff
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

fsize = 15

file_dir = os.getcwd() + '/'
master_dir = file_dir + '/../../../Results' # gets master directory
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

def process_importances(data, condition):
  # first step, sort by importance, or %IncMSE
  data = data.sort_values(by = ['%IncMSE'], ascending = False)
  # second step, get the top 10, finally sort by midpoints, ascendingly
  Dataset = data.head(10).sort_values(by = ['midpoints'])
  # drop not needed columns
  Dataset = Dataset.drop(columns = ['midpoints', 'IncNodePurity'])
  # round %IncMSE
  Dataset['%IncMSE'] = Dataset.round({'%IncMSE': 2})
  # finally, rename the columns
  Dataset.columns = ['Histogram_Bin' if x=='varnames' else x for x in Dataset.columns]
  newnames = condition + '_' + Dataset.columns
  Dataset.columns = newnames
  Dataset = Dataset.reset_index(drop = True)
  return Dataset

sheets = process_importances(Ethane_4Bar, "Ethane_4Bar")
Temp = process_importances(Ethane_20Bar, "Ethane_20Bar")
# https://pandas.pydata.org/pandas-docs/stable/user_guide/merging.html
sheets = pd.concat([sheets, Temp.reindex(sheets.index)], axis=1)

Temp = process_importances(Ethane_40Bar, "Ethane_40Bar")
sheets = pd.concat([sheets, Temp.reindex(sheets.index)], axis=1)

Temp = process_importances(Propane_1Bar, "Propane_1Bar")
sheets = pd.concat([sheets, Temp.reindex(sheets.index)], axis=1)

Temp = process_importances(Propane_5Bar, "Propane_5Bar")
sheets = pd.concat([sheets, Temp.reindex(sheets.index)], axis=1)

Temp = process_importances(Propane_10Bar, "Propane_10Bar")
sheets = pd.concat([sheets, Temp.reindex(sheets.index)], axis=1)

# finally, save it
sheets.to_csv(file_dir + '/' + 'Feature-Importance.csv', index = False)
#sheets.to_excel(file_dir + '/' + 'Feature-Importance.xlsx', index = False)
#Temp['Propane_10Bar_varnames'] = Temp['Propane_10Bar_varnames'].astype("|S")