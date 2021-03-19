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

fsize = 10

file_dir = os.getcwd()
master_dir = file_dir + '/../../../../Results' # gets master directory
f_dir = master_dir + '/'

rangemax = 20.0

#list_files = [f_dir + 'XeKr_Mix_273K_1Bar_Kr_1A_old/FitLoading/' + 'XeKr_Mix_273K_1Bar_Kr_1A_LASSO_train.csv']

list_files = [f_dir + 'XeKr_Mix_273K_1Bar_Kr_1A_old/' + 'XeKr_Mix_273K_1Bar_Kr_1A_LASSO_train.csv', 
              f_dir + 'XeKr_Mix_273K_1Bar_Kr_1A_old/' + 'XeKr_Mix_273K_1Bar_Kr_1A_RF_train_Select.csv', 
              f_dir + 'XeKr_Mix_273K_10Bar_Kr_1A_old/' + 'XeKr_Mix_273K_10Bar_Kr_1A_LASSO_train.csv', 
              f_dir + 'XeKr_Mix_273K_10Bar_Kr_1A_old/' + 'XeKr_Mix_273K_10Bar_Kr_1A_RF_train_Select.csv', 
              f_dir + 'XeKr_Mix_273K_1Bar_Kr_1A_old/' + 'XeKr_Mix_273K_1Bar_Kr_1A_LASSO_test.csv', 
              f_dir + 'XeKr_Mix_273K_1Bar_Kr_1A_old/' + 'XeKr_Mix_273K_1Bar_Kr_1A_RF_test_Select.csv', 
              f_dir + 'XeKr_Mix_273K_10Bar_Kr_1A_old/' + 'XeKr_Mix_273K_10Bar_Kr_1A_LASSO_test.csv', 
              f_dir + 'XeKr_Mix_273K_10Bar_Kr_1A_old/' + 'XeKr_Mix_273K_10Bar_Kr_1A_RF_test_Select.csv']

list_stats = [[0.42, 0.91, 0.71, 2.9], 
              [0.31, 0.94, 0.56, 3.2], 
              [0.1, 0.9, 1.1, 10.8], 
              [0.018, 0.95, 1.2, 12.8],
              [0.19, 0.92, 0.95, 5.9],
              [0.22, 0.96, 0.75, 5.8],
              [0.3, 0.87, 0.95, 2.5], 
              [0.0087, 0.93, 1.2, 7.9]]
list_names = ["1 Bar, LASSO", 
              "1 Bar, RF",
              "10 Bar, LASSO", 
              "10 Bar, RF",
              "1 Bar, LASSO", 
              "1 Bar, RF",
              "10 Bar, LASSO", 
              "10 Bar, RF"]

list_labels = ["a", "b", "c", "d", "e", "f", "g", "h"]

list_img = []
a = 0
for file in list_files:
    if(a > 0):
        break
    GCMC = pd.read_csv(file)
    if(search("LASSO", file)):
        columnname = 'lassopred'
        filename1 = "LASSO"
    else:
        columnname = 'Selectivity'
        filename1 = "RF"
    if(search("train", file)):
        color = 'orange'
        text1 = "Training data\n 1000 points"
        filename2 = "Train"
    else:
        color = 'blue'
        text1 = "Testing data\n 1000 points"
        filename2 = "Test"
    if(search("10Bar", file)):
        filename0 = "10Bar_"
    else:
        filename0 = "1Bar_"
    
    xlims = [0, rangemax]
    ylims = [0, rangemax]
    diag = np.linspace(0, xlims[1], int(rangemax/4))
    # change the x-axis
    GCMC.loc[GCMC['y_act'] > rangemax, 'y_act'] = rangemax
    # also change the y-axis
    GCMC.loc[GCMC[columnname] > rangemax, columnname] = rangemax
    fig, ax = plt.subplots()
    right = ax.spines['right']
    #right.set_visible(False)
    up = ax.spines['top']
    #up.set_visible(False)
    
    plt.plot(diag, diag, '--', color = 'black')
    line = plt.scatter(GCMC['y_act'], GCMC[columnname], color = color, linewidth = 1)
    line.set_clip_on(False)
    plt.xticks(np.arange(0, rangemax + int(rangemax/4), int(rangemax/4)), fontsize = fsize)
    plt.yticks(np.arange(0, rangemax + int(rangemax/4), int(rangemax/4)), fontsize = fsize)
    plt.xlabel(r"GCMC Selectivity", fontsize = fsize, fontweight='bold')
    plt.ylabel(r"Predicted Selectivity", fontsize = fsize, fontweight='bold')
    plt.text(0.5, 0.95, text1, ha = 'center', va = 'center', transform=ax.transAxes, fontsize = fsize, color = color)
    plt.text(0.75, 0.20, r"R$^{\rm 2}$ = " + str(list_stats[a][0]) + '\n' + 
             r"$^{\rho}$ = " + str(list_stats[a][1]) + '\n' + 
             'MAE = ' + str(list_stats[a][2]) + '\n' + 
             'RMSE = ' + str(list_stats[a][3]) + '\n', 
             ha = 'center', va = 'center', transform=ax.transAxes, fontsize = fsize, color = 'black')
    plt.axis('square')
    plt.xlim(xlims)
    plt.ylim(ylims)
    
    #plt.savefig(file_dir + '/' + filename0 + filename1 + filename2 + '_new.png', dpi=900, bbox_inches = 'tight')
    a+=1
    plt.clf()
    
    list_img.append(file_dir + '/' + filename0 + filename1 + filename2 + '_new.png')
imgs    = [ PIL.Image.open(i).convert("RGBA") for i in list_img ]
for a in range(0, len(imgs)):
  txt = Image.new("RGBA", imgs[a].size, (255,255,255,0))
  d = ImageDraw.Draw(txt)
  namefont = ImageFont.truetype("ariblk.ttf", 150)
  labelfont = ImageFont.truetype("ariblk.ttf", 200) # name for arial bold
  #https://stackoverflow.com/questions/1815165/draw-bold-italic-text-with-pil
  d.multiline_text((1200, 350),list_names[a],fill = (0,0,0,255),font=namefont)
  d.text((600, 300),list_labels[a],fill = (0,0,0,255),font=labelfont)
  
  imgs[a] = Image.alpha_composite(imgs[a], txt)
def get_concat_h(im1, im2):
    dst = Image.new('RGB', (im1.width + im2.width, im1.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    return dst

def get_concat_v(im1, im2):
    dst = Image.new('RGB', (im1.width, im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst

def get_concat_h_cut(im1, im2):
    dst = Image.new('RGB', (im1.width + im2.width, min(im1.height, im2.height)))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    return dst

def get_concat_v_cut(im1, im2):
    dst = Image.new('RGB', (min(im1.width, im2.width), im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst

left1 = get_concat_v_cut(imgs[0], imgs[1])
left2 = get_concat_v_cut(imgs[4], imgs[5])
left = get_concat_v_cut(left1, left2)
right1 = get_concat_v_cut(imgs[2], imgs[3])
right2 = get_concat_v_cut(imgs[6], imgs[7])
right = get_concat_v_cut(right1, right2)

get_concat_h_cut(left, right).save('Figure_S8.png')
