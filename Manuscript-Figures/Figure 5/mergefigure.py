import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

import PIL
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw

file_dir = os.getcwd()
master_dir = file_dir + '/../../Results' # gets master directory
img_dir = master_dir + '/'
list_im = [img_dir + 'Ethane_298K_4Bar_CH3_1A/cm3overcm3/' + 'Ethane_298K_4Bar_CH3_1A_LASSO_train.png', 
           img_dir + 'Ethane_298K_20Bar_CH3_1A/cm3overcm3/' + 'Ethane_298K_20Bar_CH3_1A_LASSO_train.png',
           img_dir + 'Ethane_298K_40Bar_CH3_1A/cm3overcm3/' + 'Ethane_298K_40Bar_CH3_1A_LASSO_train.png', 
           img_dir + 'Ethane_298K_4Bar_CH3_1A/cm3overcm3/' + 'Ethane_298K_4Bar_CH3_1A_LASSO_test.png',
           img_dir + 'Ethane_298K_20Bar_CH3_1A/cm3overcm3/' + 'Ethane_298K_20Bar_CH3_1A_LASSO_test.png', 
           img_dir + 'Ethane_298K_40Bar_CH3_1A/cm3overcm3/' + 'Ethane_298K_40Bar_CH3_1A_LASSO_test.png']
list_names = ["Ethane 4 Bar\nLASSO", 
              "Ethane 20 Bar\nLASSO", 
              "Ethane 40 Bar\nLASSO",
              "Ethane 4 Bar\nLASSO", 
              "Ethane 20 Bar\nLASSO", 
              "Ethane 40 Bar\nLASSO"]
list_labels = ["a", "b", "c", "d", "e", "f"]
imgs    = [ PIL.Image.open(i).convert("RGBA") for i in list_im ]
for a in range(0, len(imgs)):
  txt = Image.new("RGBA", imgs[a].size, (255,255,255,0))
  d = ImageDraw.Draw(txt)
  namefont = ImageFont.truetype("ariblk.ttf", 300)
  labelfont = ImageFont.truetype("ariblk.ttf", 300) # name for arial bold
  #https://stackoverflow.com/questions/1815165/draw-bold-italic-text-with-pil
  d.multiline_text((2500, 350),list_names[a],fill = (0,0,0,255),font=namefont)
  d.text((1200, 1000),list_labels[a],fill = (0,0,0,255),font=labelfont)
  
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
left = get_concat_h_cut(imgs[0], imgs[1])
up = get_concat_h_cut(left, imgs[2])
left = get_concat_h_cut(imgs[3], imgs[4])
down = get_concat_h_cut(left, imgs[5])
get_concat_v_cut(up, down).save('Figure_5.png')