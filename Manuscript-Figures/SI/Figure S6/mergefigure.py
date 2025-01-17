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
master_dir = file_dir + '/../../../Results' # gets master directory
img_dir = master_dir + '/Xe_Kr_SingleComponent_Analysis/'
list_im = [img_dir + 'Kr1bar_Loading_vs_LCD.png', 
           img_dir + 'Kr10bar_Loading_vs_LCD.png', 
           img_dir + 'Xe1bar_Loading_vs_LCD.png', 
           img_dir + 'Xe10bar_Loading_vs_LCD.png']
list_names = ["Kr 273 K, 1 Bar", "Kr 273 K, 10 Bar", "Xe 273 K, 1 Bar", "Xe 273 K, 10 Bar"]
list_labels = ["a", "b", "c", "d"]
imgs    = [ PIL.Image.open(i).convert("RGBA") for i in list_im ]
for a in range(0, len(imgs)):
  txt = Image.new("RGBA", imgs[a].size, (255,255,255,0))
  d = ImageDraw.Draw(txt)
  namefont = ImageFont.truetype("ariblk.ttf", 300)
  labelfont = ImageFont.truetype("ariblk.ttf", 500) # name for arial bold
  #https://stackoverflow.com/questions/1815165/draw-bold-italic-text-with-pil
  d.text((2000, 250),list_names[a],fill = (0,0,0,255),font=namefont)
  d.text((2500, 500),list_labels[a],fill = (0,0,0,255),font=labelfont)
  
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
right = get_concat_h_cut(imgs[2], imgs[3])

get_concat_v_cut(left, right).save('Figure_S6.png')