# -*- coding:utf-8 -*-
from module.core_mp import RAPP
import pickle
from glob import glob
import os
import numpy as np
import multiprocessing as mp

def main1(path, band,Rs=None):
    date = os.path.basename(path)
    folder = os.path.join('/data/output', date) # Output
    if not os.path.exists(folder):
        os.mkdir(folder)
    folder = os.path.join(folder, band)
    if not os.path.exists(folder):
        os.mkdir(folder)
    rapp = RAPP(targ=path+'/3C454/'+band,
                bias=path+'/bias/'+band,
                dark=path+'/dark/'+band,
                flat=path+'/twflat/'+band,
                expo_key='EXPOS',   # FITS header
                date_key='DATE',    # FITS header
                fp_size=(75,9),     # Blurred Star Enhancement (R_A,R_s)
                )
    rapp.info_init()
    rapp.match()
    img = rapp.img_combine()
    info = rapp.find_star(img, True, 10) # Measure the x stars
    rapp.ap(info, a=(1.0, 2.4, 3.6) )    # Aperture, Annular, Annular
    rapp.save(folder, True)
    rapp.draw(folder, True, img, info)

if __name__ == "__main__":
    #folders = glob('/data/phot/126CM/*')
    folders = ['data_example'] # Input
    for path in folders:
        for band in ['g','g','g']: 
            main1(path, band)
