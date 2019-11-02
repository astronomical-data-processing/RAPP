# -*- coding:utf-8 -*-
# pylint:disable=maybe-no-member
from module import core


import csv
import multiprocessing as mp
import os
from glob import glob

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from matplotlib import cm
from matplotlib import pyplot as plt
from scipy import ndimage as nd


def run1():
    app = core.APpipeline(data=r'/media/xiezhou/_dde_data/data/rawdata/3C454/i',
                          #   bias=r'/media/xiezhou/_dde_data/data/rawdata/bias/g',
                          #   dark=r'/media/xiezhou/_dde_data/data/rawdata/dark/g',
                          #   flat=r'/media/xiezhou/_dde_data/data/rawdata/twflat/g',
                          expo_key='EXPOS',
                          date_key='DATE',)
    app.info_init()
    app.match()
    img = app.img_combine()
    info0 = app.find_star(img, True, 7)
    app.match(info0)
    app.ap(info0)
    app.draw(folder=folder,
             img_ref=img,
             info0=info0)
    app.save(folder=folder)


def run2():
    app = core.APpipeline(data=r'/media/xiezhou/_dde_data/data/rawdata/3C454/i',
                          bias=r'/media/xiezhou/_dde_data/data/rawdata/bias/i',
                          dark=r'/media/xiezhou/_dde_data/data/rawdata/dark/i',
                          flat=r'/media/xiezhou/_dde_data/data/rawdata/twflat/i',
                          expo_key='EXPOS',
                          date_key='DATE',)
    app.info_init()
    app.match()
    img = app.img_combine()
    info0 = app.find_star(img, True, 7)
    app.match(info0)
    app.ap(info0)
    app.draw(folder=folder,
             img_ref=img,
             info0=info0)
    app.save(folder=folder)


if __name__ == "__main__":
    folder = 'local/result1/'
    run1()
    folder = 'local/result2/'
    run2()
