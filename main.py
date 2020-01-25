# -*- coding:utf-8 -*-
from module.core_mp import RAPP


def run():
    rapp = RAPP(targ=path+'/rawdata/3C454/'+band,
                bias=path+'/rawdata/bias/'+band,
                dark=path+'/rawdata/dark/'+band,
                flat=path+'/rawdata/twflat/'+band,
                expo_key='EXPOS',
                date_key='DATE',
                # fp_size=(75, 10),
                )
    rapp.info_init()
    rapp.match()
    img = rapp.img_combine()
    info = rapp.find_star(img, True, 10)
    rapp.ap(info)
    rapp.save(folder)
    rapp.draw(folder, True, img, info)


if __name__ == "__main__":
    folder = 'result/'
    run()
