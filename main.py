# -*- coding:utf-8 -*-
from module.core import RAPP


def run():
    rapp = RAPP(targ=r'data/',
                # bias=r'bias/',  # 可选
                # dark=r'dark/',  # 可选
                # flat=r'flat/',  # 可选
                expo_key='EXPOS',
                date_key='DATE',
                # mask=,  # 可选
                # count=,  # 可选
                # N=,
                )
    rapp.info_init()
    rapp.match()
    img = rapp.img_combine()
    info0 = rapp.find_star(raw=img,
                           ref=True,
                           count=10)
    rapp.match(info0=info0)
    rapp.ap(info0=info0,
            # a=,  # 可选
            # gain=,  # 可选
            )
    rapp.draw(folder=folder,
              # show_all=,  # 可选
              ref=img,
              info0=info0
              # a=,  # 可选
              # ont_size=,  # 可选
              )
    rapp.save(folder=folder)


if __name__ == "__main__":
    folder = 'result/'
    run()
