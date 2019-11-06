# -*- coding:utf-8 -*-
from module import core


def run():
    app = core.APpipeline(targ=r'data/',
                          # bias=r'bias/',  # 可选
                          # dark=r'dark/',  # 可选
                          # flat=r'flat/',  # 可选
                          expo_key='EXPOS',
                          date_key='DATE',
                          # mask=,  # 可选
                          # count=,  # 可选
                          # N=,
                          )
    app.info_init()
    app.match()
    img = app.img_combine()
    info0 = app.find_star(img=img,
                          ref=True,
                          count=10)
    app.match(info0=info0)
    app.ap(info0=info0,
           # a=,  # 可选
           # gain=,  # 可选
           )
    app.draw(folder=folder,
             # show_all=,  # 可选
             img_ref=img,
             info0=info0
             # a=,  # 可选
             # ont_size=,  # 可选
             )
    app.save(folder=folder)


if __name__ == "__main__":
    folder = 'local/result1/'
    run()
