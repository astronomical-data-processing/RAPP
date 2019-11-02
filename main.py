# -*- coding:utf-8 -*-
from module import core


def run():
    app = core.APpipeline(data=r'data/',
                          bias=r'bias/',
                          dark=r'dark/',
                          flat=r'flat/',
                          expo_key='EXPOS',
                          date_key='DATE',)
    app.info_init()
    app.match()
    img = app.img_combine()
    info0 = app.find_star(img, True, 10)
    app.match(info0)
    app.ap(info0)
    app.draw(folder=folder,
             img_ref=img,
             info0=info0)
    app.save(folder=folder)


if __name__ == "__main__":
    folder = 'local/result1/'
    run()
