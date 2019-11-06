更新:
================================
    整体修改了代码结构 现在更容易使用了 功能也更加直观与工具化
    增加了多进程功能 速度比原来更快了
    增加了合并图功能 现在可以用合并图作为参考图来测光了
    并修改了readme.txt 加入了点教程
********************************
# 提醒!
    请勿用IDE的调试运行程序
    因为程序使用了multiprocessing模块
********************************
# 简易测光教程:
>##      首先引入:

        from module import core

>##      随后创建测光对象:

        app = core.APpipeline(targ='xxx',
                              expo_key='EXPOSE',
                              date_key='DATE',)

        这里有三个必填参数 分别是fits文件夹路径 曝光时间的键 拍摄的精确时间的键
        这两个键要从fits文件中的header中找到
        创建好app对象就可以用获得以下参数:
            app.targ:   list    路径中所有fits文件的路径列表
            app.mask:   ndarr   蒙版图片 若初始化没填蒙版路径则默认为内切椭圆
            app.bias:   ndarr   bias图片 若初始化没填bias路径则为0
            app.dark:   ndarr   dark图片 若初始化没填dark路径则为0
            app.flat:   ndarr   flat图片 若初始化没填flat路径则为1
        这里建议能将可以填的参数都填上 完整版是这个样子的:
            app = ap.APpipeline(targ=,      必填
                                bias=,
                                dark=,
                                flat=,
                                mask=,
                                expo_key=,  必填
                                date_key=,  必填
                                count=,
                                N=)

>##      创建好对象后就可以进行信息初始化了:

        app.info_init()

        这一步会创建出app.info
        app.info是一个pandas.DataFrame
        结构为:(时间, (几何半径, 流量中心), 星)
    
>##      信息初始化后就可以开始匹配了:

        app.match()

        匹配后会给app对象增加一个参数app.shitfs
        app.shitfs是一个一维的复数ndarr
        每个复数代表了每张图对应参考图的偏移量
    
>##      匹配完后 就可以进行孔径测光了:

        app.ap()

        如此便会给app产生一个新的参数app.table
        app.table是一个dict
        结构为:(星, 时间, (星等, 误差))
        ap()中有两个重要的参数:
            a:      tuple   默认(1.2, 2.4, 3.6)
                (测光孔径比, 背景内孔径比, 背景外孔径比)
            gain:   float   默认1. 增益

>##      孔径测光完成后就可以保存我们需要的信息了:

        app.save(result='folder')
        app.draw(result='folder')

        save()用于保存csv文件与画出曲线图
        draw()用于画出参考图 可得知星在图片上的编号
********************************
# 合并图测光教程:
>##      首先:

        from module import core

>##      随后创建对象:

        app = core.APpipeline(data='xxx',
                              bias='xxx',
                              dark='xxx',
                              flat='xxx',
                              expo_key='EXPOSE',
                              date_key='DATE',
                              N=1,)

        因为使用图像合并 所以count使用默认值6 不需要太高 原因后面会说明

>##      创建好对象后就可以进行信息初始化了:

        app.info_init()
    
>##      信息初始化后就可以开始匹配了:

        app.match()
    
>##      这里开始就是关键 因为之后并不是开始测光 而是进行图像合并:

        img = app.img_combine()

>##      在合并图中找星:

        info = app.find_star(img=img,
                             ref=True,
                             count=10)

        这里说明一下 没有在创建对象时提供count 而在这里提供count的原因
        因为合并图的弱星肯定比原数据明显 所以在合并图中可以找到更多的星
        既然如此 在创建对象的时候就没必要找太多的星
    
>##      之后将合并图作为参考图重新匹配:

        app.match(info)

>##      重新匹配后就可以开始测光了 当然也需要将info作为参数丢给ap():

        app.ap(info)

>##      后面与简易版的区别在于draw()上:

        app.draw(result='folder',
                 img_ref=img,
                 info0=info)
        
        如此以来 参考图就会变成用合并图画出来的

>##     最后保存信息:

        app.save(result='folder')