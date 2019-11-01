更新:
===========
    整体修改了代码结构 现在更容易使用了 功能也更加直观与工具化
    增加了多进程功能 速度比原来更快了
    增加了合并图功能 现在可以用合并图作为参考图来测光了
    并修改了readme.txt 加入了点教程

# 简易测光教程:
    ##      首先引入:

        import APpipeline as ap

##      随后创建对象:

        app = ap.APpipeline(data='xxx',
                            expkey='EXPOS',)

###     这里有两个必填参数 分别是fits文件路径与曝光时间的键
###     曝光时间的键要从fits文件中的header中找到
###     创建好app对象就可以用获得以下参数:
            app.datap:  list    路径中所有fits文件的路径列标
            app.mask:   ndarr   蒙版图片 若初始化没填蒙版路径则默认为内切椭圆
            app.bias:   ndarr   bias图片 若初始化没填bias路径则为0
            app.dark:   ndarr   dark图片 若初始化没填dark路径则为0
            app.flat:   ndarr   flat图片 若初始化没填flat路径则为1
###     这里建议能将可以填的参数都填上 完整版是这个样子的:
            app = ap.APpipeline(data=,
                                bias=,
                                dark=,
                                flat=,
                                mask=,
                                expkey=,
                                data_key=,
                                count_star=,
                                N=)

##      创建好对象后就可以进行信息初始化了:

        app.info_init()

###     这一步会创建出app.info
###     app.info是一个pandas.DataFrame
###     一共有三个轴 分别是(时间, (几何半径, 流量中心), 星)
###     其中信息分为两个分别是
    
##      信息初始化后就可以开始匹配了:

        app.match()

###     匹配后会给app对象增加一个参数app.shitfs
###     app.shitfs是一个一维的复数ndarr
###     每个复数代表了每张图对应参考图的偏移量
    
##      匹配完后 就可以进行孔径测光了:

        app.ap()

###     如此便会给app产生一个新的参数app.table
###     app.table是一个字典结构为:(星, 时间, (星等, 误差))
###     每个流量中心实际上代表了一颗星

##      孔径测光完成后就可以保存我们需要的信息了:

        app.save(result='file')
        app.draw(result='file')

###     app.save()用于保存csv文件与画出曲线图
###     app.draw()用于画出参考图 可得知星在图片上的编号

# 合并图测光教程:
##      首先:

        import APpipeline as ap

##      随后创建对象:

        app = ap.APpipeline(data='xxx',
                            bias='xxx',
                            dark='xxx',
                            flat='xxx',
                            mask='xx_mask.png',
                            expkey='EXPOS',
                            data_key='DATE',
                            N=1,)

###     因为使用图像合并 所以count_star使用默认值 不需要太高

##      创建好对象后就可以进行信息初始化了:

        app.info_init()
    
##      信息初始化后就可以开始匹配了:

        app.match()
    
##      这里开始就是关键 因为之后并不是开始测光 而是进行图像合并:

        img = app.img_combine()
    
##      这样就得到了合并图 这时可以调用matplotlab来显示一下这张图片
##      但这就先不这么做 直接开始对合并图进行找星:

        info = app.find_star(img=img,
                             ref=True,
                             count_star=10)

###     说明一下: 这里的count_star可以根据显示出的合并图来数一下需要多少星
###     find_star()根据亮到暗找出count_star数量的星的信息并返回成info
###     info是pandas.DataFrame类型的 之后需要用到
    
##      之后将合并图作为参考图重新匹配:

        app.match(info)

##      重新匹配后就可以开始测光了 当然也需要将info作为参数丢给ap():

        app.ap(info)

##      后面与简易版的区别在于draw()上:

        app.draw(result='file',
                 img_ref=img,
                 info0=info)
        
###     如此以来 参考图就会变成用合并图画出来的
