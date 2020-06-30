# -*- coding:utf-8 -*-
# Author:weirong, zwj, Jnk_xz,
# pylint:disable=maybe-no-member
# pylint:disable=invalid-unary-operand-type
# Robust Automated Photometry Pipeline for Blurred Images
# https://ui.adsabs.harvard.edu/abs/2020PASP..132g5001H/abstract
import csv
import multiprocessing as mp
import os
from glob import glob

import numpy as np
import pandas as pd
from astropy import convolution, io
from astropy.time import Time
from matplotlib import cm
from matplotlib import pyplot as plt
from scipy import ndimage as nd


def progress_bar(i, l, title=''):
    '''
    在终端上显示进度条
    para:
        i:      int 当前循环数
        l:      int 总循环数
        title:  str 显示在进度条前的字符
    '''
    i += 1
    columns = (+ os.get_terminal_size().columns
               - len(title)
               - len(str(i))
               - len(str(l))
               - 10)
    if columns < 2:
        print(title + '[%d/%d %d%%]' % (
            i,
            l,
            int(i/l * 100)),
            end='\r')
    left = np.around(i/l * columns).astype(int)
    right = columns - left
    print(title + '[%d/%d %d%%][%s%s]' % (
        i, l, int(i/l * 100), '6' * left, 'o' * right
    ), end='\r')
    if i == l:
        print()


def remove_outliers(raw):
    '''
    异常值去除函数
    首先将图片与中值滤波后的中值图相减再取绝对值 获得残差图res
    随后以标准差运算函数作为卷积核 对残差图进行卷积 获得标准差图
    将残差图与N倍标准差图对比 就能获取异常值的位置信息
    将图片中的异常值位置 用同样位置的中值替换 并返回
    就做到了异常值去除的功能
    Para:
        img:    ndarray  数据图片
        t:      int      卷积核面积(t>=3 & t为奇数)
    return:
        img:    ndarray  去除异常值后的数据图片 datatype为float
    '''
    img = raw.copy()
    # 假如用3*3的权重 会奇异的进入一个3倍std的情况
    # 正好无法去除异常值 所以这里选用这个特别的窗口
    weights = np.array([[1., 1., 1.],
                        [1., 0., 1.],
                        [1., 1., 1.]])/8.
    med = nd.median_filter(img, 3)
    res = img - med
    std = np.sqrt(nd.convolve(res**2, weights))
    mask = (np.abs(res) > 3*std)
    img[mask] = med[mask]
    return img


def collect_data(paths, desc, expo_key=None):
    '''
    收集路径中的数据与曝光时间 并转成float ndarray返回
    '''
    expose = []
    datas = []
    lenth = len(paths)
    for i, path in enumerate(paths):
        progress_bar(i, lenth, desc)
        with io.fits.open(path, ignore_missing_end=True) as f:
            img = f[0].data
            if expo_key is not None:
                expose.append(float(f[0].header[expo_key]))
        img = img.astype(float)
        datas.append(img)
    if expo_key is None:
        return datas
    else:
        return expose, datas


def bginfo(img, mask=-1, threshold=0.05):
    '''
    自创去信号的背景信息算法
    返回 背景的sky与std
    '''
    if mask is not -1:
        img = img[mask]
    std_old = np.std(img)
    sky_old = np.median(img)
    while 1:
        jud = (img < sky_old + 3*std_old) & (img > sky_old - 3*std_old)
        std = np.std(img[jud])
        sky = np.median(img[jud])
        if np.abs(std - std_old)/std_old < threshold:
            break
        std_old = std
        sky_old = sky
    return sky, std


def img_scale(raw, top=3, bot=1):
    '''
    画图用的压缩 用于提高显示对比度
    '''
    img = raw.copy()
    sky, std = bginfo(img)
    img[img > sky+top*std] = sky+top*std
    img[img < sky-bot*std] = sky-bot*std
    return img


def circle(img, center, radius):
    '''
    画圆, 返回圆内的值的数组
    '''
    h, w = img.shape
    Y, X = center.real, center.imag
    # 这里可能有大优化 就是注释的内容 可以替换所有以下的内容
    y = np.arange(int(np.around(Y - radius)),
                  int(np.around(Y + radius)))
    y = np.intersect1d(np.arange(h), y)
    x = np.arange(int(np.around(X - radius)),
                  int(np.around(X + radius)))
    x = np.intersect1d(np.arange(w), x)
    x, y = np.meshgrid(x, y)
    mask = (((Y - y)**2 + (X - x)**2)**0.5 <= radius)
    values = img[y[mask], x[mask]]
    return values


def annular(img, center, inner, outer):
    '''
    画环, 返回环内的值的数组
    '''
    values = []
    h, w = img.shape
    Y, X = center.real, center.imag
    y = np.arange(int(np.around(Y - outer)),
                  int(np.around(Y + outer)))
    y = np.intersect1d(np.arange(h), y)
    x = np.arange(int(np.around(X - outer)),
                  int(np.around(X + outer)))
    x = np.intersect1d(np.arange(w), x)
    x, y = np.meshgrid(x, y)
    radius = ((Y - y)**2 + (X - x)**2)**0.5
    mask = (inner <= radius) * (radius <= outer)
    values = img[y[mask], x[mask]]
    return values


class RAPP(object):
    def __init__(self, targ, expo_key, date_key, N=3, mask: np.ndarray = True, fp_size=(75, 9), **kwarg):
        '''
        初始化各种路径APpipline 
        并生成bias, dark, flat, mask
        Para:
            targ:       path    数据路径
            expo_key:   str     曝光时间的键
            mask:       ndarr   (可选)蒙版路径
            data_key:   str     (可选)曝光时间关键词
            N:          int     (可选)获取多少倍背景标准差的信号信息 默认3
            fp_size:    int     (可选)背景环内径大小 与 中值滤波半径 默认(75, 9)
        kwarg:
            bias:       path    (可选)本底路径
            dark:       path    (可选)暗场路径
            flat:       path    (可选)平场路径
        说明:
            当bias为默认状态时, 则bias为0
            当dark为默认状态时, 则dark为0
            当flat为默认状态时, 则flat为1
            当mask为默认状态时, 则自动生成一个内切椭圆的mask
        '''
        self.targ = glob(os.path.join(targ, '*.fit*'))

        with io.fits.open(self.targ[0], ignore_missing_end=True) as f:
            img = f[0].data
            header = f[0].header
        targExp = header[expo_key]
        self.W, self.H = img.shape
        self.date_key = date_key
        self.N = N
        self.font_size = 24
        self.outliers = True
        bias_path = kwarg['bias'] if 'bias' in kwarg else ''
        dark_path = kwarg['dark'] if 'dark' in kwarg else ''
        flat_path = kwarg['flat'] if 'flat' in kwarg else ''

        biasp = glob(os.path.join(bias_path, '*.fit*'))
        darkp = glob(os.path.join(dark_path, '*.fit*'))
        flatp = glob(os.path.join(flat_path, '*.fit*'))

        if mask is True:
            y, x = np.meshgrid(np.arange(self.H),
                               np.arange(self.W))
            y, x = y - self.H/2, x - self.W/2
            a, b = self.H/2, self.W/2
            self.mask = (y**2 / a**2 + x**2 / b**2 <= 1)
        else:
            self.mask = mask

        self.bias = 0.0
        if biasp != []:
            self.bias = np.median(collect_data(biasp, 'bias:'), axis=0)

        self.dark = 0.0
        if darkp != []:
            self.dark = np.median(collect_data(darkp, 'dark:'), axis=0)
            self.dark -= self.bias

        self.flat = 1.0
        if flatp != []:
            flatExp, f = collect_data(flatp, 'flat:', expo_key=expo_key)
            for i in range(len(f)):
                flat = f[i]
                flat -= self.bias + self.dark / targExp * flatExp[i]
                f[i] = flat/np.median(flat)
            self.flat = np.median(f, axis=0)

        self.fpa = convolution.Ring2DKernel(fp_size[0], 1)
        self.fps = convolution.Tophat2DKernel(fp_size[1])

        self.fpa.normalize()
        self.fpa.normalize()

        arrs = np.array(self.fps).astype(bool)

        Ns = len(arrs[arrs])

        self.kernel = np.array(self.fps * Ns - self.fpa * Ns)

    def load(self, path, loop=False):
        '''
        读取路径(path)的fits文件图片 并进行预处理与
        para:
            path:       str     路径
            loop:       bool    是否用于循环 是:返回增加拍摄时间与文件名信息 否:只返回图片 默认:False
        return:
            img:        ndarr   图像
            jd:         float   时间
            name:       str     无前后缀文件名
        '''
        if loop:
            name = os.path.basename(path).split('.')[0]
        with io.fits.open(path, ignore_missing_end=True) as f:
            img = f[0].data
            if loop:
                jd = Time(f[0].header[self.date_key], format='fits').jd
        img = img.astype(float)
        img = (img - self.bias - self.dark) / self.flat
        if self.outliers:
            img = remove_outliers(img)
        if loop:
            return img, jd, name
        return img

    def find_star(self, raw, ref=False, count=8, showLBL=False, name=''):
        '''
        找星程序
        首先对图片进行中值滤波 获取滤波图img_med
        将滤波图中 滤波图 < [滤波图中值 + N * 原图标准差] 的位置定义为背景
        余下连通区进行标记 根据大小进行排序 最大的count个连通区定义为星
        随后记录星的几何半径与流量中心 并返回
        para:
            img:        ndarr       图片
            ref:        bool        True:返回pandas.DataFrame对象 False:返回dict对象
            count:      int         找多少颗星 默认:0 假如是0时, 则取对象中的count参数
        return:
            dic:        dict/df     返回的找星数据 有两个关键词
                radius:     几何半径
                centers:    流量中心
        '''
        if ref:
            # 叠加图背景去除
            img = raw - nd.median_filter(raw, footprint=np.array(self.fpa))
        else:
            # 信号增强
            raw = nd.median_filter(raw, 5)
            img = nd.convolve(raw, self.kernel)
        # 计算背景标准差
        sky, std = bginfo(img, self.mask)
        # 标记背景
        mark = img > sky + self.N*std
        # # 清除跨蒙版星
        # lbl, _ = nd.measurements.label(mark | ~self.mask)
        # mark = mark & (lbl != 1)
        # 计算连通区
        lbl, num = nd.measurements.label(mark)
        idx = np.arange(num) + 1
        # 根据连通区计算半径
        r_arr = nd.labeled_comprehension(input=lbl,
                                         labels=lbl,
                                         index=idx,
                                         func=lambda x: np.sqrt(len(x)/np.pi),
                                         out_dtype=float,
                                         default=0)
        # 第一次筛选连通区
        # jud = r_arr > 9
        # idx = idx[jud]
        # r_arr = r_arr[jud]
        if ref:
            jud = r_arr > 3
            idx = idx[jud]
            r_arr = r_arr[jud]
            # 根据R排序
            jud = np.argsort(-r_arr)
            idx = idx[jud]
            r_arr = r_arr[jud]
            if len(idx) > count:
                idx = idx[:count]
                r_arr = r_arr[:count]
        else:
            # jud = r_arr > 6
            # idx = idx[jud]
            # r_arr = r_arr[jud]
            # 获取SNR列表
            F = nd.sum(img, lbl, idx)  # OR MAX
            # # 第二次筛选连通区
            # jud = Fmax > np.pi * 9**2 * 0
            # idx = idx[jud]
            # r_arr = r_arr[jud]
            # Fmax = Fmax[jud]
            # 根据SNR最大值排序
            jud = np.argsort(-F)
            idx = idx[jud]
            r_arr = r_arr[jud]
            F = F[jud]
            # 第三次筛选连通区
            if len(idx) > 8:
                idx = idx[:8]
                r_arr = r_arr[:8]
                F = F[:8]
        if showLBL:
            save = np.zeros_like(lbl).astype(int)
            for i in idx:
                save += (lbl == i) * 1
            if not os.path.exists('local/lbl/'):
                os.mkdir('local/lbl/')
            plt.imsave('local/lbl/'+name+'_lbl.png', save)
            plt.close()
        # 计算质心
        centers = nd.measurements.center_of_mass(input=raw,
                                                 labels=lbl,
                                                 index=idx)
        centers = [complex(*center) for center in centers]
        if ref:
            return pd.DataFrame({'radius': pd.Series(r_arr), 'centers': pd.Series(centers)})
        return {'radius': pd.Series(r_arr), 'centers': pd.Series(centers)}

    def info_mp(self, path):
        '''
        多进程调用函数
        '''
        img, jd, name = self.load(path, True)
        return self.find_star(img, name=name), (jd, name)

    def info_init(self):
        '''
        初始化信号信息 创建self.info
        随后可以调用match()
        '''
        dic = {}
        pool = mp.Pool(mp.cpu_count()-2)
        i = 0
        progress_bar(i, len(self.targ), 'data:')
        for sub_dic, key in pool.imap(self.info_mp, self.targ):
            progress_bar(i, len(self.targ), 'data:')
            i += 1
            dic[key] = sub_dic
        pool.close()
        pool.join()
        self.info = pd.DataFrame(dic)

    def old_match(self, threahold=2, info0=None):
        '''
        匹配
        配准思路为 将所有可能的平移量进行尝试
        因为图片中的星必然能对齐 所以正确的平移量必然可以在其他星配准时正确
        依照这点将正确的平移量脱颖而出 并搜集
        Para:
            info0:      dataframe   参考图的信息 默认为None None时则选择当前最早拍摄的图片作为参考图
        '''
        self.shifts = []
        if info0 is None:
            info0 = self.info[self.info.columns.values[np.argmin(
                self.info.columns.codes)]]
        c_arr0 = info0['centers'].to_numpy()
        lenth = len(self.info.columns)
        for i, key in enumerate(self.info):
            progress_bar(i, lenth, 'match:')
            c_arr = self.info[key]['centers'].to_numpy()
            # 做差 得到所有可能的平移量
            shifts = c_arr[:, np.newaxis] - c_arr0
            # 再次做差 得到所有可能的平移量的残差
            res = shifts[:, :, np.newaxis, np.newaxis] - shifts
            # 将残差量化
            res = np.linalg.norm([res.real, res.imag], axis=0)
            # 设定阈值 给差打分
            points = (res < threahold) * 1
            # 搜集分数
            scores = np.sum(points, axis=2)
            scores = np.sum(scores, axis=2)
            # 因为必然有多个正确的平移量 所以大部分情况下 正确的平移量都是高分
            idxs = np.argwhere(scores >= np.max(scores)-1)
            shifts = shifts[idxs[:, 0], idxs[:, 1]]
            shift = np.average(shifts, axis=0)
            self.shifts.append(shift)

    def similar(self, D1, D2, threahold):
        res = np.abs(D2 - D1[:, np.newaxis])
        # 距离差低于阈值
        simi_arr = np.sum((res < threahold) * 1, axis=1)
        return len(simi_arr[simi_arr > 0])

    def distance(self, C):
        D = []
        for c1 in C:
            row = []
            for c2 in C:
                if c1 == c2:
                    continue
                d = c2 - c1
                row.append(np.linalg.norm([d.real, d.imag]))
            row = np.array(row)
            # 最近的6颗星
            # row = np.sort(row)[:] if len(row) >= 6 else row
            D.append(row)
        return np.array(D)

    def match(self, threahold=2, info0=None):
        self.shifts = []
        self.flag = []
        self.idx0 = []
        self.idxi = []
        if info0 is None:
            info0 = self.info[self.info.columns.values[np.argmin(
                self.info.columns.codes)]]
        C0 = info0['centers'].to_numpy()
        D0 = self.distance(C0)
        l = len(self.info.columns)
        for i, key in enumerate(self.info):
            progress_bar(i, l, 'match:')
            Ci = self.info[key]['centers'].to_numpy()
            Di = self.distance(Ci)
            Si = []
            for d0 in D0:
                row = []
                for di in Di:
                    row.append(self.similar(d0, di, threahold))
                Si.append(np.array(row))
            Si = np.array(Si)
            Smax = np.max(Si)
            if Smax < 3:
                self.shifts.append(0 + 0*1j)
                self.flag.append(False)
                self.idx0.append(None)
                self.idxi.append(None)
            else:
                idx0, idxi = np.argwhere(Si == np.max(Si))[0]
                self.shifts.append(Ci[idxi] - C0[idx0])
                self.flag.append(True)
                idx0, idxi = [], []
                for i, si in enumerate(Si):
                    if np.max(si) < 3:
                        continue
                    idx0.append(i)
                    idxi.append(np.argmax(si))
                self.idx0.append(np.array(idx0))
                self.idxi.append(np.array(idxi))

    def ap(self, info0=None, a=(1.2, 2.4, 3.6), gain=1.):
        '''
        孔径测光 创建self.table
        随后可运行draw(), save()
        Para:
            info0:      dataframe       参考图的信息 默认为None None时则选择当前最早拍摄的图片作为参考图
            a:          tuple           默认(1.2, 2.4, 3.6) (测光孔径比, 背景内孔径比, 背景外孔径比)
            gain:       float           默认1. 增益
        '''
        self.aperture, self.inner, self.outer = a
        if info0 is None:
            info0 = self.info[
                self.info.columns.values[
                    np.argmin(
                        self.info.columns.codes)]]

        c_arr0 = info0['centers'].to_numpy()
        r_arr0 = info0['radius'].to_numpy()

        r = max(r_arr0)  # 有没有更好的方法来获取这个值呢?

        table = {
            str(c0): {
                jd: [float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN')] for jd, _ in self.info
            } for c0 in c_arr0
        }

        lenth = len(self.info.columns)
        for i, items in enumerate(zip(self.targ, self.shifts, self.flag)):
            progress_bar(i, lenth, 'ap:')
            path, shift, flag = items
            if not flag:
                continue
            c_arr = c_arr0 + shift
            img, jd, _ = self.load(path, True)
            for i, items in enumerate(zip(c_arr, c_arr0)):
                # 孔径测光部分
                c, c0 = items
                adu = circle(img, c, r*self.aperture)
                if len(adu) == 0:
                    continue
                skys = annular(img, c, r*self.inner, r*self.outer)
                sky, std = bginfo(skys)
                flux = sum(adu - sky)
                if flux <= 0:
                    continue
                mag = -2.5 * np.log10(flux)
                err = 2.5 / np.log(10) / flux * np.sqrt(
                    + flux / gain
                    + len(adu) * std**2
                    + len(adu)**2 * std**2 / len(skys))
                if flux / (np.sqrt(len(adu)) * std) < 5:
                    continue
                table[str(c0)][jd] = [mag, err, flux, len(adu), std]
        self.table = table

    def draw_circle(self, filename, img, c_arr, r, flag):
        '''
        画孔径函数 参数解释同darw()
        '''
        img = img_scale(img)
        w = 10.8/2
        h = w / self.W * self.H
        fig, ax = plt.subplots(figsize=(h, w))
        ax.axis('off')
        ax.imshow(img, cmap=cm.cividis)
        if flag:
            for i, c in enumerate(c_arr):
                y, x = c.real, c.imag
                font = {
                    'family': 'serif',
                    'color': 'yellow',
                    'weight': 'normal',
                    'size': self.font_size
                }
                ax.text(x=x+r*2,
                        y=y+r*2,
                        s=i,
                        fontdict=font)
                ax.add_artist(plt.Circle((x, y),
                                         r*self.aperture,
                                         edgecolor='red',
                                         facecolor=(0, 0, 0, .0125)))
                ax.add_artist(plt.Circle((x, y),
                                         r*self.inner,
                                         edgecolor='green',
                                         facecolor=(0, 0, 0, .0125)))
                ax.add_artist(plt.Circle((x, y),
                                         r*self.outer,
                                         edgecolor='green',
                                         facecolor=(0, 0, 0, .0125)))
        plt.tight_layout()
        fig.savefig(filename + '.png')
        plt.close()

    def draw(self, folder='result', show_all=False, ref=None, info0=None):
        '''
        画图函数 将星图画上孔径进行保存 用以确定每颗星的编号
        para:
            folder:     str     保存文件夹 文件默认文件名为ref.png
            show_all:   bool    是否显示所有图片 默认:False True:显示所有
            img_ref:    ndarr   参考图路径 默认:None None时则将最早拍摄的图片作为参考图
            info0:      df      参考图对应的参考图信息 默认:None
        '''
        if not os.path.exists(folder):
            os.makedirs(folder)
        filename = os.path.join(folder, 'ref')
        if ref is not None and info0 is None:
            info0 = self.find_star(ref, True)
        elif ref is None and info0 is None:
            idx0 = np.argmin(self.info.columns.codes)
            key0 = self.info.columns.values[idx0]
            info0 = self.info[key0]
            _, name = key0
            ref = self.load(self.targ[idx0])
        c_arr0 = info0['centers'].to_numpy()
        r_arr0 = info0['radius'].to_numpy()
        centers = np.array(self.shifts)[:, np.newaxis] + c_arr0
        r = max(r_arr0)
        self.draw_circle(filename=filename,
                         img=ref,
                         c_arr=c_arr0,
                         r=r,
                         flag=True)
        if show_all:
            for path, c_arr, flag in zip(self.targ, centers, self.flag):
                img, _, name = self.load(path, True)
                filename = os.path.join(folder, name)
                self.draw_circle(filename=filename,
                                 img=img,
                                 c_arr=c_arr,
                                 r=r,
                                 flag=flag)

    def save(self, folder='result', getD=False):
        '''
        将数据保存为csv与生成光变曲线
        Para:
            folder:    str  默认为相对路径中的result 结果路径
        '''
        if not os.path.exists(folder):
            os.makedirs(folder)
        print('outputing plot and csv')
        jds = self.info.columns.levels[0]
        w = 10.8
        h = w / self.W * self.H
        plt.figure(figsize=(h, w))
        for i, fc in enumerate(self.table):
            name = str(i).zfill(np.log10(len(self.table)).astype(int) + 1)
            filename = os.path.join(folder, name)
            with open(filename + '.csv', 'w', newline='') as f:
                writer = csv.writer(f)
                mag_lst = []
                err_lst = []
                for jd in jds:
                    mag = self.table[str(fc)][jd][0]
                    err = self.table[str(fc)][jd][1]
                    flux = self.table[str(fc)][jd][2]
                    number = self.table[str(fc)][jd][3]
                    std = self.table[str(fc)][jd][4]
                    snr = flux / (np.sqrt(number) * std)
                    writer.writerow([jd, mag, err, snr, flux, number, std])
                    mag_lst.append(mag)
                    err_lst.append(err)
                # if np.isnan(np.sum(mag_lst)):
                #     continue
                # plt.errorbar(jds, mag_lst, err_lst, label=i)
                plt.plot(jds, mag_lst, label=i)
        plt.ylim(-20, -10)
        plt.legend()
        filename = os.path.join(folder, 'plot')
        plt.savefig(filename + '.eps')
        plt.close()
        if getD:
            idx0_in = self.idx0[0].copy()
            for i in range(len(self.idx0)):
                idx0 = self.idx0[i]
                if self.flag[i]:
                    idx0_in = np.intersect1d(idx0_in, idx0)
            for i in range(len(idx0_in)):
                D = []
                for j, key in enumerate(self.info):
                    if not self.flag[j]:
                        continue
                    Ci = self.info[key]['centers'].to_numpy()
                    idx0 = self.idx0[j]
                    idxi = self.idxi[j]
                    _, _, jud_i = np.intersect1d(
                        idx0_in, idx0, return_indices=True)
                    idxi = idxi[jud_i]
                    Di = Ci[idxi] - Ci[idxi[i]]
                    norm = np.linalg.norm([Di.real, Di.imag], axis=0)
                    norm = np.sort(norm)
                    D.append(norm)
                D = np.array(D)
                D_std = np.std(D, axis=0)
                filename = os.path.join(folder, str(i).zfill(2)+'_D.txt')
                np.savetxt(filename, D_std, delimiter=', ', fmt='%.4f')

    def img_combine_big(self):
        '''
        合并大图像 将所有图都合并到一张大图中 返回这张图
        与img_combine的区别在于 img_combine_big是含有所有信息的
        img_combine是仅截取与原图同样像素的图片
        但big仅可以用来看 没法用来作为参考图
        '''
        shifts = np.array(self.shifts)
        ys_arr, xs_arr = shifts.real, shifts.imag
        ymin, ymax = np.min(ys_arr), np.max(ys_arr)
        xmin, xmax = np.min(xs_arr), np.max(xs_arr)
        h, w = self.W, self.H
        img_total = np.zeros((h + int(abs(ymin)) + int(abs(ymax)),
                              w + int(abs(xmin)) + int(abs(xmax))))
        iys_arr, ixs_arr = ys_arr.astype(int), xs_arr.astype(int)
        L = len(iys_arr)
        Y = int(abs(ymax))
        X = int(abs(xmax))
        for i, items in enumerate(zip(iys_arr, ixs_arr, self.targ)):
            progress_bar(i, L, 'combine:')
            iys, ixs, path = items
            img = self.load(path)
            sky, std = bginfo(img, mask=self.mask)
            img -= sky
            img /= std
            y = Y - iys
            x = X - ixs
            img_total[y:y+h, x:x+w] += img * self.mask
        img_total /= np.sqrt(L)
        return img_total

    def combine(self, path):
        '''
        合并图的多进程调用函数
        '''
        if not self.flag[self.targ.index(path)]:
            return False, None, None, None
        shift = self.shifts[self.targ.index(path)]
        img = self.load(path)
        sky, std = bginfo(img, mask=self.mask)
        img -= sky
        img /= std
        h, w = img.shape

        Y, X = shift.real, shift.imag

        yT = np.arange(h) - int(Y)
        y = np.arange(h) + int(Y)

        xT = np.arange(w) - int(X)
        x = np.arange(w) + int(X)

        xT = np.intersect1d(np.arange(w), xT)
        x = np.intersect1d(np.arange(w), x)

        yT = np.intersect1d(np.arange(h), yT)
        y = np.intersect1d(np.arange(h), y)

        xT, yT = np.meshgrid(xT, yT)
        x, y = np.meshgrid(x, y)

        return True, yT, xT, img[y, x] * self.mask[y, x]

    def img_combine(self):
        '''
        合并所有图 将所有图都合并到一张图中 返回这张图
        合并图可作为参考图进行匹配与测光 效果应该会更好
        '''
        img_total = np.zeros((self.W, self.H))
        lenth = len(self.targ)
        pool = mp.Pool(mp.cpu_count()-2)
        i = 0
        N = 0
        for res in pool.imap(self.combine, self.targ):
            progress_bar(i, lenth, 'combine:')
            i += 1
            flag, y, x, img = res
            if flag:
                N += 1
                img_total[y, x] += img
        pool.close()
        pool.join()
        return img_total / np.sqrt(N)
