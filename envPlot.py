# -*- coding: utf-8 -*-
"""
    @ProjectName: pyenv
    @File: envPlot.py
    @Author: Chaos
    @Date: 2023/3/30
    @Description: env模块相关的绘图模块，基本照着 matlab 版本写的
"""
import os
import readwrite
import numpy as np
from scipy.interpolate import interp1d, CubicSpline
import matplotlib.pyplot as plt


def plotssp(envfil):
    """
    Plots the sound speed profile\n

    Usage:
       plotssp(envfile)

    The envfile should be given without any extension

    :param: envfile: ENV file
    """
    if not envfil.endswith('.env'):
        envfil = envfil + '.env'
    # 检查文件是否存在
    if not os.path.exists(envfil):
        # 文件不存在就报FileNotFoundError，同时停止运行程序
        raise FileNotFoundError("Please check that the file name \"%s\" and path are correct! " % envfil)

    ENVCore, _ = readwrite.read_env_core(envfil)
    SSP = ENVCore.ssp
    SSPType = ENVCore.Top.Opt[0]

    # 吐槽：
    # 我写代码写的的少，看不懂这群科研巨佬这个逻辑是个什么玩意儿
    # 我就搞不懂，介质数量大于一的话，你直接拓展对应向量的维度不就行了
    # 为毛单独拉出来一个再 m 文件中单独拉出来一个 ssp.raw 的结构体
    # 搞得变量一大堆，ssp对象里面的东西也用不了
    for medium in range(ENVCore.NMedia):
        npts = round(SSP.raw[medium].z[-1] - SSP.raw[medium].z[0] + 1)

        if 'S' in SSPType or 'P' in SSPType:
            npts = 2 * npts  # need to see internal points for cubic interpolation

        z_eval = np.linspace(SSP.raw[medium].z[0], SSP.raw[medium].z[-1], npts)

        # plot the compression wave speed
        if SSPType == 'N':  # n^2 Linear
            f = interp1d(SSP.raw[medium].z, 1.0/np.square(SSP.raw[medium].cp))
            n2_eval = f(z_eval)
            c_eval = 1/np.sqrt(n2_eval)
        elif SSPType == 'P':
            raise ValueError("学艺不精，不知道用什么函数，姑且先报个错提醒一下")
        elif SSPType == 'S':  # Cubic Spline using not-a-knot boundary condition
            f = CubicSpline(SSP.raw[medium].z, SSP.raw[medium].cp)
            c_eval = f(z_eval)
        else:  # % piecewise Linear
            f = interp1d(SSP.raw[medium].z, SSP.raw[medium].cp)
            c_eval = f(z_eval)

        plt.figure(num=1, layout='constrained')
        plt.plot(c_eval.real, z_eval, color='b', linestyle='-', linewidth=2)
        plt.scatter(SSP.raw[medium].cp.real, SSP.raw[medium].z, marker='o', color='b')

        # plot the shear wave speed (if any non-zero values were supplied)
        if np.any(SSP.raw[medium].cs):
            if SSPType == 'N':
                f = interp1d(SSP.raw[medium].z, 1.0 / np.square(SSP.raw[medium].cs))
                n2_eval = f(z_eval)
                c_eval = 1 / np.sqrt(n2_eval)
            elif SSPType == 'P':
                raise ValueError("学艺不精，不知道用什么函数，姑且先报个错提醒一下")
            elif SSPType == 'S':
                f = interp1d(SSP.raw[medium].z, SSP.raw[medium].cs, kind='cubic', fill_value='extrapolate')
                c_eval = f(z_eval).real
            else:  # % piecewise Linear
                f = interp1d(SSP.raw[medium].z, SSP.raw[medium].cs)
                c_eval = f(z_eval)
            plt.close()
            plt.figure(num=1, layout='constrained')
            plt.scatter(SSP.raw[medium].cs, SSP.raw[medium].z, marker='o', color='b')
            plt.plot(c_eval.real, z_eval, color='b', linestyle='-', linewidth=2)

    plt.gca().invert_yaxis()  # because view messes up the zoom feature
    # axis IJ
    plt.xlabel('Sound Speed (m/s)')
    plt.ylabel('Depth (m)')
    plt.show()


if __name__ == '__main__':
    print("==================== 测试部分 ====================")

    envfile = 'C:\\Users\\Chaos\\Desktop\\Acoustics-Toolbox-Release_2022\\teatChaos\\Munk.env'
    envCore = plotssp(envfile)

    print('==================== 测试结束 ====================')




