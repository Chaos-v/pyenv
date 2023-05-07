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


def plotssp(envfil):
    """
    Plots the sound speed profile\n

    Usage:
       plotssp(envfile)

    The envfile should be given without any extension

    :param: envfile: ENV file
    :return: none
    """
    if not envfil.endswith('.env'):
        envfil = envfil + '.env'
    # 检查文件是否存在
    if not os.path.exists(envfil):
        # 文件不存在就报FileNotFoundError，同时停止运行程序
        raise FileNotFoundError("Please check that the file name \"%s\" and path are correct! " % envfil)

    print("========== Program end. Debug End ==========")


if __name__ == '__main__':
    print("========================================")
    envfile = 'C:\\Users\\Chaos\\Desktop\\Acoustics-Toolbox-Release_2022\\teatChaos\\Munk.env'
    envCore, fid = readwrite.read_env_core(envfile)


