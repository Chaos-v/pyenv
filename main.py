# -*- coding: utf-8 -*-
"""
    @ProjectName: pyenv
    @File: main.py
    @Author: Chaos
    @Date: 2023/4/3
    @Description: 重写声学工具箱matlab版本的代码以实现在 Python 中的画图功能。
                  目前根据项目需要完成声速剖面部分的绘制。
"""
from envPlot import plotssp

if __name__ == '__main__':
    print("==================== Main ====================")

    # ENV 文件路径
    envfile = 'C:\\Users\\Chaos\\Desktop\\Acoustics-Toolbox-Release_2022\\teatChaos\\Munk.env'

    envCore = plotssp(envfile)

    print('==================== 测试结束 ====================')
