# -*- coding: utf-8 -*-
"""
    @ProjectName: pyenv
    @File: readwrite.py
    @Author: Chaos
    @Date: 2023/3/30
    @Description: 将给定的参数写入 ENV 文件中
        ！！！ 一般而言，给定 env 后读取的情况大多余直接给定参数写入 env 文件，该功能有待完善。
"""
from env import *
def read_env_core(envfil):
    """
    Read the core of the environmental file
    This is the part used by all the models in the Acoustics Toolbox
    """
    global HV, NFirstAcoustic, NLastAcoustic
    global alphaR, betaR, rhoR, alphaI, betaI

    alphaR = 1500  # defaults
    betaR = 0
    rhoR = 1
    alphaI = 0
    betaI = 0

    NFirstAcoustic = 0

    Top = TopOption()
    Bot = BottomOption()




    return TitleEnv, freq, ssp, bdry, lines, line_ind





if __name__ == '__main__':
    print("Hello World")
