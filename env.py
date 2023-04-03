# -*- coding: utf-8 -*-
"""
    @ProjectName: pyenv
    @File: env.py
    @Author: Chaos
    @Date: 2023/3/30
    @Description: 对 title 环境文件中的各个组成部分进行抽象
"""
from enum import Enum
from scipy.interpolate import interp1d
import numpy as np


class Model(str, Enum):
    """枚举类"""
    BELLHOP = 'bellhop'
    KRAKEN = 'kraken'
    SCOOTER = 'scooter'
    SPARC = 'sparc'


class Title:
    """
    Env 文件标题部分，没有输入时直接给一个默认值
    """

    def __init__(self, sTitle="\'Munk profile\'"):
        self.__title = sTitle

    # 字段访问
    def getTitle(self):
        return self.__title

    # 字段修改
    def setTitle(self, sTitle):
        self.__title = sTitle


class Freq:
    """
    Frequency in Hz.

    :param sFreq: 频率
    """
    def __init__(self, sFreq=50):
        self.__freq = sFreq

    @property
    def freq(self):
        return self.__freq

    @freq.setter
    def freq(self, sFreq):
        self.__freq = sFreq


class SSP:
    """
    声速剖面部分

    :param sz: 声速剖面中的水深，numpy的一维向量
        Depths the ssp is taken at
    :param sAlphaR: 声速剖面中的声速，numpy一维向量
    :param sBetaR: 剪切波速度，numpy一维向量
    :param sRho: 对应深度的水密度，numpy一维向量
    :param sAlphaI: 衰减 (p-wave)，numpy一维向量
    :param sBetaI: 剪切衰减，numpy一维向量
    """
    def __init__(self, sz, sAlphaR, sBetaR, sRho, sAlphaI, sBetaI):
        self.z = sz  # 深度
        self.alphaR = sAlphaR  # 声速
        self.betaR = sBetaR  # 剪切波速度
        self.rho = sRho  # 密度
        self.alphaI = sAlphaI  # 衰减 (alpha (z))
        self.betaI = sBetaI  # 剪切衰减

        self.sspf = None
        self.betaI_f = None
        self.betaR_f = None
        self.rho_f = None
        self.alphaR_f = None
        self.alphaI_f = None

    def make_sspf(self):
        self.sspf = interp1d(self.z, self.alphaR)

    def interpAll(self):
        self.betaI_f = interp1d(self.z, self.betaI)
        self.betaR_f = interp1d(self.z, self.betaR)
        self.rho_f = interp1d(self.z, self.rho)
        self.alphaR_f = interp1d(self.z, self.alphaR)
        self.alphaI_f = interp1d(self.z, self.alphaI)
        return self.alphaR_f, self.betaR_f, self.rho_f, self.alphaI_f, self.betaI_f


class SSPPart:
    """
    还没构思好
    """
    def __init__(self):
        pass


if __name__ == '__main__':
    print("========================================")

    z1 = [0.0, 200.0, 250.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0,
          2600.0,
          2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0, 4800.0, 5000.0]

    alphaR = [1548.52, 1530.29, 1526.69, 1517.78, 1509.49, 1504.30, 1501.38, 1500.14, 1500.12, 1501.02, 1502.57,
              1504.62,
              1507.02, 1509.69, 1512.55, 1515.56, 1518.67, 1521.85, 1525.10, 1528.38, 1531.70, 1535.04, 1538.39,
              1541.76,
              1545.14, 1548.52, 1551.91]
    pw = 1
    aw = 0
    betaR = 0.0 * np.array([1] * len(z1))
    rho = pw * np.array([1] * len(z1))
    alphaI = aw * np.array([1] * len(z1))
    betaI = 0.0 * np.array([1] * len(z1))
    ssp1 = SSPraw(z1, alphaR, betaR, rho, alphaI, betaI)
