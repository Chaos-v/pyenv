# -*- coding: utf-8 -*-
"""
    @ProjectName: pyenv
    @File: env.py
    @Author: Chaos
    @Date: 2023/3/30
    @Description: 对环境文件中的各个组成部分进行抽象，目前适配 bellhop
"""
from enum import Enum
from scipy.interpolate import interp1d
import numpy as np


class EnvFile:
    def __init__(self):
        pass


class Model(str, Enum):
    """
    枚举类
    """
    BELLHOP = 'bellhop'
    KRAKEN = 'kraken'
    SCOOTER = 'scooter'
    SPARC = 'sparc'


# ==================== bellhop及kraken通用部分 ====================
# (1)
class Title:
    """
    Env 文件标题部分，没有输入时直接给一个默认值

    :param sTitle: Title
    """
    def __init__(self, sTitle="\'Munk profile\'"):
        self.__title = sTitle

    # 字段访问
    def getTitle(self):
        return self.__title

    # 字段修改
    def setTitle(self, sTitle):
        self.__title = sTitle


# (2)
class Frequency:
    """
    Frequency in Hz.

    :param sFreq: 频率 (要是不写那就是50)
    """
    def __init__(self, sFreq=50):
        self.__freq = sFreq

    @property
    def freq(self):
        return self.__freq

    @freq.setter
    def freq(self, sFreq):
        self.__freq = sFreq


# (3)
class NumberOfMedia:
    """
    The problem is divided into media within which it is assumed that the material properties vary smoothly.
    A new medium should be used at fluid/elastic interfaces or at interfaces where the density changes discontinuously.
    The number of media in the problem is defined excluding the upper and lower half-space.
    BELLHOP is limited to one medium (NMedia=1) and actually ignores this parameter.

    :param aNMedia: Number of media.
    """
    def __init__(self, aNMedia):
        self.__nMedia = aNMedia

    @property
    def NMedia(self):
        return self.__nMedia

    @NMedia.setter
    def NMedia(self, aNMedia):
        self.__nMedia = aNMedia


# (4)
class TopOption:
    def __init__(self, sTopOpt):
        self.__topOpt = sTopOpt

    def checkTopOpt(self):
        pass  # 检查什么的太麻烦了，不想写了，等死吧


# (5)
class SSP:
    """
    声速剖面部分

    :param nmesh: list，包括[NMESH, SIGMA, Z(nSSP)]
    :param cZ: 声速剖面中的水深，numpy的一维向量，Z(nSSP)
    :param cCP: 声速剖面中的声速，numpy一维向量，CP(nSSP)
    :param cCS: 剪切波速度，numpy一维向量
    :param cRho: 对应深度的水密度，numpy一维向量
    :param cAP: 衰减 (p-wave)，numpy一维向量
    :param cAS: 剪切衰减，numpy一维向量
    """
    def __init__(self, nmesh, cZ, cCP, cCS=None, cRho=None, cAP=None, cAS=None):
        self.__nMesh = nmesh
        self.__z = cZ  # 深度
        self.__cp = cCP  # 声速
        self.__cs = cCS  # 剪切波速度
        self.__rho = cRho  # 密度
        self.__ap = cAP  # 衰减 (alpha (z))
        self.__as = cAS  # 剪切衰减

        self.__sspF = None
        self.betaI_f = None
        self.betaR_f = None
        self.rho_f = None
        self.alphaR_f = None
        self.alphaI_f = None

    @property
    def NMesh(self):
        return self.__nMesh

    @property
    def Depth(self):
        return self.__z

    @property
    def CP(self):
        return self.__cp

    @property
    def CS(self):
        return self.__cs

    @property
    def RHO(self):
        return self.__rho

    @property
    def AP(self):
        return self.__ap

    @property
    def AS(self):
        return self.__as

    @property
    def sspFigure(self):
        return self.__sspF

    def makeSSPf(self):
        self.__sspF = interp1d(self.Depth, self.CP)

    def interpAll(self):
        self.betaI_f = interp1d(self.Depth, self.AS)
        self.betaR_f = interp1d(self.Depth, self.CS)
        self.rho_f = interp1d(self.Depth, self.RHO)
        self.alphaR_f = interp1d(self.Depth, self.CP)
        self.alphaI_f = interp1d(self.Depth, self.AP)
        return self.alphaR_f, self.betaR_f, self.rho_f, self.alphaI_f, self.betaI_f


# (6)
class BottomOption:
    """
    Bottom Option.

    Syntax:
         BOTOPT  SIGMA

         or, if the power-law attenuation option 'm' was selected:

         BOTOPT  SIGMA BETA fT

    :param bOption: Type of bottom boundary condition.
    :param bHalfspace: Bottom Halfspace Properties from geoAcoustic values or from grain size
    :param bSigma: Interfacial roughness (m).
    :param bBeta: Power for the power law
    :param bfT: Transition frequency (Hz)
    """
    def __init__(self, bOption, bHalfspace=None, bSigma=None, bBeta=None, bfT=None):
        self.__bottomOption = bOption.upper()
        if len(bOption) == 1:
            self.__checkBotOpt1(bOption[0])
            if bOption[0].upper() == 'A' or bOption[0].upper() == 'G':
                self.__bottomHalfspace = bHalfspace
                print(type(bHalfspace))
        elif len(bOption) == 2:
            self.__checkBotOpt1(bOption[0])
            self.__checkBotOpt2(bOption[1])
        else:
            raise ValueError("Bottom Option Parameter Input Error!")

        self.__sigma = bSigma
        self.__beta = bBeta
        self.__fT = bfT

    @staticmethod
    def __checkBotOpt1(s):
        botOpt1 = ['V', 'A', 'R', 'G', 'F', 'P']
        if s[0] not in botOpt1:
            raise ValueError("BotOpt(1:1) Parameter Error!")

    @staticmethod
    def __checkBotOpt2(s):
        botOpt2 = ['~', '_', '*']
        if s not in botOpt2:
            raise ValueError("BotOpt(2:2) Parameter Error!")

    @property
    def BottomOpt(self):
        return self.__bottomOption

    @property
    def Sigma(self):
        return self.__sigma

    @property
    def Beta(self):
        return self.__beta

    @property
    def fT(self):
        return self.__fT


# (6a)
class BottomHalfspace:
    """
    Bottom Halfspace Properties from geoAcoustic values.

    :param bDepth: Depth (m)
    :param bAlphaR: Bottom P-wave speed (m/s).
    :param bBetaR: Bottom S-wave speed (m/s).
    :param bRho: Bottom density (g/cm3).
    :param bAlphaI: Bottom P-wave attenuation. (units as given by TOPOPT(3:3) )
    :param bBetaI: Bottom S-wave attenuation.
    """
    def __init__(self, bDepth, bAlphaR=np.array([]), bBetaR=np.array([]), bRho=np.array([]), bAlphaI=np.array([]),
                 bBetaI=np.array([])):
        self.__depth = bDepth
        self.__alphaR = np.array(bAlphaR)
        self.__betaR = np.array(bBetaR)
        self.__rho = np.array(bRho)
        self.__alphaI = np.array(bAlphaI)
        self.__betaI = np.array(bBetaI)

    @property
    def Depth(self):
        return self.__depth

    @property
    def AlphaR(self):
        return self.__alphaR

    @property
    def BetaR(self):
        return self.__betaR

    @property
    def Rho(self):
        return self.__rho

    @property
    def AlphaI(self):
        return self.__alphaI

    @property
    def BetaI(self):
        return self.__betaI


# (6b)
class BottomHalfspaceGrainSize:
    """
    Bottom Halfspace Properties from grain size.

    This line should only be included if BOTOPT(1:1)='G', i.e. if the user has specified a homogeneous halfspace for
    the bottom BC defined by grain size. The bottom sound speed, attenuation, and density is calculated using formulas
    from the UW_APL High-Frequency handbook.

    :param: bDepth: Depth(m)
    :param: bGrainSize: Grain size (phi units).
    """
    def __init__(self, bDepth, bGrainSize):
        self.__depth = bDepth
        self.__grainSize = bGrainSize

    @property
    def Depth(self):
        return self.__depth

    @property
    def GrainSize(self):
        return self.__grainSize


# ==================== bellhop部分 ====================


# ==================== Kraken部分 ====================


if __name__ == '__main__':
    print("========================================")
    alphaR = 1600  # p wave speed in sediment
    betaR = 0  # no shear wave
    alphaI = .5  # p wave atten
    betaI = 0  # s wave atten
    rhob = 1600

    hs = BottomHalfspace(5000, alphaR, betaR, rhob, alphaI, betaI)
    b = BottomOption('A*', hs)
    print(b.BottomOpt, hs.Depth)
