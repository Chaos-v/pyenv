# -*- coding: utf-8 -*-
"""
    @ProjectName: pyenv
    @File: env.py
    @Author: Chaos
    @Date: 2023/3/30
    @Description: 对环境文件中的各个组成部分进行抽象，目前适配 bellhop
"""
# from enum import Enum
# from scipy.interpolate import interp1d
import numpy as np


class EnvObj:
    """
    ENV 文件对象

    :param oTitle: Title 对象
    :param oFreq: Frequency 对象
    :param oNMedia: NumberOfMedia 对象
    :param oBdry: Boundary 对象，包括顶部和底部两个部分的参数
    :param oSSP: SSP 对象

    """
    def __init__(self, oTitle, oFreq, oNMedia, oBdry, oSSP):
        self.__title = oTitle
        self.__freq = oFreq
        self.__nmedia = oNMedia
        self.__topOpt = oBdry.Top
        self.__botOpt = oBdry.Bot
        self.__ssp = oSSP

    @property
    def Title(self):
        return self.__title.getTitle()

    @property
    def Freq(self):
        return self.__freq.freq

    @property
    def NMedia(self):
        return self.__nmedia.NMedia

    @property
    def Top(self):
        return self.__topOpt

    @property
    def Bot(self):
        return self.__botOpt

    @property
    def ssp(self):
        return self.__ssp


# class Model(str, Enum):
#     """
#     枚举类
#     """
#     BELLHOP = 'bellhop'
#     KRAKEN = 'kraken'
#     SCOOTER = 'scooter'
#     SPARC = 'sparc'


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
    Number of media.
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
    """
    Top Option.

    :param: aTopOpt: 长度小于7的字符串
    :param: aTopHS: 半空间参数细节
    """
    def __init__(self, aTopOpt="", aTopHS=None):
        if aTopOpt == "":
            self.__topOpt = aTopOpt
        else:
            self.__topOpt = str(aTopOpt).upper()  # 边界条件类型

        if aTopHS is None:
            self.__topHS = TopHalfspace()
        else:
            self.__topHS = aTopHS

        self.BC = None  # 我也不知道这什么逼玩意儿，函数中有就得写

    @property
    def Opt(self):
        return self.__topOpt

    @Opt.setter
    def Opt(self, aOpt):
        self.__topOpt = aOpt

    @property
    def HS(self):
        return self.__topHS

    @HS.setter
    def HS(self, oHalfspace):
        self.__topHS = oHalfspace


# (4a)
class TopHalfspace:
    """
    Top Halfspace Properties
    """
    def __init__(self):
        self.__betaI = None
        self.__alphaI = None
        self.__rho = None
        self.__betaR = None
        self.__alphaR = None
        self.__depth = None

    @property
    def depth(self):
        return self.__depth

    @depth.setter
    def depth(self, aDepth):
        self.__depth = aDepth

    @property
    def cp(self):
        return self.__alphaR

    @cp.setter
    def cp(self, aAlphaR):
        self.__alphaR = aAlphaR

    @property
    def cs(self):
        return self.__betaR

    @cs.setter
    def cs(self, aBetaR):
        self.__betaR = aBetaR

    @property
    def rho(self):
        return self.__rho

    @rho.setter
    def rho(self, aRho):
        self.__rho = aRho

    @property
    def apt(self):
        return self.__alphaI

    @apt.setter
    def apt(self, aAlphaI):
        self.__alphaI = aAlphaI

    @property
    def ast(self):
        return self.__betaI

    @ast.setter
    def ast(self, aBetaI):
        self.__betaI = aBetaI


# (4b)
class BioLayerParam:
    """
    Biological Layer Parameters (for attenuation due to fish)

    不用看了，就是没写，谁需要谁写
    """
    pass


# (5)
class SoundSpeedProfile:
    """
    声速剖面部分，包括读取 ENV 文件的核心以及按照介质数量划分数据的 SSPRaw 对象

    :param cNmesh: list，包括[NMESH, SIGMA, Z(nSSP)]
    :param cZ: 声速剖面中的水深，numpy的一维向量，Z(nSSP)
    :param cCP: 声速剖面中的声速，numpy一维向量，CP(nSSP)
    :param cCS: 剪切波速度，numpy一维向量
    :param cRho: 对应深度的水密度，numpy一维向量
    :param cAP: 衰减 (p-wave)，numpy一维向量
    :param cAS: 剪切衰减，numpy一维向量
    """
    def __init__(self, cNmesh=None, cSigma=None, cDepth=None, cZ=None, cCP=None, cCS=None, cRho=None, cAP=None, cAS=None):
        # ssp部分第一行三个参数
        self.__nMesh = cNmesh
        self.__sigma = cSigma
        self.__depth = cDepth
        # ssp部分剖面部分
        self.__z = cZ  # 深度
        self.__cp = cCP  # 声速
        self.__cs = cCS  # 剪切波速度
        self.__rho = cRho  # 密度
        self.__ap = cAP  # 衰减 (alpha (z))
        self.__as = cAS  # 剪切衰减

        self.__sspRaw = []  # 包含 SSPRaw 对象的列表

        # The Francois-Garrison formula depends on salinity (S), temperature (T), pH, and depth (z_bar).
        # That information is then provided on the line immediately following
        # 为什么不私有化了：懒，还有就是写烦了
        self.francoisGarrisonParam = {}  # 仅当 TOPOPT(4:4)=F 时


    @property
    def NMESH(self):
        return self.__nMesh

    @property
    def DEPTH(self):
        return self.__depth

    @property
    def SIGMA(self):
        return self.__sigma

    @property
    def Z(self):
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
    def raw(self):
        return self.__sspRaw

    @raw.setter
    def raw(self, nList):
        self.__sspRaw = nList

    # @property
    # def sspFigure(self):
    #     return self.__sspF

    # def makeSSPf(self):
    #     self.__sspF = interp1d(self.DEPTH, self.CP)
    #
    # def interpAll(self):
    #     self.betaI_f = interp1d(self.Depth, self.AS)
    #     self.betaR_f = interp1d(self.Depth, self.CS)
    #     self.rho_f = interp1d(self.Depth, self.RHO)
    #     self.alphaR_f = interp1d(self.Depth, self.CP)
    #     self.alphaI_f = interp1d(self.Depth, self.AP)
    #     return self.alphaR_f, self.betaR_f, self.rho_f, self.alphaI_f, self.betaI_f


class SSPRaw:
    """
    每一个 SSPRaw 对象存储一个 NMedia 的 SSP 数据，包括：z, cp, cs, rho, ap, as

    :param z: list Depth，之后转化为一维 numpy 矩阵
    :param alphaR: list SSP value，之后转化为一维 numpy 矩阵
    :param betaR: list shear speeds，之后转化为一维 numpy 矩阵
    :param rho: list density value at the points，之后转化为一维 numpy 矩阵
    :param alphaI: list attenuation (p-wave)，之后转化为一维 numpy 矩阵
    :param betaI: list shear attenuation，之后转化为一维 numpy 矩阵
    """
    def __init__(self, z, alphaR, betaR, rho, alphaI, betaI):
        self.__z = np.array(z)
        self.__cp = np.array(alphaR)  # cp
        self.__cs = np.array(betaR)  # cs
        self.__rho = np.array(rho)  # rho
        self.__ap = np.array(alphaI)  # ap
        self.__as = np.array(betaI)  # as

    @property
    def z(self):
        return self.__z

    @property
    def cp(self):
        return self.__cp

    @property
    def cs(self):
        return self.__cs

    @property
    def rho(self):
        return self.__rho

    @property
    def AP(self):
        return self.__ap

    @property
    def AS(self):
        return self.__as


# (6)
class BottomOption:
    """
    Bottom Option.

    Syntax:
         BOTOPT  SIGMA

         or, if the power-law attenuation option 'm' was selected:

         BOTOPT  SIGMA BETA fT

    :param bOption: String(length 1 or 2). Type of bottom boundary condition.('V''A''R''G''F')
    :param bSigma: Interfacial roughness (m).
    :param bHalfspace: BottomHalfspace('A') or BottomHalfspaceGrainSize('G') Object.
                       Bottom Halfspace Properties from geoAcoustic values or from grain size，
    :param bBeta: Power for the power law
    :param bfT: Transition frequency (Hz)
    """
    def __init__(self, bOption, bSigma=None, bHalfspace=None, bBeta=None, bfT=None):
        if bOption is None:
            self.__bottomOption = bOption
            if bHalfspace is None:
                # 如果没有指定半空间细节，默认创建一个对象
                self.__bottomHalfspace = BottomHalfspace()

        else:
            self.__bottomOption = str(bOption).upper()  # 边界条件类型
            self.__bottomHalfspace = self.__createHalfspaceObject(self.__checkBotOpt(bOption), bHalfspace)

        self.__sigma = bSigma
        self.__beta = bBeta
        self.__fT = bfT

    @staticmethod
    def __checkBotOpt(s):
        """
        检查参数是否存在

        :param s: BotOpt string
        :return: 1 BottomHalfspace 对象，对应参数 'A'
        :return: 2 BottomHalfspaceGrainSize 对象，对应参数 'G'
        :return: 0 什么都没有
        """
        botOpt1 = ['V', 'A', 'R', 'G', 'F', 'P']
        botOpt2 = ['~', '_', '*', ' ']
        if 0 < len(s) <= 2:
            if len(s) == 2 and s[1] not in botOpt2:
                raise ValueError("BotOpt(2:2) Parameter Error!")
            if s[0] not in botOpt1:
                raise ValueError("BotOpt(1:1) Parameter Error!")
            else:
                if s[0] == 'A':
                    return 1
                elif s[0] == 'G':
                    return 2
                else:
                    return 0
        else:
            raise ValueError("Bottom Option Parameter Error!")

    @staticmethod
    def __createHalfspaceObject(n, sHS):
        """根据输入值创建不同的Halfspace对象"""
        if n == 0:
            return sHS
        elif n == 1:
            return BottomHalfspace()
        elif n == 2:
            BottomHalfspaceGrainSize()
        else:
            raise ValueError("Parameter Error!")

    @property
    def Opt(self):
        return self.__bottomOption

    @Opt.setter
    def Opt(self, aBotOpt: str):
        # 检查底部选项字符串的长度，然后检查根据情况对半空间赋值
        sOpt = aBotOpt.upper()
        self.__checkBotOpt(sOpt)
        self.__bottomOption = sOpt

    @property
    def HS(self):
        """
        底部声学半空间，
        """
        return self.__bottomHalfspace

    @HS.setter
    def HS(self, botHS):
        self.__bottomHalfspace = botHS

    @property
    def Sigma(self):
        return self.__sigma

    @Sigma.setter
    def Sigma(self, n):
        self.__sigma = n


    @property
    def Beta(self):
        return self.__beta

    @Beta.setter
    def Beta(self, n):
        self.__beta = n

    @property
    def fT(self):
        return self.__fT

    @fT.setter
    def fT(self, n):
        self.__fT = n


# (6a)
class BottomHalfspace:
    def __init__(self, bDepth=None, bAlphaR=np.array([]), bBetaR=np.array([]), bRho=np.array([]), bAlphaI=np.array([]),
                 bBetaI=np.array([])):
        """
            Bottom Halfspace Properties from geoAcoustic values.

            :param bDepth: Depth (m)
            :param bAlphaR: Bottom P-wave speed (m/s).
            :param bBetaR: Bottom S-wave speed (m/s).
            :param bRho: Bottom density (g/cm3).
            :param bAlphaI: Bottom P-wave attenuation. (units as given by TOPOPT(3:3) )
            :param bBetaI: Bottom S-wave attenuation.
        """
        self.__depth = bDepth
        self.__alphaR = np.array(bAlphaR)
        self.__betaR = np.array(bBetaR)
        self.__rho = np.array(bRho)
        self.__alphaI = np.array(bAlphaI)
        self.__betaI = np.array(bBetaI)

    @property
    def depth(self):
        return self.__depth

    @depth.setter
    def depth(self, aDepth):
        self.__depth = aDepth

    @property
    def cp(self):
        return self.__alphaR

    @cp.setter
    def cp(self, aAlphaR):
        self.__alphaR = aAlphaR

    @property
    def cs(self):
        return self.__betaR

    @cs.setter
    def cs(self, aBetaR):
        self.__betaR = aBetaR

    @property
    def rho(self):
        return self.__rho

    @rho.setter
    def rho(self, aRho):
        self.__rho = aRho

    @property
    def apb(self):
        return self.__alphaI

    @apb.setter
    def apb(self, aAlphaI):
        self.__alphaI = aAlphaI

    @property
    def asb(self):
        return self.__betaI

    @asb.setter
    def asb(self, aBetaI):
        self.__betaI = aBetaI


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
    def __init__(self, bDepth=None, bGrainSize=None):
        self.__depth = bDepth
        self.__grainSize = bGrainSize

    @property
    def Depth(self):
        return self.__depth

    @Depth.setter
    def Depth(self, nDepth):
        self.__depth = nDepth

    @property
    def GrainSize(self):
        return self.__grainSize

    @GrainSize.setter
    def GrainSize(self, nGrainSize):
        self.__grainSize = nGrainSize


# 无分类
class Boundary:
    """
    边界部分，包括顶部边界和底部边界两个对象
    """
    def __init__(self):
        self.Top = TopOption()
        self.Bot = BottomOption(None)


# ==================== bellhop部分 ====================


# ==================== Kraken部分 ====================


if __name__ == '__main__':
    print("========================================")

