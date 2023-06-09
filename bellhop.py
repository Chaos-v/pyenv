# -*- coding: utf-8 -*-
"""
    @ProjectName: main.py
    @File: bellhop.py
    @Author: Chaos
    @Date: 2023/5/4
    @Description: 
"""
from env import *
import numpy as np
from math import pi


def AddArr():
    pass


def bellhopM():
    pass


def Contents():
    pass


def crci(c, alpha, freq, AttenUnit):
    """
     Convert real wave speed and attenuation to a single complex wave speed
     5 cases:
     N for Nepers/meter
     M for dB/meter      (M for Meters)
     F for dB/m-kHZ      (F for frequency dependent)
     W for dB/wavelength (W for Wavelength)
     Q for Q
     T for Thorpe

     :param: c
     :param: alpha
     :param: freq
     :param: AttenUnit
    """
    global T, S, pH, z_bar  # 这些好像也没参与计算，也没有明显的调用

    omega = 2.0 * pi * freq

    # *** Convert to Nepers/m ***
    alphaT = 0.0

    if AttenUnit[0] == 'N':  # Nepers/m
        alphaT = alpha
    elif AttenUnit[0] == 'M':  # dB/meter
        alphaT = alpha / 8.6858896
    elif AttenUnit[0] == 'F':  # dB/m-kHZ
        alphaT = alpha * freq / 8685.8896
    elif AttenUnit[0] == 'W':  # dB/wavelength
        if c != 0.0:
            alphaT = alpha * freq / (8.6858896 * c)
    elif AttenUnit[0] == 'Q':
        if c * alpha != 0.0:
            alphaT = omega / (2.0 * c * alpha)
    else:
        print('Unknown attenuation unit')

    # added volume attenuation
    if AttenUnit[1] == 'T':  # Thorp
        f2 = np.square((freq / 1000.0))
        # Original Thorp (1967) formula
        #    alphaT = 40.0 * f2 / ( 4100.0 + f2 ) + 0.1 * f2 / ( 1.0 + f2 );
        #    alphaT = alphaT / 914.4;     # dB / m
        #    alphaT = alphaT / 8.6858896; # Nepers / m

        # Updated formula from JKPS Eq. 1.34
        Thorp = 3.3e-3 + 0.11 * f2 / (1.0 + f2) + 44.0 * f2 / (4100.0 + f2) + 3e-4 * f2  # dB/km
        Thorp = Thorp / 8685.8896  # Nepers / m
        alphaT = alphaT + Thorp
    #     *** Convert Nepers/m to equivalent imaginary sound speed ***
    alphaT = alphaT * c * c / omega
    crciO = complex(c, alphaT)
    return crciO


def InfluenceGeoGaussian():
    pass


def InfluenceGeoHat():
    pass


def makeshdarr():
    pass


def Munk_interp_tests():
    pass


def Munk_rd_axis():
    pass


def Munk():
    pass


def Munkr():
    pass


def reducestep():
    pass


def reflect():
    pass


def scalep():
    pass


def ssp_copy():
    pass


def ssp_new():
    pass


def ssp_new2():
    pass


def ssp_binterp():
    pass


def ssp_cubic():
    pass


def ssp():
    pass


def step():
    pass


def topbot(fid, freq, BCType, AttenUnit):
    """
    Handles top and bottom boundary conditions

    Input:
        ENVFIL: Environmental file
        freq:   frequency
        BCType: Boundary condition type

    Output:
        Bdry.cp:    P-wave speed in halfspace
        Bdry.cs:    S-wave speed in halfspace
        Bdry.rho:   density in halfspace
    """
    global alphaR, betaR, rhoR, alphaI, betaI

    # 不加就会报错，我也很绝望啊，可能时global没搞清楚吧
    alphaR = 1500  # defaults
    betaR = 0
    rhoR = 1
    alphaI = 0
    betaI = 0

    # *** Echo to PRTFIL user's choice of boundary condition ***
    if BCType == 'V':
        print('    VACUUM')
    elif BCType == 'R':
        print('    Perfectly RIGID')
    elif BCType == 'A':
        print('    ACOUSTO-ELASTIC half-space')
    elif BCType == 'F':
        print('    FILE used for reflection loss')
    elif BCType == 'W':
        print('    Writing an IRC file')
    elif BCType == 'P':
        print('    reading PRECALCULATED IRC')
    else:
        raise ValueError('Fatal error: Unknown boundary condition type')

    # ****** Read in BC parameters depending on particular choice ******
    cp = 0.0
    cs = 0.0
    rho = 0.0

    HS = TopHalfspace()  # Top Halfspace

    # ACOUSTO-ELASTIC half-space.
    # Requires another line with the halfspace parameters as described in env block (4a).
    # *** Half-space properties ***
    if BCType == 'A':
        tmp = next(fid)
        # tmp = [x for x in tmp if x[0].isnumeric()]  # filter out non numbers
        try:
            if not tmp.find('/') == -1:
                endIndex = tmp.index('/')
                tmp = tmp[:endIndex].split()
                tmp = [x for x in tmp if x[0].isnumeric()]  # filter out non numbers
            else:
                tmp = tmp.split()
                tmp = [x for x in tmp if x[0].isnumeric()]
        except ValueError:
            print("ENVFILE Format Error!")
        except Exception as err:
            print("Fatal Error! Error Type:\n\t" + repr(err))

        num_vals = len(tmp)
        if num_vals == 6:
            ztmp, alphaR, betaR, rhoR, alphaI, betaI = [float(x) for x in tmp[0:6]]
        elif num_vals == 5:
            ztmp, alphaR, betaR, rhoR, alphaI = [float(x) for x in tmp[0:5]]
        elif num_vals == 4:
            ztmp, alphaR, betaR, rhoR = [float(x) for x in tmp[0:4]]
        elif num_vals == 3:
            ztmp, alphaR, betaR = [float(x) for x in tmp[0:3]]
        elif num_vals == 2:
            ztmp, alphaR = [float(x) for x in tmp[0:2]]
        elif num_vals == 1:
            ztmp = [float(x) for x in tmp][0]
        else:  # there were no vals to read in so defaults will be used
            pass

        print('%10.2f    %10.2f    %10.2f    %10.2f    %10.4f    %10.4f \n' % (ztmp, alphaR, betaR, rhoR, alphaI, betaI))

        cp = crci(alphaR, alphaI, freq, AttenUnit)
        cs = crci(betaR, betaI, freq, AttenUnit)
        rho = rhoR

        HS.cp = alphaR
        HS.apt = alphaI
        HS.cs = betaR
        HS.ast = betaI
        HS.rho = rhoR

    return cp, cs, rho, HS, fid


def trace():
    pass


if __name__ == '__main__':
    print("Hello World")
