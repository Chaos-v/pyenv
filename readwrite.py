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
from bellhop import *
from re import findall


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

    Bdry = Boundary()
    Bdry.Top.HS.cp = 0.0
    Bdry.Top.HS.cs = 0.0
    Bdry.Top.HS.rho = 0.0
    Bdry.Bot.HS.cp = 2000.0
    Bdry.Bot.HS.cs = 0.0
    Bdry.Bot.HS.rho = 2.0

    try:
        with open(envfil, 'r') as f:
            lines = f.readlines()
            fid = iter(lines)
    except FileNotFoundError:
        print('Unable to open environmental file')
    except Exception as err:
        print("Fatal Error! Error Type:\n\t"+repr(err))

    # ========== (1) ==========
    # Extract letters between the quotes
    TitleEnv = findall('\'(.*)\'', next(fid))[0]
    print(TitleEnv)
    TITLE = Title(TitleEnv)  # ENV文件标题，=========>指向 class Title

    # ========== (2) ==========
    freq = float((next(fid).split())[0])
    print('Frequency = %d Hz ' % freq)
    FREQ = Frequency(freq)  # 频率，=========>指向 class Frequency

    # ========== (3) ==========
    NMedia = int(next(fid).split()[0])
    print('Number of media = %i \n' % NMedia)
    NMEDIA = NumberOfMedia(NMedia)  # Number of media. =========>指向 class NumberOfMedia

    # ========== (4) ==========
    # Extract option letters between the quotes
    TopOpt = findall('\'(.*)\'', next(fid))[0]
    TopOpt = TopOpt + ' '*(7 - (len(TopOpt)+1))
    # convert the deprecated '*' option to '~'
    if TopOpt[4] == '*':
        TopOpt = TopOpt[: 4] + '~' + TopOpt[5:]
    Bdry.Top.Opt = TopOpt

    SSPType = Bdry.Top.Opt[0]
    Bdry.Top.BC = Bdry.Top.Opt[1]  # 别问，问就不知道这个BC是什么逼玩意儿
    AttenUnit = Bdry.Top.Opt[2:4]

    # *** SSP approximation options ***
    if SSPType == 'N':
        print('    N2-Linear approximation to SSP')
    elif SSPType == 'C':
        print('    C-Linear approximation to SSP')
    elif SSPType == 'P':
        print('    PCHIP approximation to SSP')
    elif SSPType == 'S':
        print('    Spline approximation to SSP')
    elif SSPType == 'Q':
        print('    Quadrilateral approximation to range-dependent SSP')
    elif SSPType == 'H':
        print('    Hexahedral approximation to range and depth dependent SSP')
    elif SSPType == 'A':
        print('    Analytic SSP option')
    else:
        raise ValueError('Fatal error: Unknown option for SSP approximation')

    # *** Attenuation options ***
    attenUnit = AttenUnit[0]
    if attenUnit == 'N':
        print('    Attenuation units: nepers/m')
    elif attenUnit == 'F':
        print('    Attenuation units: dB/mkHz')
    elif attenUnit == 'M':
        print('    Attenuation units: dB/m')
    elif attenUnit == 'W':
        print('    Attenuation units: dB/wavelength')
    elif attenUnit == 'Q':
        print('    Attenuation units: Q')
    elif attenUnit == 'L':
        print('    Attenuation units: Loss tangent')
    else:
        raise ValueError('Fatal error: Unknown attenuation units')

    # *** optional addition of volume attenuation using standard formulas
    if len(Bdry.Top.Opt) >= 4:
        if Bdry.Top.Opt[3] == 'T':
            print('    THORP attenuation added')
        elif Bdry.Top.Opt[3] == 'F':
            print('    Francois-Garrison attenuation added')
            # 创建字典，=========>指向 class SSP.francoisGarrisonParam
            fGParam = {"Temperature": float(next(fid)), "Salinity": next(fid), "pH": next(fid),
                       "z_bar": next(fid)}

            T = fGParam.get('Temperature')
            S = fGParam.get('Salinity')
            pH = fGParam.get('pH')
            z_bar = fGParam.get('z_bar')

            print('        T =  %4.1f degrees   S = %4.1f psu   pH = %4.1f   z_bar = %6.1f m \n' % (T, S, pH, z_bar))

    if len(Bdry.Top.Opt) >=5:
        if Bdry.Top.Opt[4] == '*':
            print('    Development options enabled')

    Bdry.Top.cp, Bdry.Top.cs, Bdry.Top.rho, Bdry.Top.HS, fid = topbot(fid, freq, Bdry.Top.BC, AttenUnit)

    # main loop to readin in SSP
    print('       z          alphaR         betaR           rho        alphaI         betaI');
    print('      (m)          (m/s)         (m/s)         (g/cm^3)      (m/s)         (m/s) ');

    sspDict = {}
    sspDict['z'] = []
    sspDict['c'] = []
    sspDict['cs'] = []
    sspDict['rho'] = []

    sspDict['N'] = [0] * NMedia
    sspDict['sigma'] = [0] * NMedia
    sspDict['depth'] = [0] * (NMedia + 1)
    # ==========> 指向 class SSP.NMesh
    sspDict['raw'] = []
    sspDict['cz'] = []
    sspDict['Npts'] = [0] * NMedia

    Loc = [0] * NMedia

    for medium in range(NMedia):
        if medium == 0:
            Loc[medium] = 0
        else:
            Loc[medium] = Loc[medium - 1] + SSP.Npts(medium - 1)

        # NMesh 部分，这部分只要读的时候有标点符号就会报错，不管了
        tmp = next(fid).split()
        sspDict['N'][medium], sspDict['sigma'][medium], sspDict['depth'][medium+1] = int(float(tmp[0])), float(tmp[1]), float(tmp[2])

        print('    ( Number of points = %d  Roughness = %6.2f  Depth = %8.2f )' % (sspDict['N'][medium], sspDict['sigma'][medium], sspDict['depth'][medium+1]))

        # read in the SSP

    print("function end")







if __name__ == '__main__':
    print("========================================")
    envfile = 'C:\\Users\\Chaos\\Desktop\\Acoustics-Toolbox-Release_2022\\teatChaos\\Munk.env'
    read_env_core(envfile)
