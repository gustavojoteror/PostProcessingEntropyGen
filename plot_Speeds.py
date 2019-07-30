import os
import sys
from sys import argv
import logging
import csv
import numpy as np
import matplotlib
#matplotlib.use("TkAgg")            ####### added for MAC
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams['figure.figsize'] = [8,10]

from src.IO import ReadUserInput
from src.function import readPointCart
from src.function import IntegrateProperty
from src.point import Point



if __name__ == "__main__":
    """Plot function to compared Q3D and 3D e.g. 'python plot_numDoamin.py XXXXX' where XXX is the case to be plotted"""

    typeNumDom  = argv[1]
    propName   = argv[2]
    my_dpi = 96


    DIR         = os.getcwd()+'/'
    INFile      = DIR+'param.in'
    IN          = ReadUserInput(INFile)

    typeAVG     = IN['typeAVG']
    if typeNumDom=='Q3D':
        speeds      = ['26krpm', '28krpm', '32krpm'] #['26krpm', '28krpm', '32krpm']
    else:
        speeds      = ['26krpm', '28krpm','32krpm'] #] #['26krpm', '28krpm', '32krpm']
    plotformat  = ['b--', 'k-', 'r-.', 'k-']
    scale       = 1
    savefig     = 1
    Intgrate    = 1

    colorVertLines= 'grey'
    formatVertLines='--'

    Pt_sle     = readPointCart(IN, 'STATOR_LE')
    Pt_thr     = readPointCart(IN, 'THROAT_A')
    Pt_ste     = readPointCart(IN, 'STATOR_TE')
    Pt_intL    = readPointCart(IN, 'INTERFACE_L')

    radSLE  = Pt_sle.rad
    radthr  = Pt_thr.rad+0.0005
    radSTE  = Pt_ste.rad
    radintL = Pt_intL.rad
    radRLE  =  float(IN['radLERotor'])
    radRTE  =  float(IN['radTERotor'])

    lineSLE  = np.ones(2)*radSLE
    lineSTE  = np.ones(2)*radSTE
    linethr  = np.ones(2)*radthr
    lineintL = np.ones(2)*radintL
    lineRLE  = np.ones(2)*radRLE
    lineRTE  = np.ones(2)*radRTE
    entropy  = np.ones(2)

    if typeAVG=='Vol':
        if propName=='s_t':
            entropy[0]= 0
            entropy[1]= 1200
        elif propName=='s_v' or propName=='s_g':
            if Intgrate == 1:
                entropy[0]= 0
                entropy[1]= 25000
            else:
                entropy[0]= 0
                entropy[1]= 1600000
        elif propName=='Be':
            entropy[0]= 0
            entropy[1]= 0.25
    elif typeAVG=='Mass':
        if propName=='s_t':
            entropy[0]= 0
            entropy[1]= 50
        elif propName=='s_v' or propName=='s_g':
            if Intgrate == 1:
                entropy[0]= 0
                entropy[1]= 1000
            else:
                entropy[0]= 0
                entropy[1]= 3800000 #1200000
        elif propName=='Be':
            entropy[0]= 0
            entropy[1]= 0.25

    value=1.0


    for i in range(len(speeds)):

        if typeNumDom=='Q3D':
            file = "DATA/{0}/entropy{1}_{2}.txt".format(typeNumDom, typeAVG, speeds[i])
        else:
            file = "DATA/{0}/entropy{1}_{2}.txt".format(typeNumDom, typeAVG, speeds[i])
        print('reading the data from: ', file)

        property=[]
        radius=[]

        radius, property= IntegrateProperty(file, propName,Intgrate)

        if scale == 1:
            radius= radius/radintL


        ax = plt.subplot(1, 1, 1)
        plt.plot(radius, property,  plotformat[i], linewidth=2, markersize=6)


    if scale == 1:
        lineSLE  *= 1/radintL
        linethr  *= 1/radintL
        lineSTE  *= 1/radintL
        lineintL *= 1/radintL
        lineRLE  *= 1/radintL
        lineRTE  *= 1/radintL

    # plotting constant radius line
    textLocationY = 0.75
    dlocationY = 0.03
    dlocationX = -0.0

    plt.plot(lineSLE, entropy,  formatVertLines, color=colorVertLines, linewidth=1.5)
    # plt.text(lineSLE[0]+dlocationX, entropy[0]+((entropy[1]-entropy[0])*textLocationY), "SLE")
    plt.text(lineSLE[0]+dlocationX, entropy[0]+((entropy[1]-entropy[0])*0.4), "SLE", rotation=-90, fontsize=14)
    plt.plot(linethr, entropy, formatVertLines, color=colorVertLines, linewidth=1.5)
    plt.text(linethr[0]+dlocationX, entropy[0]+((entropy[1]-entropy[0])*0.5), "Throat", rotation=-90, fontsize=14)
    plt.plot(lineSTE, entropy, formatVertLines, color=colorVertLines, linewidth=1.5)
    plt.text(lineSTE[0]+dlocationX, entropy[0]+((entropy[1]-entropy[0])*0.6), "STE", rotation=-90, fontsize=14)
    plt.plot(lineintL, entropy, formatVertLines, color=colorVertLines, linewidth=1.5)
    plt.text(lineintL[0]+dlocationX, entropy[0]+((entropy[1]-entropy[0])*0.7), "S/R", rotation=-90, fontsize=14)
    plt.plot(lineRLE, entropy, formatVertLines, color=colorVertLines, linewidth=1.5)
    plt.text(lineRLE[0]+dlocationX, entropy[0]+((entropy[1]-entropy[0])*0.8), "RLE", rotation=-90, fontsize=14)
    plt.plot(lineRTE, entropy, formatVertLines, color=colorVertLines, linewidth=1.5)
    plt.text(lineRTE[0]+dlocationX, entropy[0]+((entropy[1]-entropy[0])*0.9), "RTE", rotation=-90, fontsize=14)

    # formatting
    plt.xlim(min(radius),radSLE/radintL+0.05)
    plt.ylim(entropy[0],entropy[1])


    plt.xlabel('Scaled radius, R/R$_{S/R} [-]$', fontsize=16)
    if typeAVG=='Vol':
        lableee = "Accumulated entropy generation (volume averaged), ${0} [W/m^3 K]$".format(propName)
        plt.ylabel(lableee, fontsize=16)
    elif typeAVG=='Mass':
        lableee = "Accumulated rate of entropy generation (mass averaged), ${0} [W/m^3 K]$".format(propName)
        plt.ylabel(lableee, fontsize=16)

    if propName=='Be':
        lableee = "Bejan number, ${0} [-]$".format(propName)
        plt.ylabel(lableee, fontsize=16)

    # ax.set_yscale('log')
    if typeNumDom=='Q3D':
        plt.legend((speeds[0], speeds[1], speeds[2]))
    else:
        plt.legend((speeds[0], speeds[1], speeds[2]))
    # ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    if savefig==1:
        file = "DATA/{0}/{1}_{2}.png".format(typeNumDom, propName, typeAVG)
        plt.savefig(file)
        print('imaged stored: ', file)

    plt.show()