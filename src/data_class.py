import sys
import os
import pdb
import numpy as np
import tecplot as tp
import matplotlib
matplotlib.use("TkAgg")            ####### added for MAC
from matplotlib import pyplot as plt





class Data(object):
    # -----------------------------------------------------------------------------------------------#
    # INITIALIZE THE CLASS WITH A TECPLOT FILE
    # -----------------------------------------------------------------------------------------------#
    def __init__(self, IN):
        "Automatic Initialization Function: Takes in the read file"
        self.IN         = IN
        self.fileName   = self.IN['fileName']

        self.DATA = tp.data.load_tecplot(self.fileName + ".plt")
        print('Loading data from tecplot file: ', self.fileName + ".plt. The following are the available variable.")
        print(self.DATA.variable_names)

        self.InputsFromFile()
        #self.test()





    # -------------------------------------------------------------------------------------------------------#
    #       CLASS VARIABLES INITIALIZATION
    # -------------------------------------------------------------------------------------------------------#

    def InputsFromFile(self):

        ########################################################################
        #       GEOMETRICAL VARIABLES
        # cartesian coordinates
        self._x             = self.DATA.variable('X').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._y             = self.DATA.variable('Y').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._z             = self.DATA.variable('Z').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        #polar coordinates with respect to Z
        self._radZ          = self.DATA.variable('radZ').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._thetaZ        = self.DATA.variable('thetaZ').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._spanZ         = self.DATA.variable('span_z').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()

        self._radIn         = max(self._radZ)
        self._radOut        = min(self._radZ)

        self._radDistri     = np.linspace(max(self._radZ), min(self._radZ), 50)


        print("Inflow radius: ", self._radIn)
        print("Outflow radius: ", self._radOut)

        ########################################################################
        #       THERMODYNAMIC VARIABLES
        self._rho           = self.DATA.variable('RHO').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._T             = self.DATA.variable('temp').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._h             = self.DATA.variable('enthalpy').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._h_tot         = self.DATA.variable('hTotal').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()

        self._p             = self.DATA.variable('press').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._p_tot         = self.DATA.variable('pTotal').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._p_aver        = self.DATA.variable('aver_press').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()

        self._s             = self.DATA.variable('entropy').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._s_aver        = self.DATA.variable('aver_entropy').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._s_gen         = self.DATA.variable('sgen').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()

        ########################################################################
        #       DYNAMIC VARIABLES
        # absolute velocity, cartesian coordinates
        self._velX          = self.DATA.variable('vel-X').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._velY          = self.DATA.variable('vel-Y').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._velZ          = self.DATA.variable('vel-Z').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        # relative velocity, cartesian coordinates
        self._vrelX         = self.DATA.variable('vrel-X').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._vrelY         = self.DATA.variable('vrel-Y').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._vrelZ         = self.DATA.variable('vrel-Z').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        # mach numbers
        self._Mach          = self.DATA.variable('Mach').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._MachRel       = self.DATA.variable('MachRel').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._Mach_aver     = self.DATA.variable('aver_mach').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._MachRel_aver  = self.DATA.variable('aver_machRel').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()

        ########################################################################
        #       TURBULENCE MODELING RELATED VARIABLES (Spallart-Allmaras)
        self._nuSA          = self.DATA.variable('nuSA').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()
        self._muT           = self.DATA.variable('muT').values('writeFlaggedCvsTecplot, np=14').as_numpy_array()


        self._n             = len(self._x)
        self._index         = np.arange(1., 1., self._n)

        print('the size of the array is: ', self._n)

    # -------------------------------------------------------------------------------------------------------#
    #       test function
    # -------------------------------------------------------------------------------------------------------#

    def test(self):


        x = np.arange(0., 5.0, 0.2)
        plt.figure(1)
        #plt.plot(x, x, 'bs', x, x**2, 'r--', x, x**3, 'g^')
        plt.show()

    def GetMassAverage(xx, yy, prop, angle):
        LINE = tp.data.extract.extract_line(zip(xx, yy))
        a = abs(yy[0] - yy[-1]) / (len(yy) - 1)
        PRESSURE = (LINE.values('Pressure').as_numpy_array())
        DENSITY = (LINE.values('Density').as_numpy_array())
        ENTROPY = (LINE.values(prop).as_numpy_array())
        MOMENTUM_X = (LINE.values('X-Momentum').as_numpy_array()) * np.sin(np.radians(angle))
        MOMENTUM_Y = (LINE.values('Y-Momentum').as_numpy_array()) * np.cos(np.radians(angle))
        MOMENTUM = (MOMENTUM_X ** 2 + MOMENTUM_Y ** 2) ** 0.5
        M_ENTROPY = np.sum(MOMENTUM * ENTROPY * a)
        M_MASS = np.sum(MOMENTUM * a)
        # print(M_ENTROPY)
        # print(M_MASS)
        return (M_ENTROPY / M_MASS)




