import os
import sys
import logging
import csv
import numpy as np
import tecplot as tp
import copy

from src.point import Point

#=====================================================
#               functions
#=====================================================
def IntegrateProperty(file, propName,integrate):
    step=[]
    radius=[]
    s_v=[]
    s_t=[]
    s_gen=[]
    dVA=[]
    property=[]


    list_curve=[]
    tmp = []

    f=open(file,"r")
    lines=f.readlines()

    for x in lines:
        step.append(float(x.split(' ')[0]))
        radius.append(float(x.split(' ')[1]))
        s_v.append(float(x.split(' ')[2]))
        s_t.append(float(x.split(' ')[3]))
        s_gen.append(float(x.split(' ')[4]))
        dVA.append(float(x.split(' ')[5]))
        tmp.append(x)
    f.close()



    if propName=='s_t':
        property =  s_t
    elif propName=='s_v':
        property =  s_v
    elif propName=='s_g':
        property =  s_gen
    elif propName=='Be':
        property =  copy.deepcopy(s_t)



    if integrate == 1:
        if propName=='Be':
            # for j in range(len(property)-1):
            #     s_t[j+1] = (s_t[j+1]*dVA[j+1] + s_t[j]*dVA[j]) /(dVA[j+1] + dVA[j])
            #     s_v[j+1] = (s_v[j+1]*dVA[j+1] + s_v[j]*dVA[j]) /(dVA[j+1] + dVA[j])
            #     dVA[j+1] = dVA[j+1] + dVA[j]
            for j in range(len(property)):
                property[j] =s_t[j]/(s_t[j]+s_v[j])
        else:
            for j in range(len(property)-1):
                property[j+1] = (property[j+1]*dVA[j+1] + property[j]*dVA[j]) /(dVA[j+1] + dVA[j])
                dVA[j+1] = dVA[j+1] + dVA[j]
    else:
        if propName=='Be':
            for j in range(len(property)):
                property[j] =s_t[j]/(s_t[j]+s_v[j])
        else:
            for j in range(len(property)-1):
                property[j+1] = property[j+1] + property[j]

    return radius , property

def readPointCart(IN, name):
    ###############################
    # reading a single point from the input file
    # input:   IN   (input file)
    #          name (name of the point)
    # output:  point

    x0, y0, z0 = IN[name].split(",")
    return Point(x=float(x0), y=float(y0), z=float(z0))


def readPointsFile(name, type='XY', delimeter="\t"):
    ###############################
    # reading an array of points from a file
    # input:   name     (name of the file with the points)
    #          type     (how the points are defined with: XYZ, ZY, RZ (radius and z-coordinate)
    # output:  pts      (array of points)

    pts =[]
    f = open(name, "r")
    lines = f.readlines()
    for x in lines:
        if type=='XYZ':
            tmp =  Point(x=float(x.split(delimeter)[0]), y=float(x.split(delimeter)[1]), z=float(x.split(delimeter)[2]))
        elif type=='RZ':
            tmp =  Point(rad=float(x.split(delimeter)[0]), z=float(x.split(delimeter)[1]), theta=float(x.split(delimeter)[2]))
        elif type=='RTHETA':
            tmp =  Point(rad=float(x.split(delimeter)[0]), z=float(x.split(delimeter)[2]), theta=float(x.split(delimeter)[1]))
        elif type == 'XY':
            tmp = Point(x=float(x.split(delimeter)[0]), y=float(x.split(delimeter)[1]))
        else:
            print("In readPointsFile, %s is an invalid type... (XYZ, RZ, XY) Exiting... " % (type))
            sys.exit()
        # pts= np.append(pts, tmp)
        pts.append(tmp)
    f.close()

    return pts


def readHeightDistribution(IN, typeNumDom):
    ###############################
    # reading the height distributions
    # input:   IN           (input file)
    #          typeNumDom   (how the domain is defined: 3D, Q3D
    # output:  heightDist   (height distribution)
    ###############################################################################
    # Reading

    fileheight  = IN['fileHD']
    filemidspan = IN['fileMS']

    hd = readPointsFile(fileheight, 'RZ')
    ms = readPointsFile(filemidspan, 'RZ', ' ')  # extra input is because it has a different delimiter (space)

    # scaling points as they are in [mm]
    for i in range(len(hd)):
        hd[i].scale(1 / 1000)
    for i in range(len(ms)):
        ms[i].scale(1 / 1000)

    # depending of the numerical domain different height is used
    if typeNumDom == '3D':
        heightDist = hd
    elif typeNumDom == 'Q3D':
        heightDist = ms
    else:
        print("%s is an invalid numerical domain type... (3D, Q3D) Exiting... " % (typeNumDom))
        sys.exit()

    return heightDist

def getHeight(pt, HD):
    ###############################
    # calculating the height of point (z component) depending on the radius
    # input:   pt (point which height needs to be calculated)
    #          HD (height distribution)
    # output:  height

    if pt.rad > HD[len(HD)-1].rad:
        # for the stator (constant height)
        height = HD[len(HD)-1].Z

    else:
        # for the rotor (variable height)
        j = 0
        size = len(HD)
        while (pt.rad > HD[j].rad) and (j < size):
            j += 1

        dZ = HD[j].Z   - HD[j-1].Z
        dR = HD[j].rad - HD[j-1].rad

        height = HD[j-1].Z + (pt.rad - HD[j-1].rad)*(dZ/dR)

        if j==0:
            # extrapolating near the outflow of the rotor
            dZ = HD[j+1].Z   - HD[j].Z
            dR = HD[j+1].rad - HD[j].rad
            height = HD[j].Z + (pt.rad - HD[j].rad) * (dZ / dR)

    return height


def createCurve(pt, n, ang=-20, HD = None):
    ###############################
    # circular curve consisting of n points
    # input:   pt (starting point: x,y,z)
    #          n  (numb of points in the curve)
    #          ang(angle in degrees till the curve goes)
    #          HD (height distribution
    # output:  x  (x array of the points)
    #          y  (y array of the points)
    #          z  (z array of the points)
    # description: the curve start in pt and then
    # it is rotated by "angle" degrees.

    x   = np.zeros(n)
    y   = np.zeros(n)
    z   = np.zeros(n)
    for i in range(n):
        tmp = Point(x=pt.X, y=pt.Y, z=pt.Z)
        tmp.rotZ(i*ang/n)
        x[i] = tmp.X
        y[i] = tmp.Y
        z[i] = tmp.Z
        if HD:
            z[i] = getHeight(tmp, HD)

    return x, y, z


def createLine(pt1, pt2, n,  HD = None):
    ###############################
    # line consisting of n points
    # input:   pt1  (Starting point of the line: x,y,z)
    #          pt2  (Ending point of the line: x,y,z)
    #          n    (numb of points in the curve)
    # output:  x    (x array of the points)
    #          y    (y array of the points)
    #          z  (z array of the points)
    # observations: this function can only be used for the stator due to the height distribution!!

    pt  = pt2 - pt1
    vec = np.array([float(pt.X), float(pt.Y)])/n
    x   = np.zeros(n)
    y   = np.zeros(n)
    z   = np.zeros(n)
    for i in range(n):
        tmp = Point(x=float(pt1.X) + (i*vec[0]), y=float(pt1.Y) + (i*vec[1]), z=pt1.Z)
        x[i] = tmp.X
        y[i] = tmp.Y
        z[i] = tmp.Z
        if HD:
            z[i] = getHeight(tmp, HD)

    return x, y, z

def getArcLenght(IN, radius, angle):
    ###############################
    #  calculating the arc length for the area calculation
    # input:   IN       (input file)
    #          radius   (radius at which to calculated the arc)
    #          angle    (angle of the arc)
    # output:  area     (arc lenght, not actually an area)
    # observations: this function was made specially for the rotor channel were the arc is not complete due to the blades!

    radRotorLE =  float(IN['radLERotor'])
    radRotorTE =  float(IN['radTERotor'])

    if (radius > radRotorTE) and (radius < radRotorLE):
        # This only applies for the rotor channel because of the blade the arc length is not completely dtheta.
        radArc = readPointsFile(IN['fileRadArc'], 'RTHETA')

        #  interpolating in the distribution to find the actual arc angle
        j = 0
        size = len(radArc)
        while (radius > radArc[j].rad) and (j < size):
            j += 1

        dt = radArc[j].theta - radArc[j - 1].theta
        dR = radArc[j].rad   - radArc[j - 1].rad

        angle = radArc[j - 1].theta + (radius - radArc[j - 1].rad) * (dt / dR)

        angle = 2*angle         # the arc in the distribution is give for one blade passage, but actually we have two.
        area = ((4 * np.arctan(1.0) / 180. * angle) * radius)
    else:

        area=  ((4 * np.arctan(1.0) / 180. * angle) * radius)

    return area

def getLineDistance(IN, pt1, pt2, n):
    ###############################
    #  calculating the distance for the area calculation
    # input:   IN       (input file)
    #          pt1  (Starting point of the line: x,y,z)
    #          pt2  (Ending point of the line: x,y,z)
    # output:  area     (arc lenght, not actually an area)
    # observations: this function was made specially for the rotor channel were the arc is not complete due to the blades!

    radRotorLE =  float(IN['radLERotor'])
    radRotorTE =  float(IN['radTERotor'])

    radius= (pt1.rad+pt2.rad)*0.5

    if (radius > radRotorTE) and (radius < radRotorLE):
        # This only applies for the rotor channel because the height changes in a line
        hd = readPointsFile(IN['fileHD'], 'RZ')

        pt   = pt2 - pt1
        vec  = np.array([float(pt.X), float(pt.Y), float(pt.Z)]) / n
        magn = pt.magnitude()
        area = np.zeros(n)
        for i in range(n):
            tmp     = Point(x=float(pt1.X) + (i * vec[0]), y=float(pt1.Y) + (i * vec[1]), z=float(pt1.Y) + (i * vec[2]))
            area[i] = getHeight(tmp, hd) * magn/n
    else:

        area= pt1.distance2D(pt2)/ n

    return area


def getAverage(IN, x, y, z, propName='aver_entropy', typeline='ARC', typeNumDom='Q3D', typeavg='MASS', nameSpan=None, dtheta=20):
    ###############################
    # calculate the mass/area average of a property over a line/arc
    # input:   IN           (input file)
    #          x, y, z      (coordinates of the curve (can be a line or arc))
    #          propName     (name of the property that is being averaged)
    #          typeline     (type of curve, can be either line or arc)
    #          typeNumDom   (type of numerical domian: 3D or Q3D)
    #          typeAvg      (type of averaging method: MASS, AREA or other)
    #          nameSpan     (name of the span: hub, mid, shr)
    #          dtheta       (only for ARC curves, is the total arc lenght in degrees)
    # output:  average

    if typeNumDom == '3D':
        mprop, mass = getAverage3D(IN, x, y, z, propName, typeline, typeavg, nameSpan, dtheta)
    elif typeNumDom == 'Q3D':
        mprop, mass = getAverage2D(IN, x, y, z, propName, typeline, typeavg, dtheta)
    else:
        print("%s is an invalid numerical domain type... (3D, Q3D) Exiting... " % (typeNumDom))
        sys.exit()

    return mprop / mass


def getAverage2D(IN, x, y, z, propName, typeline, typeavg, dtheta, dz=1.0):
    ###############################
    # calculate the mass/area average of a property over a line/arc
    # input:   IN       (input file)
    #          x, y, z  (coordinates of the curve (can be a line or arc))
    #          propName (name of the property that is being averaged)
    #          typeline (type of curve, can be either line or arc)
    #          typeavg  (type of average, can be either area or mass based)
    #          dtheta   (only for ARC curves, is the total arc lenght in degrees)
    #          dz       (only for 3D, is discretization in z-direction)
    # output:  average
    # observation: this function is used in getAverage3D() which is called at each z-coordinate. That is why dz
    # comes into the function. For quasi-3D we assumed dz=1.0

    # extracting the line from the tecplot file
    LINE    = tp.data.extract.extract_line(zip(x, y, z))
    n       = len(y) - 1

    # extracting variables
    rho     = (LINE.values('RHO').as_numpy_array())
    prop    = (LINE.values(propName).as_numpy_array())
    vel_x   = (LINE.values('vrel-X').as_numpy_array())
    vel_y   = (LINE.values('vrel-Y').as_numpy_array())
    vel_z   = (LINE.values('vrel-Z').as_numpy_array())
    angle   = np.arctan((LINE.values('Y').as_numpy_array())/(LINE.values('X').as_numpy_array()))/ (4 * np.arctan(1.0) / 180.)


    # to check if all the line is include in the domain. In the case it is not, the lenght is changed
    #  which corresponds to the amount of intersection with the tecplot file
    if (n+1)!=len(rho):
        n = len(rho) - 1

    momentum = np.zeros(n + 1)
    mom      = np.sqrt(vel_x ** 2 + vel_y ** 2 + vel_z ** 2) * rho
    m        = np.ones(n + 1)

    # calculating the total and partial area
    if typeline == 'ARC':
        pt1     = Point(x=x[0], y=y[0], z=z[0])
        radius  = pt1.rad
        dArea   = getArcLenght(IN, radius, dtheta) * dz / n

        for i in range(n+1):
            # calculating the normal unit vector to the surface
            normVect = Point(rad=radius, theta=angle[i], z=0)       # unit vector normal to the surface
            # normVect.rotZ(180)                                      # the flow goes inwards
            normVect.scale(1 / normVect.magnitude())                # making it unit vector
            # projecting the velocity to the normal vector
            momentum[i] = np.abs((vel_x[i] * normVect.X * rho[i]) + (vel_y[i] * normVect.Y * rho[i])
                                 + (vel_z[i] * normVect.Z * rho[i]))
    else:
        pt1     = Point(x=x[0], y=y[0], z=z[0])
        pt2     = Point(x=x[-1], y=y[-1], z=z[-1])
        dArea   = getLineDistance(IN,pt1, pt2,n+1) *dz

        # calculating the normal unit vector to the surface
        normVect = pt2 - pt1
        normVect.rotZ(-90)                                  # making it normal to the line
        normVect.scale(1/normVect.magnitude())              # making it unit vector
        # projecting the velocity to the normal vector
        for i in range(n + 1):
            momentum[i] = np.abs((vel_x[i] * normVect.X * rho[i]) + (vel_y[i] * normVect.Y * rho[i])
                                 + (vel_z[i] * normVect.Z * rho[i]))

     # calculating the average
    if typeavg == 'MASS':
        m_prop = np.sum(momentum * prop * dArea)
        m_mass = np.sum(momentum * dArea)
    elif typeavg == 'AREA':
        m_prop = np.sum(m * prop * dArea)
        m_mass = np.sum(m * dArea)
    else:
        m_prop = np.sum(mom*prop * dArea)
        m_mass = np.sum(mom*dArea)


    return m_prop, m_mass

def getAverage3D(IN, x, y, z, propName, typeline, typeavg, nameSpan, dtheta):
    ###############################
    # calculate the mass/area average of a property over a line/arc
    # input:   IN       (input file)
    #          x, y, z  (coordinates of the curve (can be a line or arc))
    #          propName (name of the property that is being averaged)
    #          typeline (type of curve, can be either line or arc)
    #          typeNumDom   (type of numerical domian: 3D or Q3D)
    #          typeAvg      (type of averaging method: MASS, AREA or other)
    #          nameSpan (name of the span: hub, mid, shr)
    #          dtheta   (only for ARC curves, is the total arc length in degrees)
    #---------------------------------------------------------------------------
    #          lowSpan  (lower bound of the spanwise range averaging)
    #          topSpan  (upper bound of the spanwise range averaging)
    #          nSpan    (discretization in the 3rd direction)
    # output:  average
    # observations: this functions assumes that the numerical domain has as lower z-coordinate bound: zmin=0
    #   because the top and lower boundaries can't be intersected a tolerance has been added.

    # reading from IN file
    nSpan   = int(IN['nSPAN'])
    # nameSpan=IN['nameSPAN']


    tol     = z[0]*2/100                        # included because I do not want to be in the border of the domain
    if nameSpan == 'hub':
        zmin = 0.0 + tol
        zmax = z[0] * 0.35
    elif nameSpan=='mid':
        zmin = z[0] * 0.35
        zmax = z[0] * 0.65
    elif nameSpan=='shr':
        zmin = z[0] * 0.65
        zmax = z[0] - tol
    else:
        zmin = 0.0  + tol                             # assuming the numerical domain starts at z=0
        zmax = z[0] - tol



    # zmin = zmin + tol
    # zmax = zmax - tol
    dz = (zmax - zmin) / (nSpan - 1)

    m_prop  = np.zeros(nSpan)
    m_mass  = np.zeros(nSpan)


    for k in range(nSpan):

        # updating the z-coordinate of the line
        for j in range(len(y)):
            z[j] = zmin + k*dz

        # print('Height z= %.6f with an zmax=%.6f and zmin=%.6f at the span discretization %d' % (z[0], zmax, zmin, k))
        m_prop[k], m_mass[k] = getAverage2D(IN, x, y, z, propName, typeline, typeavg, dtheta, dz)


    m_total_prop = np.sum(m_prop)
    m_total_mass = np.sum(m_mass)

    return m_total_prop,  m_total_mass


def getAverageInLine(IN, nameLine, heightDist, propName, n=100, typeNumDom='Q3D', typeavg='MASS', nameSpan=None):
    ###############################
    # calculate the average value of a property over a line
    # input:   IN           (input file)
    #          nameLine     (name related to the line)
    #          n            (discretization of the curve for the averaging)
    #          heightDist   (height distribution)
    #          propName     (name of the property that is being averaged)
    #          typeNumDom   (type of numerical domian: 3D or Q3D)
    #          typeAvg      (type of averaging method: MASS, AREA or other)
    #          nameSpan     (name of the span: hub, mid, shr)
    # output:  prop     (averaged value of the property)
    #          radi     (radius, is an approximation, where the avg is calculated)

    #   reading the points that generate the line
    Pt1     = readPointCart(IN, nameLine + "_A")
    Pt2     = readPointCart(IN, nameLine + "_B")

    # generating the line
    x, y, z = createLine(Pt1, Pt2, n, heightDist)

    # calculating the average and the radius
    prop = float(getAverage(IN, x, y, z, propName, 'LINE', typeNumDom, typeavg, nameSpan))
    radi = (Pt1.rad + Pt2.rad)*0.5

    # printing results
    # print(nameLine, ': ', propName, ': ', prop, ' radius: ', radi)

    return prop, radi


def getAverageInLine_loop(IN, name, nlines, heightDist, propName, n=100, typeNumDom='Q3D', typeavg='MASS', nameSpan=None):
    ###############################
    # calculate the average value of a property over several lines
    # input:   IN           (input file)
    #          name         (name related to the line)
    #          nlines       (number of lines)
    #          n            (discretization of the curve for the averaging)
    #          heightDist   (height distribution)
    #          propName     (name of the property that is being averaged)
    #          typeNumDom   (type of numerical domian: 3D or Q3D)
    #          typeAvg      (type of averaging method: MASS, AREA or other)
    #          nameSpan     (name of the span: hub, mid, shr)
    # output:  prop     (averaged value of the property, is an array)
    #          radi     (radius, is an approximation, where the avg is calculated, is an array)
    #          nameLine (name related to the arc, is an array)
    # observations: the outputs are arrays,

    tmplist = []
    tmpProp = []
    tmpRad = []
    nameLine = np.array(tmplist, dtype=np.str_)
    property = np.array(tmpProp, dtype=np.float32)
    radius = np.array(tmpRad, dtype=np.float32)

    for i in range(nlines):

        nameL = name+str(i+1)
        prop, radi = getAverageInLine(IN, nameL, heightDist, propName, n, typeNumDom, typeavg, nameSpan)

        property = np.append(property, prop)
        radius = np.append(radius, radi)
        nameLine = np.append(nameLine, nameL)

    return property, radius, nameLine


def getAverageInArc_loop(IN, name, radIn, radOut, n, heightDist, propName, theta,
                         typeNumDom='Q3D', typeavg='MASS', nameSpan=None, nLine=200, rotAngle=-40.0):
    ###############################
    # calculating the average of a property between two radius with a certain discretization
    # input:   IN           (input file)
    #          name         (name related to the curve)
    #          radIn        (starting radius)
    #          radOut       (ending radius)  note:  radIn > radOut
    #          n            (discretization in radius, between radIn & radOut)
    #          heightDist   (height distribution)
    #          propName     (name of the property that is being averaged)
    #          theta        (angle at which to define the points with the radius)
    #          typeNumDom   (type of numerical domian: 3D or Q3D)
    #          typeAvg      (type of averaging method: MASS, AREA or other)
    #          nameSpan     (name of the span: hub, mid, shr)
    #          nLine        (discretization of the curve)
    #          rotAngle     (length of the generate arc)
    # output:  prop     (averaged value of the property, is an array)
    #          radi     (radius, is an approximation, where the avg is calculated, is an array)
    #          nameArc  (name related to the arc, is an array)
    # observations: the outputs are arrays, and the inputs: radIn > radOut

    dr = (radIn - radOut) / n

    tmplist = []
    tmpProp = []
    tmpRad = []
    nameLine     = np.array(tmplist, dtype=np.str_)
    property    = np.array(tmpProp, dtype=np.float32)
    radius      = np.array(tmpRad, dtype=np.float32)

    for i in range(n):
        
        # print('Radius: %.3f with an %d out of %d' % ((radIn - i * dr), i, n))
        Pt = Point(rad=(radIn - i * dr), theta=theta)
        x, y, z = createCurve(Pt, nLine, rotAngle, heightDist)
        prop = float(getAverage(IN, x, y, z, propName, 'ARC', typeNumDom, typeavg, nameSpan))

        property    = np.append(property, prop)
        radius      = np.append(radius, float(Pt.rad))
        nameLine    = np.append(nameLine, name+'_'+str(i))

        # print(name, "_", i, ' ', propName, ': ', prop, ' radius: ', float(Pt.rad))

    return property, radius, nameLine


def getAverageStator(IN, n, heightDist, propName, typeNumDom='Q3D', typeavg='MASS', nameSpan=None):
    ###############################
    # calculating the average of a property in the whole stator
    # input:   IN           (input file)
    #          n            (discretization in radius, for when the function getAverageInArc_loop() is used)
    #          heightDist   (height distribution)
    #          propName     (name of the property that is being averaged)
    #          typeNumDom   (type of numerical domian: 3D or Q3D)
    #          typeAvg      (type of averaging method: MASS, AREA or other)
    #          nameSpan     (name of the span: hub, mid, shr)
    # output:  property     (averaged value of the property, is an array)
    #          radius       (radius, is an approximation, where the avg is calculated, is an array)
    #          curveNames   (name related to the arc, is an array)
    # observations: the outputs are arrays

    tmplist = []
    tmpProp = []
    tmpRad = []
    curveNames  = np.array(tmplist, dtype=np.str_)
    property    = np.array(tmpProp, dtype=np.float32)
    radius      = np.array(tmpRad, dtype=np.float32)

    # reading points needed
    pt_inflow   = readPointCart(IN, 'INFLOW')
    pt_sle      = readPointCart(IN, 'STATOR_LE')
    pt_ste      = readPointCart(IN, 'STATOR_TE')
    pt_intL     = readPointCart(IN, 'INTERFACE_L')

    ################################################
    # from inflow to Leading edge of the stator (curves are defined as ARC)
    prop, radi, nameLine = getAverageInArc_loop(IN, 'StatorLE', pt_inflow.rad, pt_sle.rad, n,
                                                heightDist, propName, pt_inflow.theta, typeNumDom, typeavg, nameSpan, 150, -30)
    property = np.append(property, prop)
    radius = np.append(radius, radi)
    curveNames = np.append(curveNames, nameLine)


    ################################################
    #  in the stator channel  (curves are defined as LINE)

    prop, radi, nameLine = getAverageInLine_loop(IN, "TSLE", 3, heightDist, propName, 100, typeNumDom, typeavg, nameSpan)
    property = np.append(property, prop)
    radius = np.append(radius, radi)
    curveNames = np.append(curveNames, nameLine)

    prop, radi = getAverageInLine(IN, "THROAT", heightDist, propName, 100, typeNumDom, typeavg, nameSpan)
    property = np.append(property, prop)
    radius = np.append(radius, radi)
    curveNames = np.append(curveNames, "THROAT")

    prop, radi, nameLine = getAverageInLine_loop(IN, "TSTE", 7, heightDist, propName, 100, typeNumDom, typeavg, nameSpan)
    property = np.append(property, prop)
    radius = np.append(radius, radi)
    curveNames = np.append(curveNames, nameLine)

    ################################################
    # from interface to Leading edge of the rotor (curves are defined as ARC)
    prop, radi, nameLine = getAverageInArc_loop(IN, 'StatorTE', pt_ste.rad, pt_intL.rad, n,
                                                heightDist, propName, pt_ste.theta, typeNumDom, typeavg, nameSpan, 100, -20)
    property = np.append(property, prop)
    radius = np.append(radius, radi)
    curveNames = np.append(curveNames, nameLine)

    return property, radius, curveNames


def getAverageRotor(IN, n, heightDist, propName, typeNumDom='Q3D', typeavg='MASS', nameSpan=None):
    ###############################
    # calculating the average of a property in the whole rotor
    # input:   IN           (input file)
    #          n            (discretization in radius, for when the function getAverageInArc_loop() is used)
    #          heightDist   (height distribution)
    #          propName     (name of the property that is being averaged)
    #          typeNumDom   (type of numerical domian: 3D or Q3D)
    #          typeAvg      (type of averaging method: MASS, AREA or other)
    #          nameSpan     (name of the span: hub, mid, shr)
    # output:  property     (averaged value of the property, is an array)
    #          radius       (radius, is an approximation, where the avg is calculated, is an array)
    #          curveNames   (name related to the arc, is an array)
    # observations: the outputs are arrays

    tmplist = []
    tmpProp = []
    tmpRad = []
    curveNames = np.array(tmplist, dtype=np.str_)
    property = np.array(tmpProp, dtype=np.float32)
    radius = np.array(tmpRad, dtype=np.float32)

    pt1 = readPointCart(IN, 'INTERFACE_L')

    # from interface to Leading edge of the rotor
    prop, radi, nameLine = getAverageInArc_loop(IN, 'RotorLE', float(IN['radInterface']), float(IN['radLERotor']), n,
                                                heightDist, propName, pt1.theta, typeNumDom, typeavg, nameSpan)
    property = np.append(property, prop)
    radius = np.append(radius, radi)
    curveNames = np.append(curveNames, nameLine)


    # from Leading edge to trailing edge of the rotor
    prop, radi, nameLine = getAverageInArc_loop(IN, 'RotorChannel', float(IN['radLERotor']),
                                                float(IN['radTERotor']), n * 4, heightDist, propName, pt1.theta,
                                                typeNumDom, typeavg, nameSpan)
    property = np.append(property, prop)
    radius = np.append(radius, radi)
    curveNames = np.append(curveNames, nameLine)

    # from trailing edge of the rotor to the outflow
    prop, radi, nameLine = getAverageInArc_loop(IN, 'RotorTE', float(IN['radTERotor']), float(IN['radOutflow']), n,
                                                heightDist, propName, pt1.theta, typeNumDom, typeavg, nameSpan)
    property = np.append(property, prop)
    radius = np.append(radius, radi)
    curveNames = np.append(curveNames, nameLine)

    return property, radius, curveNames


def writedat(filename, property, radius, list_curve):
    with open(filename,'w') as f:
        for a, b, c in zip(property, radius, list_curve):
            print("%8.5e\t%8.5e\t%8.5e\t%s" % (a, b, c), file=f)

    print(filename, " has been written.")

########################################################################################################################
