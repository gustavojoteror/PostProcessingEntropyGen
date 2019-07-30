import numpy as np


class Point(object):

    def __init__(self, x=None, y=None, z=None, rad=None, theta=None):

        # constructor with cartesian coordinates
        if x:
            self.X = x
            self.Y = y
            self.calcPolarCoord()
            #print("Point initialized with cartesian coordinates")

        # constructor with polar coordinates
        if rad:
            self.rad    = rad
            self.theta  = theta
            self.calcCartesian()
            #print("Point initialized with polar coordinates")

        if z:
            self.Z = z
        else:
            self.Z = 0

    ######################################
    # Function to calculate cartesian and polar coordinates from each other
    def calcPolarCoord(self):
        self.rad   = np.sqrt(self.X**2 + self.Y**2)
        self.theta = np.arctan(self.Y / self.X) / (4 * np.arctan(1.0) / 180.)

    def calcCartesian(self):
        self.X = self.rad * np.cos(self.theta * (4 * np.arctan(1.0) / 180.))
        self.Y = self.rad * np.sin(self.theta * (4 * np.arctan(1.0) / 180.))

    ######################################
    # Function to print variables out of the object
    def printCart(self):
        print("x: ", self.X, " y: ", self.Y, " z: ", self.Z)

    def printPolar(self):
        print("radius: ", self.rad, " theta: ", self.theta)

    ######################################
    # Function to get variables out of the object
    def getCartesian2D(self):
        return self.X, self.Y

    def getCartesian3D(self):
        return self.X, self.Y, self.Z

    def getPolar(self):
        return self.rad, self.theta

    ######################################
    # Function to translate the point
    def moveX(self, dx):
        self.X = self.X + dx
        self.calcPolarCoord()

    def moveY(self, dy):
        self.Y = self.Y + dy
        self.calcPolarCoord()

    def moveZ(self, dz):
        self.Z = self.Z + dz
        self.calcPolarCoord()

    def moveXY(self, dx, dy):
        self.X = self.X + dx
        self.Y = self.Y + dy
        self.calcPolarCoord()

    def moveXYZ(self, dx, dy, dz):
        self.X = self.X + dx
        self.Y = self.Y + dy
        self.Z = self.Z + dz
        self.calcPolarCoord()

    ######################################
    # Magnitude of the point
    def magnitude(self):
        return np.sqrt(self.X ** 2 + self.Y ** 2 + self.Z ** 2)

    # Function to rotate the point
    def rotZ(self, angle):
        self.theta=self.theta+angle
        self.calcCartesian()

    ######################################
    # Function to change the point
    def scale(self, fact):
        self.X = fact * self.X
        self.Y = fact * self.Y
        self.Z = fact * self.Z
        self.calcPolarCoord()

    def __add__(self, other):
        return Point(x=self.X + other.X, y=self.Y + other.Y, z=self.Z + other.Z)

    def __sub__(self, other):
        return Point(x=self.X - other.X, y=self.Y - other.Y, z=self.Z - other.Z)

    def distance2D(self, other):
        dx = self.X - other.X
        dy = self.Y - other.Y
        return np.sqrt(dx ** 2 + dy ** 2)

    def dotProduct(self, other):
        return (self.X * other.X)+(self.Y * other.Y)+(self.Z * other.Z)

    def ProjectTo(self, other):
        return other.scale(((self.X * other.X)+(self.Y * other.Y)+(self.Z * other.Z))/other.magnitude())

