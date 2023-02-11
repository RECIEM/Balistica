# Red Ciudadana de Estaciones Meteorologicas
#
# Copyright @ 2021
#
# Authors: Santiago Nunez-Corrales <snunezcr@gmail.com>
#          Jose Brenes-Andre <jbrenes54@gmail.com>
import numpy as np
import pandas as pd
from balistica.PhysicsEngine.PhysicsHandler import PhysicsHandler


class IdealPhysicsHandler(PhysicsHandler):
    def __init__(self, v0=0, theta=0, h=0, distance=-1):
        self.v0 = v0
        self.theta = theta
        self.b = 0
        self.height= h
        self.distance = distance
        self.data = None
        self.barrier = False

    def __x(self, t):
        return self.v0 * np.cos(self.theta) * t

    def __y(self, t):
        return (self.v0 * np.sin(self.theta) * t) - 0.5 * (self.g * np.power(t, 2))

    def __vx(self, t):
        return self.v0 * np.cos(self.theta) + t*0       # Horrible vectorial hack

    def __vy(self, t):
        return (self.v0 * np.sin(self.theta)) - (self.g * t)

    def compute(self):
        tstart = 0
        tend = 200
        tsamples = 10001
        trng = np.linspace(tstart, tend, tsamples)
        xrng = self.__x(trng)
        yrng = self.__y(trng)
        vxrng = self.__vx(trng)
        vyrng = self.__vy(trng)
        vrng = np.sqrt(np.power(vxrng, 2) + np.power(vyrng, 2))
        darray = np.transpose(np.array([trng, xrng, yrng, vxrng, vyrng, vrng]))
        self.data = pd.DataFrame({'t':darray[:,0], 'x':darray[:,1], 'z':darray[:,2], 'vx':darray[:,3], 'vy':darray[:,4], 'v':darray[:,5]})

        if self.barrier:
            self.data = self.data[self.data['x'] <= self.distance]

    def save_csv(self, filename):
        if (filename == '') or (self.data is None):
            return
        else:
            self.data.to_csv(filename)

    def maxT(self):
        return (self.v0 * np.sin(self.theta))/self.g

    def maxH(self):
        return (np.power(self.v0*np.sin(self.theta), 2))/(2*self.g)

    def totalR(self):
        if self.data is None:
            return 0.0
        else:
            adjdata = self.data[self.data['z'] >= np.min([0, self.height])]
            return adjdata.tail(1)['x'].values[0]

    def maxDistance(self):
        if self.data is None:
            return 0.0
        else:
            return self.data['x'].max()

    def totalT(self):
        if self.data is None:
            return 0.0
        else:
            adjdata = self.data[self.data['z'] >= np.min([0, self.height])]
            return adjdata.tail(1)['t'].values[0]

    def finalTheta(self):
        if self.data is None:
            return 0.0
        else:
            adjdata = self.data[self.data['z'] >= np.min([0, self.height])]
            return -1*np.rad2deg(np.arctan(adjdata.tail(1)['vy'].values[0]/adjdata.tail(1)['vx'].values[0]))

    def finalV(self):
        if self.data is None:
            return 0.0
        else:
            adjdata = self.data[self.data['z'] >= np.min([0, self.height])]
            return adjdata.tail(1)['v'].values[0]
