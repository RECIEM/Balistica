# Red Ciudadana de Estaciones Meteorologicas
#
# Copyright @ 2019
#
# Author: Santiago Nunez-Corrales <snunezcr@gmail.com>

import numpy as np
import pandas as pd
from scipy.integrate import odeint
from PhysicsEngine import PhysicsHandler


class NumericalVPhysicsHandler(PhysicsHandler):

    def __init__(self, v0=0, theta=0, b=1, height=-1, distance=-1):
        self.v0 = v0
        self.theta = theta
        self.b = b
        self.height = height
        self.distance = distance
        self.data = None

    def compute(self):
        tstart = 0
        tend = 200
        tsamples = 5001
        trng = np.linspace(tstart, tend, tsamples)

        def vx(x, t, b, v0, theta):
            return v0 * np.cos(theta) * np.exp(-b * t)

        def vy(y, t, g, b, v0, theta):
            return (((g / b) + (v0 * np.sin(theta))) * np.exp(-b * t)) - (g / b)

        vx0 = self.v0*np.cos(self.theta)
        vy0 = self.v0*np.sin(self.theta)

        vyrng = vx(vx0, trng, self.b, self.v0, self.theta)
        vxrng = vy(vy0, trng, self.g, self.b, self.v0, self.theta)

        # Integrate
        xrng = odeint(vx, 0.0, trng, args=(self.b, self.v0, self.theta)).flatten()
        yrng = odeint(vy, 0.0, trng, args=(self.g, self.b, self.v0, self.theta)).flatten()

        vrng = np.sqrt(np.power(vxrng, 2) + np.power(vyrng, 2))
        darray = np.transpose(np.array([trng, xrng, yrng, vxrng, vyrng, vrng]))
        self.data = pd.DataFrame(
            {'t': darray[:, 0], 'x': darray[:, 1], 'y': darray[:, 2], 'vx': darray[:, 3], 'vy': darray[:, 4],
             'v': darray[:, 5]})
        self.data = self.data[self.data['y'] >= 0.0]

        if self.height >= 0:
            self.data = self.data[self.data['x'] <= self.distance]

    def save_csv(self, filename):
        if (filename == '') or (self.data is None):
            return
        else:
            self.data.to_csv(filename)

    def maxT(self):
        return (1.0/self.b)*np.log(1.0 + ((self.b*self.v0*np.sin(self.theta))/self.g))

    def maxH(self):
        return ((self.v0 * np.sin(self.theta))/ self.b) - (self.g/np.power(self.b,2.0)) * np.log(1.0 + ((self.b * self.v0 * np.sin(self.theta)) / self.g))

    def totalR(self):
        if self.data is None:
            return 0.0
        else:
            return self.data.tail(1)['x'].values[0]

    def totalT(self):
        if self.data is None:
            return 0.0
        else:
            return self.data.tail(1)['t'].values[0]

    def finalTheta(self):
        if self.data is None:
            return 0.0
        else:
            return -1 * np.rad2deg(np.arctan(self.data.tail(1)['vy'].values[0] / self.data.tail(1)['vx'].values[0]))

    def finalV(self):
        if self.data is None:
            return 0.0
        else:
            return self.data.tail(1)['v'].values[0]
