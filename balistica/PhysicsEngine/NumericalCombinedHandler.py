# Red Ciudadana de Estaciones Meteorologicas
#
# Copyright @ 2021
#
# Authors: Santiago Nunez-Corrales <snunezcr@gmail.com>
#          Jose Brenes-Andre <jbrenes54@gmail.com>
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.integrate import cumtrapz
from balistica.PhysicsEngine.PhysicsHandler import PhysicsHandler


class NumericalCombinedHandler(PhysicsHandler):
    sa_norm = 1.6075

    def __init__(self, v0=0, theta=0, dens=0.7, rho=1, a=0.05, b=0.05, c=0.05, Cd=1, height=0, distance=0):
        self.v0 = v0
        self.theta = theta
        self.rho = rho
        self.a = a
        self.b = b
        self.c = c
        self.Cd = Cd
        self.height = height
        self.distance = distance
        self.windx = 0
        self.data = None
        self.computeIdeal = False
        self.barrier = False
        self.dens = dens
        self.m = 0

    @property
    def sphericity(self):
        radius = np.power(self.a * self.b * self.c, 1.0 / 3)
        sphArea = 4 * np.pi * (radius ** 2)
        parArea = self.surfArea
        return sphArea / parArea

    @property
    def sph_volume(self):
        radius = np.power(self.a * self.b * self.c, 1.0 / 3)
        return (4.0 / 3.0) * np.pi * (radius ** 3)

    @property
    def surfArea(self):
        projareap = np.power(self.a * self.b, self.sa_norm) + np.power(self.a * self.c, self.sa_norm) + np.power(
            self.b * self.c, self.sa_norm)
        return 4 * np.pi * np.power(projareap / 3.0, 1.0 / self.sa_norm)

    def setMass(self):
        self.m = self.dens * self.sph_volume

    @staticmethod
    def norm(a, b):
        return np.sqrt(np.power(a, 2) + np.power(b, 2))


    def E(self, v):
        return (self.rho*v*v*self.surfArea)/(2*self.m)

    def compute(self):
        tstart = 0
        tend = 200
        tsamples = 10001
        trng = np.linspace(tstart, tend, tsamples)

        vx0 = self.v0 * np.cos(self.theta)
        vy0 = self.v0 * np.sin(self.theta)

        def acc(t, v):
            vx = v[0]
            vy = v[1]

            v = self.norm(vx, vy)
            Enow = self.E(v)

            dvxdt = -Enow * self.Cd * (vx - self.windx)/v
            dvydt = -Enow * self.Cd * (vy / v) - self.g
            return [dvxdt, dvydt]

        # Integrate velocities
        vel0 = [vx0, vy0]
        vel = solve_ivp(acc, [0, 200], vel0, method='RK45', t_eval=trng).y
        vxrng = vel[0]
        vyrng = vel[1]

        # Integrate positions
        xrng = cumtrapz(vxrng, trng, initial=0)
        yrng = cumtrapz(vyrng, trng, initial=0)

        vrng = np.sqrt(np.power(vxrng, 2) + np.power(vyrng, 2))
        darray = np.transpose(np.array([trng, xrng, yrng, vxrng, vyrng, vrng]))
        self.data = pd.DataFrame(
            {'t': darray[:, 0], 'x': darray[:, 1], 'z': darray[:, 2], 'vx': darray[:, 3], 'vz': darray[:, 4],
             'v': darray[:, 5]})

        if self.barrier:
            self.data = self.data[self.data['x'] <= self.distance]

    def save_csv(self, filename):
        if (filename == '') or (self.data is None):
            return
        else:
            self.data.to_csv(filename)

    def maxT(self):
        if self.data is None:
            return 0.0
        else:
            return self.data[self.data['z'] == self.data['z'].max()]['t'].values[0]

    def maxH(self):
        if self.data is None:
            return 0.0
        else:
            return self.data[self.data['z'] == self.data['z'].max()]['z'].values[0]

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
            adjdata = self.data[self.data['z'] >= np.min([0, self.height])]
            return adjdata['x'].max()

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

            if adjdata.tail(1)['vx'].values[0] == 0:
                return 90.0
            else:
                return -1 * np.rad2deg(np.arctan(adjdata.tail(1)['vz'].values[0] / adjdata.tail(1)['vx'].values[0]))

    def finalV(self):
        if self.data is None:
            return 0.0
        else:
            adjdata = self.data[self.data['z'] >= np.min([0, self.height])]
            return adjdata.tail(1)['v'].values[0]

