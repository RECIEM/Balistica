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


class NumericalVSqWindPhysicsHandler(PhysicsHandler):
    # Sphericity
    sa_norm = 1.6075

    def __init__(self, v0=0, theta=0, d=1, height=0, distance=-1, a=0.01, b=0.01, c=0.01):
        self.v0 = v0
        self.theta = theta
        self.D = d
        self.height = height
        self.distance = distance
        self.data = None
        self.windx = 0
        self.barrier = False

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
            dvxdt = -self.D * vx * np.sqrt(np.power(vx, 2) + np.power(vy, 2)) - self.windx
            dvydt = -self.g - self.D * vy * np.sqrt(np.power(vx, 2) + np.power(vy, 2))
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
            return -1 * np.rad2deg(np.arctan(adjdata.tail(1)['vz'].values[0] / adjdata.tail(1)['vx'].values[0]))

    def finalV(self):
        if self.data is None:
            return 0.0
        else:
            adjdata = self.data[self.data['z'] >= np.min([0, self.height])]
            return adjdata.tail(1)['v'].values[0]
