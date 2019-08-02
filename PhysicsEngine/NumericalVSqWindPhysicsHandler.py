# Red Ciudadana de Estaciones Meteorologicas
#
# Copyright @ 2019
#
# Author: Santiago Nunez-Corrales <snunezcr@gmail.com>

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.integrate import cumtrapz
from PhysicsEngine import PhysicsHandler


class NumericalVSqWindPhysicsHandler(PhysicsHandler):

    def __init__(self, v0=0, theta=0, b=1, height=-1, distance=-1):
        self.v0 = v0
        self.theta = theta
        self.b = b
        self.height = height
        self.distance = distance
        self.data = None
        self.windx = 0

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
            dvxdt = -self.b * vx * np.sqrt(np.power(vx, 2) + np.power(vy, 2)) - self.windx
            dvydt = -self.g - self.b * vy * np.sqrt(np.power(vx, 2) + np.power(vy, 2))
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
            {'t': darray[:, 0], 'x': darray[:, 1], 'y': darray[:, 2], 'vx': darray[:, 3], 'vy': darray[:, 4],
             'v': darray[:, 5]})
        self.data = self.data[self.data['y'] >= 0.0]

        if self.height >= 0:
            self.data = self.data[self.data['x'] <= self.distance]

        if self.height >= 0:
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
            return self.data[self.data['y'] == self.data['y'].max()]['t'].values[0]

    def maxH(self):
        if self.data is None:
            return 0.0
        else:
            return self.data[self.data['y'] == self.data['y'].max()]['y'].values[0]

    def totalR(self):
        if self.data is None:
            return 0.0
        else:
            adjdata = self.data[self.data['y'] >= self.height]
            return adjdata.tail(1)['x'].values[0]

    def totalT(self):
        if self.data is None:
            return 0.0
        else:
            adjdata = self.data[self.data['y'] >= self.height]
            return adjdata.tail(1)['t'].values[0]

    def finalTheta(self):
        if self.data is None:
            return 0.0
        else:
            adjdata = self.data[self.data['y'] >= self.height]
            return -1 * np.rad2deg(np.arctan(adjdata.tail(1)['vy'].values[0] / adjdata.tail(1)['vx'].values[0]))

    def finalV(self):
        if self.data is None:
            return 0.0
        else:
            adjdata = self.data[self.data['y'] >= self.height]
            return adjdata.tail(1)['v'].values[0]

