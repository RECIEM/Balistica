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


class NumericalThermoFluidHandler(PhysicsHandler):
    # Constants for all instances of the solver
    ###########################################

    # Sphericity
    sa_norm = 1.6075

    # Sutherland's viscosity model
    tref = 273.15
    tsuth = 110.4
    muref = 1.716e-5

    def __init__(self, v0=0, theta=0, dens=0.7, a=0.05, b=0.05, c=0.05, rho = 1.2754, temp = 293.0, h = 0, d = 0):
        self.v0 = v0
        self.theta = theta
        self.a = a
        self.b = b
        self.c = c
        self.rho = rho
        self.T = temp
        self.height = h
        self.distance = d
        self.windx = 0
        self.windz = 0
        self.data = None
        self.computeIdeal = False
        self.barrier = False
        # Intermediate sphericity-related values, compute as early as possible to avoid constant recomputation
        self.phi = self.sphericity
        self.A = np.exp(2.3288 - (6.4581 * self.phi) + 2.4486 * (self.phi ** 2))
        self.B = 0.0964 + (0.5565 * self.phi)
        self.C = np.exp(4.905 - (13.8944 * self.phi) + (18.4222 * (self.phi ** 2)) - (10.2599 * (self.phi ** 3)))
        self.D = np.exp(1.4681 + (12.2584 * self.phi) - (20.7322 * (self.phi ** 2)) + (15.8855 * (self.phi ** 3)))
        self.dSph = np.power(self.a * self.b * self.c, 1.0/3)
        # Compute kinematic viscosity only once
        self.mu = self.compMu()
        self.dens = dens
        self.m = 0

    @property
    def sphericity(self):
        radius = np.power(self.a * self.b * self.c, 1.0/3)
        sphArea = 4 * np.pi * (radius ** 2)
        parArea = self.surfArea
        return sphArea/parArea

    @property
    def sph_volume(self):
        radius = np.power(self.a * self.b * self.c, 1.0 / 3)
        return (4.0/3.0) * np.pi * (radius ** 3)

    def compMu(self):
        mu = self.muref * np.power(self.T / self.tref, 3.0/2) * ((self.tref + self.tsuth) / (self.T + self.tsuth))
        return mu

    def Re(self, v):
        return (self.dSph * v * self.rho) / self.mu

    def Cd(self, v):
        return (24.0 / self.Re(v))*(1 + (self.A * (self.Re(v)**self.B))) + (self.C / (1 + (self.D / self.Re(v))))

    @staticmethod
    def norm(a, b):
        return np.sqrt(np.power(a, 2) + np.power(b, 2))

    @property
    def surfArea(self):
        projareap = np.power(self.a*self.b, self.sa_norm) + np.power(self.a*self.c, self.sa_norm) + np.power(self.b*self.c, self.sa_norm)
        return 4*np.pi*np.power(projareap/3.0, 1.0/self.sa_norm)

    def E(self, v):
        return (self.rho*v*v*self.surfArea)/(2 * self.m)

    def setMass(self):
        self.m = self.dens * self.sph_volume

    def compute(self):
        # Update mu value
        self.phi = self.sphericity
        self.A = np.exp(2.3288 - (6.4581 * self.phi) + 2.4486 * (self.phi ** 2))
        self.B = 0.0964 + (0.5565 * self.phi)
        self.C = np.exp(4.905 - (13.8944 * self.phi) + (18.4222 * (self.phi ** 2)) - (10.2599 * (self.phi ** 3)))
        self.D = np.exp(1.4681 + (12.2584 * self.phi) - (20.7322 * (self.phi ** 2)) + (15.8855 * (self.phi ** 3)))
        self.dSph = np.power(self.a * self.b * self.c, 1.0/3)
        self.mu = self.compMu()

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

            dvxdt = -Enow * self.Cd(v) * (vx - self.windx)/v
            dvydt = -Enow * self.Cd(v) * (vy / v) - self.g
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

        # Record Reynolds number
        cdrng = self.Cd(vrng)
        rerng = self.Re(vrng)

        darray = np.transpose(np.array([trng, xrng, yrng, vxrng, vyrng, vrng, cdrng, rerng]))

        self.data = pd.DataFrame(
            {'t': darray[:, 0], 'x': darray[:, 1], 'z': darray[:, 2], 'vx': darray[:, 3], 'vz': darray[:, 4],
             'v': darray[:, 5], 'cd': darray[:, 6], 're': darray[:, 7]})

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
