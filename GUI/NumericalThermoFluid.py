# Red Ciudadana de Estaciones Meteorologicas
#
# Copyright @ 2019
#
# Author: Santiago Nunez-Corrales <snunezcr@gmail.com>

import sys
import numpy as np
import matplotlib
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
from PhysicsEngine import NumericalThermoFluidHandler
from PhysicsEngine import IdealPhysicsHandler
from tkinter import filedialog


class NumericalV2WindThermoGUI(tk.Frame):
    def __init__(self, master=None):
        self.physicshandler = NumericalThermoFluidHandler()
        self.idealphysicshandler = IdealPhysicsHandler()
        self.bounds = None
        self.goodparams = False

        tk.Frame.__init__(self, master)
        self.grid()

        # Top level panel structure
        self.panels = tk.Frame(self)
        self.panels.pack(fill=tk.BOTH, expand=1)

        # Left and right panels
        self.leftpanel = tk.Frame(self.panels, relief=tk.GROOVE)
        self.leftpanel.pack(side=tk.LEFT)
        self.rightpanel = tk.Frame(self.panels)
        self.rightpanel.pack(side=tk.RIGHT)

        # Controls grid for upper left pannel
        self.ulpanel = tk.LabelFrame(self.leftpanel, text='Parameters')
        self.ulpanel.pack(side=tk.TOP)

        # Control for density
        self.densitylable = tk.Label(self.ulpanel, text='Projectile density (kg/m^3)')
        self.densitylable.grid(row=0, column=0)
        self.densityinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.densityinput.grid(row=0, column=1)
        self.densityinput.insert(0, '3400.0')

        # Control for sphericity
        self.pspherlabel = tk.Label(self.ulpanel, text='Sphericity')
        self.pspherlabel.grid(row=1, column=0, columnspan=4)

        self.alabel = tk.Label(self.ulpanel, text='a (m)')
        self.alabel.grid(row=2, column=0)
        self.blabel = tk.Label(self.ulpanel, text='b (m)')
        self.blabel.grid(row=2, column=1)
        self.clabel = tk.Label(self.ulpanel, text='c (m)')
        self.clabel.grid(row=2, column=2)

        self.ainput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.ainput.grid(row=3, column=0)
        self.ainput.insert(0, '0.05')
        self.binput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.binput.grid(row=3, column=1)
        self.binput.insert(0, '0.05')
        self.cinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.cinput.grid(row=3, column=2)
        self.cinput.insert(0, '0.05')

        # Control for density
        self.rholable = tk.Label(self.ulpanel, text='Air density (kg/m^3)')
        self.rholable.grid(row=4, column=0)
        self.rhoinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.rhoinput.grid(row=4, column=1)
        self.rhoinput.insert(0, '1.2754')

        # Control for temperature
        self.templabel = tk.Label(self.ulpanel, text='Air temperature (K)')
        self.templabel.grid(row=5, column=0)
        self.tempinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.tempinput.insert(0, '293.0')
        self.tempinput.grid(row=5, column=1)

        # Control for angle
        self.anglelable = tk.Label(self.ulpanel, text='Initial angle (degrees)') # TODO: change all
        self.anglelable.grid(row=6, column=0)
        self.angleinput = tk.Scale(self.ulpanel, from_=0, to=90, resolution=1, length=170, orient=tk.HORIZONTAL)
        self.angleinput.grid(row=6, column=1)

        # Control for velocity
        self.velocitylabel = tk.Label(self.ulpanel, text='Initial velocity (m/s)')
        self.velocitylabel.grid(row=7, column=0)
        self.velocityinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.velocityinput.grid(row=7, column=1)
        self.velocityinput.insert(0, '125')

        self.latIlabel = tk.Label(self.ulpanel, text='I. Lat (m)')
        self.latIlabel.grid(row=8, column=0)
        self.lonIlabel = tk.Label(self.ulpanel, text='I. Lon (m)')
        self.lonIlabel.grid(row=8, column=1)
        self.heightIlabel = tk.Label(self.ulpanel, text='I. Height (m)')
        self.heightIlabel.grid(row=8, column=2)

        self.latIinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.latIinput.grid(row=9, column=0)
        self.latIinput.insert(0, '0')
        self.lonIinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.lonIinput.grid(row=9, column=1)
        self.lonIinput.insert(0, '0')
        self.heightIinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.heightIinput.grid(row=9, column=2)
        self.heightIinput.insert(0, '0')

        self.latFlabel = tk.Label(self.ulpanel, text='F. Lat (m)')
        self.latFlabel.grid(row=10, column=0)
        self.lonFlabel = tk.Label(self.ulpanel, text='F. Lon (m)')
        self.lonFlabel.grid(row=10, column=1)
        self.heightFlabel = tk.Label(self.ulpanel, text='F. Height (m)')
        self.heightFlabel.grid(row=10, column=2)

        self.latFinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.latFinput.grid(row=11, column=0)
        self.latFinput.insert(0, '100')
        self.lonFinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.lonFinput.grid(row=11, column=1)
        self.lonFinput.insert(0, '100')
        self.heightFinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.heightFinput.grid(row=11, column=2)
        self.heightFinput.insert(0, '0')


        self.barrierset = tk.BooleanVar()
        self.barriercheck = tk.Checkbutton(self.ulpanel, justify=tk.RIGHT, variable=self.barrierset, onvalue=True,
                                           offvalue=False, text='Show barrier', command=self.setBarrierHandler)
        self.barriercheck.grid(row=12, column=0)

        self.idealset = tk.BooleanVar()
        self.idealcheck = tk.Checkbutton(self.ulpanel, justify=tk.RIGHT, variable=self.idealset, onvalue=True,
                                           offvalue=False, text='Show ideal')
        self.idealcheck.grid(row=12, column=1)

        self.pwindlabel = tk.Label(self.ulpanel, text='Wind settings:')
        self.pwindlabel.grid(row=13, column=0, columnspan=2)

        self.windanlabel = tk.Label(self.ulpanel, text='Azimuth (degrees):')
        self.windanlabel.grid(row=14, column=0)
        self.windangle = tk.Scale(self.ulpanel, from_=0, to=359, resolution=1, length=200, orient=tk.HORIZONTAL)
        self.windangle.grid(row=14, column=1, columnspan=2)

        self.windmglabel = tk.Label(self.ulpanel, text='Magnitude (m/s):')
        self.windmglabel.grid(row=15, column=0)
        self.windmag = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.windmag.grid(row=15, column=1)

        self.windmag.insert(0, '0')

        # Controls grid for upper left pannel
        self.blpanel = tk.Frame(self.leftpanel)
        self.blpanel.pack(side=tk.BOTTOM)

        # Buttons for various functions
        self.blanklabel= tk.Label(self.blpanel, text="")
        self.blanklabel.grid(row=0, column=0, columnspan=2)

        self.computebutton = tk.Button(self.blpanel, text="Compute", width=20, command=self.compute, default=tk.NORMAL)
        self.computebutton.grid(row=1, column=0, columnspan=3)

        self.graphxt = tk.Button(self.blpanel, text="x(t) vs. t", width=10, command=self.txGraph, default=tk.NORMAL)
        self.graphxt.grid(row=2, column=0)

        self.graphzt = tk.Button(self.blpanel, text="z(t) vs. t", width=10, command=self.tyGraph, default=tk.NORMAL)
        self.graphzt.grid(row=2, column=1)

        self.graphvt = tk.Button(self.blpanel, text="v(t) vs. t", width=10, command=self.tvGraph, default=tk.NORMAL)
        self.graphvt.grid(row=2, column=2)

        self.graphcdt = tk.Button(self.blpanel, text="Cd(v) vs. t", width=10, command=self.cdtGraph, default=tk.NORMAL)
        self.graphcdt.grid(row=2, column=3)

        self.graphzx = tk.Button(self.blpanel, text="z(t) vs. x(t)", width=10, command=self.xyGraph, default=tk.NORMAL)
        self.graphzx.grid(row=3, column=0)

        self.graphvx = tk.Button(self.blpanel, text="v(t) vs. x(t)", width=10, command=self.xvGraph, default=tk.NORMAL)
        self.graphvx.grid(row=3, column=1)

        self.graphvz = tk.Button(self.blpanel, text="v(t) vs. z(t)", width=10, command=self.vyGraph, default=tk.NORMAL)
        self.graphvz.grid(row=3, column=2)

        self.graphcdv = tk.Button(self.blpanel, text="Cd(v) vs. v", width=10, command=self.cdvGraph, default=tk.NORMAL)
        self.graphcdv.grid(row=3, column=3)

        self.userlabel = tk.Label(self.blpanel, text="", fg="red")
        self.userlabel.grid(row=4, column=0, columnspan=3)

        self.csvbutton= tk.Button(self.blpanel, text="Save to CSV", command=self.saveCSV, default=tk.NORMAL)
        self.csvbutton.grid(row=5, column=0)

        self.pngbutton = tk.Button(self.blpanel, text="Save to PNG", command=self.savePNG, default=tk.NORMAL)
        self.pngbutton.grid(row=5, column=1)

        self.quitbutton = tk.Button(self.blpanel, text="Quit", command=self.bye, default=tk.NORMAL)
        self.quitbutton.grid(row=5, column=2)

        self.physicshandler.v0 = 0
        self.physicshandler.theta = 0
        self.physicshandler.b = 1

        fig, axs = plt.subplots(1, 1, figsize=(7, 6), dpi=80)
        axs.set_xlabel('Distance (m)')
        axs.set_ylabel('Height (m)')
        axs.set_xlim(0, 100)
        axs.set_ylim(0, 100)
        axs.set_title('Projectile ballistics')
        canvas = FigureCanvasTkAgg(fig, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()
        self.mostrecentfig = fig

    def setBarrierHandler(self):
        if self.barrierset.get():
            self.physicshandler.barrier = True
            self.idealphysicshandler.barrier = True
        else:
            self.physicshandler.barrier = False
            self.idealphysicshandler.barrier = False

    @staticmethod
    def norm(a, b):
        return np.sqrt(np.power(a, 2) + np.power(b, 2))

    def consumebounds(self):
        latI = 0.0
        try:
            latI = float(self.latIinput.get())
        except:
            self.userlabel['text'] = "Initial latitude format incorrect"
            self.bounds = None

        latF = 0.0
        try:
            latF = float(self.latFinput.get())
        except:
            self.userlabel['text'] = "Final latitude format incorrect"
            self.bounds = None

        lonI = 0.0
        try:
            lonI = float(self.lonIinput.get())
        except:
            self.userlabel['text'] = "Initial longitude format incorrect"
            self.bounds = None

        lonF = 0.0
        try:
            lonF = float(self.lonFinput.get())
        except:
            self.userlabel['text'] = "Final longitude format incorrect"
            self.bounds = None

        heightI = 0.0
        try:
            heightI = float(self.heightIinput.get())
        except:
            self.userlabel['text'] = "Initial latitude format incorrect"
            self.bounds = None

        heightF = 0.0
        try:
            heightF = float(self.heightFinput.get())
        except:
            self.userlabel['text'] = "Initial latitude format incorrect"
            self.bounds = None

        self.bounds = (latI, latF, lonI, lonF, heightI, heightF)

    def consumeparams(self):
        try:
            self.physicshandler.a = float(self.ainput.get())
        except:
            self.userlabel['text'] = "Sphericity for a, format incorrect"
            return

        try:
            self.physicshandler.b = float(self.binput.get())
        except:
            self.userlabel['text'] = "Sphericity for b, format incorrect"
            return

        try:
            self.physicshandler.c = float(self.cinput.get())
        except:
            self.userlabel['text'] = "Sphericity for c, format incorrect"
            return

        try:
            self.physicshandler.dens = float(self.densityinput.get())
        except:
            self.userlabel['text'] = "Density format incorrect"
            return

        # Compute mass explicitly
        self.physicshandler.setMass()

        try:
            self.physicshandler.rho = float(self.rhoinput.get())
        except:
            self.userlabel['text'] = "Air density, format incorrect"
            return

        try:
            # Convert to g/cm^3 to kg/m^3
            self.physicshandler.T = float(self.tempinput.get())
        except:
            self.userlabel['text'] = "Air temperature, format incorrect"
            return

        self.physicshandler.theta = np.deg2rad(float(self.angleinput.get()))

        if self.idealset.get():
            self.idealphysicshandler.theta = self.physicshandler.theta

        try:
            self.physicshandler.v0 = float(self.velocityinput.get())

            if self.idealset.get():
                self.idealphysicshandler.v0 = self.physicshandler.v0
        except:
            self.userlabel['text'] = "Velocity, format incorrect"
            return

        # Process lat, lon and height data
        self.consumebounds()

        # With that, compute the contribution of wind in x, z
        if self.bounds is not None:
            try:
                windmag = float(self.windmag.get())
            except:
                self.userlabel['text'] = "Wind speed magnitude, format incorrect"
                return

            try:
                windtheta = np.deg2rad(float(self.windangle.get()))
            except:
                self.userlabel['text'] = "Wind angle, format incorrect"
                return

            # Second, compute the linear transformation that computes the inner product of
            # wind contribution with respect to the direction of the trajectory
            latI, latF, lonI, lonF, _, _ = self.bounds

            dx = lonF - lonI
            dy = latF - latI

            if (dx <= 0) or (dy <= 0):
                self.userlabel['text'] = "Distance must be greater than zero"
                return

            beta = np.arctan(np.abs(dy / dx))

            if (dx > 0) and (dy > 0):
                azimuth = np.pi / 2 - windtheta
            elif (dx > 0) and (dy < 0):
                azimuth = np.pi / 2 + windtheta
            elif (dx < 0) and (dy < 0):
                azimuth = np.pi * (3.0 / 2) - windtheta
            else:
                azimuth = np.pi * (3.0 / 2) + windtheta

            self.physicshandler.windx = windmag * np.cos(azimuth - beta)
            self.goodparams = True

    def compute(self):
        self.userlabel['text'] = ""
        self.consumeparams()

        if not self.goodparams:
            return

        distance, height = self.geofigbounds()
        self.physicshandler.distance = distance
        self.physicshandler.height = height

        self.physicshandler.compute()

        if self.idealset.get():
            self.idealphysicshandler.compute()

        self.xyGraph()

    def geofigbounds(self):
        if self.bounds is not None:
            latI, latF, lonI, lonF, heightI, heightF = self.bounds
            distance = np.sqrt(np.power((latF - latI), 2) + np.power((lonF - lonI), 2))
            height = heightF - heightI
            return (distance, height)
        else:
            return (0, 0)

    def txGraph(self):
        for s in self.rightpanel.grid_slaves():
            s.destroy()

        figtx, axs = plt.subplots(1, 1, figsize=(7, 6), dpi=80)
        selected = self.physicshandler.data[self.physicshandler.data['t'] <= self.physicshandler.totalT()]
        axs.plot(selected['t'], selected['x'], '-', linewidth=2, color='b', label='With drag ~ v^2')

        if self.idealset.get():
            iselected = self.idealphysicshandler.data[
                self.idealphysicshandler.data['t'] <= self.idealphysicshandler.totalT()]
            axs.plot(iselected['t'], iselected['x'], '--', linewidth=1, color='g', label='Ideal')

        axs.set_xlabel('Time (s)')
        axs.set_ylabel('Distance (m)')
        axs.set_title('Projectile ballistics')

        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles, labels)

        canvas = FigureCanvasTkAgg(figtx, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()
        self.mostrecentfig = figtx

    def tyGraph(self):
        for s in self.rightpanel.grid_slaves():
            s.destroy()

        figty, axs = plt.subplots(1, 1, figsize=(7, 6), dpi=80)
        selected = self.physicshandler.data[self.physicshandler.data['t'] <= self.physicshandler.totalT()]
        axs.plot(selected['t'], selected['y'], '-', linewidth=2, color='b', label='With drag ~ v^2')

        if self.idealset.get():
            iselected = self.idealphysicshandler.data[
                self.idealphysicshandler.data['t'] <= self.idealphysicshandler.totalT()]
            axs.plot(iselected['t'], iselected['y'], '--', linewidth=1, color='g', label='Ideal')

        axs.set_xlabel('Time (s)')
        axs.set_ylabel('Height (m)')
        axs.set_title('Projectile ballistics')

        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles, labels)

        canvas = FigureCanvasTkAgg(figty, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()
        self.mostrecentfig = figty

    def tvGraph(self):
        for s in self.rightpanel.grid_slaves():
            s.destroy()

        figtv, axs = plt.subplots(1, 1, figsize=(7, 6), dpi=80)
        selected = self.physicshandler.data[self.physicshandler.data['t'] <= self.physicshandler.totalT()]
        axs.plot(selected['t'], selected['v'], '-', linewidth=2, color='b', label='With drag ~ v^2')

        if self.idealset.get():
            iselected = self.idealphysicshandler.data[
                self.idealphysicshandler.data['t'] <= self.idealphysicshandler.totalT()]
            axs.plot(iselected['t'], iselected['v'], '--', linewidth=1, color='g', label='Ideal')

        axs.set_xlabel('Time (s)')
        axs.set_ylabel('Velocity (m/s)')
        axs.set_title('Projectile ballistics')

        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles, labels)

        canvas = FigureCanvasTkAgg(figtv, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()
        self.mostrecentfig = figtv

    def xyGraph(self):
        for s in self.rightpanel.grid_slaves():
            s.destroy()

        distance, height = self.geofigbounds()
        figxy, axs = plt.subplots(1, 1, figsize=(7, 6), dpi=80)
        axs.plot(self.physicshandler.data['x'], self.physicshandler.data['y'], '-', linewidth=2, color='b', label='With drag ~ v^2')

        if self.idealset.get():
            axs.plot(self.idealphysicshandler.data['x'], self.idealphysicshandler.data['y'], '--', linewidth=1, color='g', label='Ideal')

        axs.set_xlabel('Distance (m)')
        axs.set_ylabel('Relative Height (m)')

        boundarySize = 2

        if self.idealset.get():
            if self.barrierset.get():
                maxax = np.max([self.physicshandler.totalR() + boundarySize, self.physicshandler.maxH() + boundarySize,
                            self.idealphysicshandler.totalR() + boundarySize,
                                self.idealphysicshandler.maxH() + boundarySize,  distance + boundarySize])
            else:
                maxax = np.max([self.physicshandler.totalR() + boundarySize, self.physicshandler.maxH() + boundarySize,
                                self.idealphysicshandler.totalR() + boundarySize,
                                self.idealphysicshandler.maxH() + boundarySize])
        else:
            if self.barrierset.get():
                maxax = np.max([self.physicshandler.totalR() + boundarySize, self.physicshandler.maxH() + boundarySize,
                                distance + boundarySize])
            else:
                maxax = np.max([self.physicshandler.totalR() + boundarySize, self.physicshandler.maxH() + boundarySize])

        axs.set_xlim(0, maxax)
        axs.set_ylim(np.min([0, height]), maxax)
        axs.set_title('Projectile ballistics')

        if self.barrierset.get():
            axs.axvline(x=distance, linewidth=1, color='red', linestyle='--')
            axs.plot([distance], [height], marker='P', color='green')

        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles, labels)

        canvas = FigureCanvasTkAgg(figxy, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()
        self.mostrecentfig = figxy

    def xvGraph(self):
        for s in self.rightpanel.grid_slaves():
            s.destroy()

        figxv, axs = plt.subplots(1, 1, figsize=(7, 6), dpi=80)
        selected = self.physicshandler.data[self.physicshandler.data['x'] <= self.physicshandler.totalR()]
        axs.plot(selected['x'], selected['v'], '-', linewidth=2, color='b', label='With drag ~ v^2')

        if self.idealset.get():
            iselected = self.idealphysicshandler.data[
                self.idealphysicshandler.data['x'] <= self.idealphysicshandler.totalR()]
            axs.plot(iselected['x'], iselected['v'], '--', linewidth=1, color='g', label='Ideal')

        axs.set_xlabel('Distance (m)')
        axs.set_ylabel('Velocity (m/s)')
        axs.set_title('Projectile ballistics')

        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles, labels)

        canvas = FigureCanvasTkAgg(figxv, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()
        self.mostrecentfig = figxv

    def vyGraph(self):
        for s in self.rightpanel.grid_slaves():
            s.destroy()

        figyv, axs = plt.subplots(1, 1, figsize=(7, 6), dpi=80)
        selected = self.physicshandler.data[self.physicshandler.data['y'] >= self.physicshandler.height]
        axs.plot(selected['y'], selected['v'], '-', linewidth=2, color='b', label='With drag ~ v^2')

        if self.idealset.get():
            iselected = self.idealphysicshandler.data[
                self.idealphysicshandler.data['y'] >= self.idealphysicshandler.height]
            axs.plot(iselected['y'], iselected['v'], '--', linewidth=1, color='g', label='Ideal')

        axs.set_xlabel('Height (m)')
        axs.set_ylabel('Velocity (m/s)')
        axs.set_title('Projectile ballistics')

        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles, labels)

        canvas = FigureCanvasTkAgg(figyv, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()
        self.mostrecentfig = figyv

    def cdvGraph(self):
        for s in self.rightpanel.grid_slaves():
            s.destroy()

        figvre, axs = plt.subplots(1, 1, figsize=(7, 6), dpi=80)
        axs.plot(self.physicshandler.data['v'], self.physicshandler.data['cd'], '-', linewidth=2, color='b', label='With drag ~ v^2')

        if self.idealset.get():
            axs.plot(self.idealphysicshandler.data['y'], self.idealphysicshandler.data['v'], '--', linewidth=1, color='g', label='Ideal')

        axs.set_xlabel('Velocity (m/s)')
        axs.set_ylabel('Cd')
        axs.set_title('Projectile ballistics')

        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles, labels)

        canvas = FigureCanvasTkAgg(figvre, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()
        self.mostrecentfig = figvre

    def cdtGraph(self):
        for s in self.rightpanel.grid_slaves():
            s.destroy()

        figtre, axs = plt.subplots(1, 1, figsize=(7, 6), dpi=80)
        axs.plot(self.physicshandler.data['t'], self.physicshandler.data['cd'], '-', linewidth=2, color='b', label='With drag ~ v^2')

        if self.idealset.get():
            axs.plot(self.idealphysicshandler.data['y'], self.idealphysicshandler.data['v'], '--', linewidth=1, color='g', label='Ideal')

        axs.set_xlabel('Time (s)')
        axs.set_ylabel('Cd')
        axs.set_title('Projectile ballistics')

        handles, labels = axs.get_legend_handles_labels()
        axs.legend(handles, labels)

        canvas = FigureCanvasTkAgg(figtre, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()
        self.mostrecentfig = figtre

    def addStatistics(self):
        stats = tk.LabelFrame(self.rightpanel, text='Results')
        stats.grid(row=1, column=0)

        sphrLabel = tk.Label(stats, text=f'Sphericity: {self.physicshandler.sphericity:.4f}')
        sphrLabel.grid(row=0, column=0)

        surfLabel = tk.Label(stats, text=f'Surface area: {self.physicshandler.surfArea:.4f} m^2')
        surfLabel.grid(row=0, column=1)

        surfLabel = tk.Label(stats, text=f'Volume: {self.physicshandler.sph_volume:.4f} m^3')
        surfLabel.grid(row=1, column=0)

        massLabel = tk.Label(stats, text=f'Mass: {self.physicshandler.m:.4f} kg')
        massLabel.grid(row=1, column=1)

        rangeLabel = tk.Label(stats, text=f'Range: {self.physicshandler.totalR():.1f} m')
        rangeLabel.grid(row=2, column=0)

        rangeLabel = tk.Label(stats, text=f'Max height: {self.physicshandler.maxH():.1f} m')
        rangeLabel.grid(row=2, column=1)

        mheightLabel = tk.Label(stats, text=f'Time to max height: {self.physicshandler.maxT():.1f} s')
        mheightLabel.grid(row=3, column=0)

        mheightLabel = tk.Label(stats, text=f'Time of flight: {self.physicshandler.totalT():.1f} s')
        mheightLabel.grid(row=3, column=1)

        mheightLabel = tk.Label(stats, text=f'Velocity of impact: {self.physicshandler.finalV():.1f} m/s')
        mheightLabel.grid(row=4, column=0)

        mheightLabel = tk.Label(stats, text=f'Angle of impact: {self.physicshandler.finalTheta():.1f} degrees')
        mheightLabel.grid(row=4, column=1)

    def saveCSV(self):
        if self.physicshandler.data is None:
            self.userlabel['text'] = "No computed data exists"
        else:
            fname = filedialog.asksaveasfilename(initialdir = ".", title = "Select file",filetypes = (("CSV files","*.csv"),("all files","*.*")))
            self.physicshandler.save_csv(fname)
            self.userlabel['text'] = "File saved"

    def savePNG(self):
        if self.physicshandler.data is None:
            self.userlabel['text'] = "No computed data exists"
        else:
            fname = filedialog.asksaveasfilename(initialdir=".", title="Select file",
                                                 filetypes=(("PNG files", "*.png"), ("all files", "*.*")))
            self.mostrecentfig.savefig(fname)

    def bye(self):
        self.quit()
        self.destroy()
        sys.exit()


if __name__ == "__main__":
    app = NumericalV2WindThermoGUI()
    app.mainloop()