# Red Ciudadana de Estaciones Meteorologicas
#
# Copyright @ 2021
#
# Authors: Santiago Nunez-Corrales <snunezcr@gmail.com>
#          Jose Brenes-Andre <jbrenes54@gmail.com>
import numpy as np
import matplotlib
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
from matplotlib.figure import Figure
from balistica.PhysicsEngine.IdealPhysicsHandler import IdealPhysicsHandler
from tkinter import filedialog


class IdealGUI(tk.Frame):
    def __init__(self, master=None):
        self.physicshandler = IdealPhysicsHandler(0, 0)

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

        # Control for angle
        self.anglelable = tk.Label(self.ulpanel, text='Initial angle (degrees)')
        self.anglelable.grid(row=0, column=0)
        self.angleinput = tk.Scale(self.ulpanel, from_=0, to=90, resolution=1, length=170,orient=tk.HORIZONTAL)
        self.angleinput.grid(row=0, column=1)

        # Control for velocity
        self.velocitylabel = tk.Label(self.ulpanel, text='Initial velocity (m/s)')
        self.velocitylabel.grid(row=1, column=0)
        self.velocityinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.velocityinput.grid(row=1, column=1)

        self.velocityinput.insert(0, '125')

        self.latIlabel = tk.Label(self.ulpanel, text='I. Lat (m)')
        self.latIlabel.grid(row=2, column=0)
        self.lonIlabel = tk.Label(self.ulpanel, text='I. Lon (m)')
        self.lonIlabel.grid(row=2, column=1)
        self.heightIlabel = tk.Label(self.ulpanel, text='I. Height (m)')
        self.heightIlabel.grid(row=2, column=2)

        self.latIinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.latIinput.grid(row=3, column=0)
        self.lonIinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.lonIinput.grid(row=3, column=1)
        self.heightIinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.heightIinput.grid(row=3, column=2)

        self.latIinput.insert(0, '0')
        self.lonIinput.insert(0, '0')
        self.heightIinput.insert(0, '0')

        self.pblanklabel = tk.Label(self.ulpanel, text='')
        self.pblanklabel.grid(row=4, column=0, columnspan=2)

        self.latFlabel = tk.Label(self.ulpanel, text='F. Lat (m)')
        self.latFlabel.grid(row=5, column=0)
        self.lonFlabel = tk.Label(self.ulpanel, text='F. Lon (m)')
        self.lonFlabel.grid(row=5, column=1)
        self.heightFlabel = tk.Label(self.ulpanel, text='F. Height (m)')
        self.heightFlabel.grid(row=5, column=2)

        self.latFinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.latFinput.grid(row=6, column=0)
        self.lonFinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.lonFinput.grid(row=6, column=1)
        self.heightFinput = tk.Entry(self.ulpanel, justify=tk.RIGHT, width=10)
        self.heightFinput.grid(row=6, column=2)

        self.latFinput.insert(0, '100')
        self.lonFinput.insert(0, '100')
        self.heightFinput.insert(0, '0')

        self.barrierset = tk.BooleanVar()
        self.barriercheck = tk.Checkbutton(self.ulpanel, justify=tk.RIGHT, variable=self.barrierset, onvalue=True,
                                           offvalue=False, text='Show barrier')
        self.barriercheck.grid(row=7, column=0)

        # Controls grid for upper left pannel
        self.blpanel = tk.Frame(self.leftpanel)
        self.blpanel.pack(side=tk.BOTTOM)

        # Buttons for various functions
        # Buttons for various functions
        self.blanklabel= tk.Label(self.blpanel, text="")
        self.blanklabel.grid(row=0, column=0, columnspan=2)

        self.computebutton = tk.Button(self.blpanel, text="Compute", width=20, command=self.compute, default=tk.NORMAL)
        self.computebutton.grid(row=1, column=0, columnspan=3)

        self.computebutton = tk.Button(self.blpanel, text="x(t) vs. t", width=10, command=self.txGraph, default=tk.NORMAL)
        self.computebutton.grid(row=2, column=0)

        self.computebutton = tk.Button(self.blpanel, text="z(t) vs. t", width=10, command=self.tyGraph, default=tk.NORMAL)
        self.computebutton.grid(row=2, column=1)

        self.computebutton = tk.Button(self.blpanel, text="v(t) vs. t", width=10, command=self.tvGraph, default=tk.NORMAL)
        self.computebutton.grid(row=2, column=2)

        self.computebutton = tk.Button(self.blpanel, text="z(t) vs. x(t)", width=10, command=self.xyGraph, default=tk.NORMAL)
        self.computebutton.grid(row=3, column=0)

        self.computebutton = tk.Button(self.blpanel, text="v(t) vs. x(t)", width=10, command=self.xvGraph, default=tk.NORMAL)
        self.computebutton.grid(row=3, column=1)

        self.computebutton = tk.Button(self.blpanel, text="v(t) vs. z(t)", width=10, command=self.yvGraph, default=tk.NORMAL)
        self.computebutton.grid(row=3, column=2)

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
        axs.set_title('Ballistics - Ideal case')
        canvas = FigureCanvasTkAgg(fig, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()

        self.mostrecentfig : Figure = fig

    def geography(self):
        latI = 0.0
        try:
            latI = float(self.latIinput.get())
        except:
            self.userlabel['text'] = "Initial latitude format incorrect"

        latF = 0.0
        try:
            latF = float(self.latFinput.get())
        except:
            self.userlabel['text'] = "Final latitude format incorrect"

        lonI = 0.0
        try:
            lonI = float(self.lonIinput.get())
        except:
            self.userlabel['text'] = "Initial longitude format incorrect"

        lonF = 0.0
        try:
            lonF = float(self.lonFinput.get())
        except:
            self.userlabel['text'] = "Final longitude format incorrect"

        heightI = 0.0
        try:
            heightI = float(self.heightIinput.get())
        except:
            self.userlabel['text'] = "Initial latitude format incorrect"

        heightF = 0.0
        try:
            heightF = float(self.heightFinput.get())
        except:
            self.userlabel['text'] = "Initial latitude format incorrect"

        distance = np.sqrt(np.power((latF - latI), 2) + np.power((lonF - lonI), 2))
        height = heightF - heightI

        return (distance, height)

    def compute(self):
        self.userlabel['text'] = ""

        vel0 = 0.0
        try:
            vel0 = float(self.velocityinput.get())
        except:
            self.userlabel['text'] = "Velocity format incorrect"
            return

        theta = np.deg2rad(float(self.angleinput.get()))

        self.physicshandler.v0 = vel0
        self.physicshandler.theta = theta

        distance, height = self.geography()

        self.physicshandler.distance = distance

        if self.barrierset.get():
            self.physicshandler.height = height
            self.physicshandler.barrier = True
        else:
            self.physicshandler.height = -1
            self.physicshandler.barrier = False

        self.physicshandler.compute()
        self.xyGraph()

    def txGraph(self):
        for s in self.rightpanel.grid_slaves():
            s.destroy()

        figtx, axs = plt.subplots(1, 1, figsize=(7, 6), dpi=80)
        selected = self.physicshandler.data[self.physicshandler.data['t'] <= self.physicshandler.totalT()]
        axs.plot(selected['t'], selected['x'], '-', linewidth=2, color='b')
        axs.set_xlabel('Time (s)')
        axs.set_ylabel('Distance (m)')
        axs.set_title('Ballistics with constant drag (b) proportional to v')
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
        axs.plot(selected['t'], selected['z'], '-', linewidth=2, color='b')
        axs.set_xlabel('Time (s)')
        axs.set_ylabel('Height (m)')
        axs.set_title('Ballistics with drag (b) proportional to v')
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
        axs.plot(selected['t'], selected['v'], '-', linewidth=2, color='b')
        axs.set_xlabel('Time (s)')
        axs.set_ylabel('Velocity (m/s)')
        axs.set_title('Ballistics with drag (b) proportional to v')
        canvas = FigureCanvasTkAgg(figtv, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()
        self.mostrecentfig = figtv

    def xyGraph(self):
        for s in self.rightpanel.grid_slaves():
            s.destroy()

        distance, height = self.geography()
        figxy, axs = plt.subplots(1, 1, figsize=(7, 6), dpi=80)
        axs.plot(self.physicshandler.data['x'], self.physicshandler.data['z'], '-', linewidth=2, color='b')
        axs.set_xlabel('Distance (m)')
        axs.set_ylabel('Height (m)')

        if self.barrierset.get():
            maxax = np.max([self.physicshandler.totalR() + 10, self.physicshandler.maxH() + 10, distance + 20])
            minay = np.min([0, self.physicshandler.height - 10])
        else:
            maxax = np.max([self.physicshandler.totalR() + 10, self.physicshandler.maxH() + 10])
            minay = 0

        axs.set_xlim(np.min([0, self.physicshandler.totalR()]), maxax)
        axs.set_ylim(minay, maxax)

        axs.set_title('Ballistics with drag (b) proportional to v')

        if self.barrierset.get():
            axs.axvline(x=distance, color='red', linestyle='--')
            axs.plot([distance], [height], marker='P', color='green')

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
        axs.plot(selected['x'], selected['v'], '-', linewidth=2, color='b')
        axs.set_xlabel('Distance (m)')
        axs.set_ylabel('Velocity (m/s)')
        axs.set_title('Ballistics with drag (b) proportional to v')
        canvas = FigureCanvasTkAgg(figxv, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()
        self.mostrecentfig = figxv

    def yvGraph(self):
        for s in self.rightpanel.grid_slaves():
            s.destroy()

        figyv, axs = plt.subplots(1, 1, figsize=(7, 6), dpi=80)
        selected = self.physicshandler.data[self.physicshandler.data['z'] >= self.physicshandler.height]
        axs.plot(selected['z'], selected['v'], '-', linewidth=2, color='b')
        axs.set_xlabel('Height (m)')
        axs.set_ylabel('Velocity (m/s)')
        axs.set_title('Ballistics with drag (b) proportional to v')
        axs.invert_xaxis()
        canvas = FigureCanvasTkAgg(figyv, master=self.rightpanel)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0)
        self.addStatistics()
        self.mostrecentfig = figyv

    def addStatistics(self):
        stats = tk.LabelFrame(self.rightpanel, text='Results')
        stats.grid(row=1, column=0)

        rangeLabel = tk.Label(stats, text=f'Range: {self.physicshandler.totalR():.1f} m')
        rangeLabel.grid(row=0, column=0)

        rangeLabel = tk.Label(stats, text=f'Max height: {self.physicshandler.maxH():.1f} m')
        rangeLabel.grid(row=1, column=0)

        mheightLabel = tk.Label(stats, text=f'Time to max height: {self.physicshandler.maxT():.1f} s')
        mheightLabel.grid(row=2, column=0)

        mheightLabel = tk.Label(stats, text=f'Time of flight: {self.physicshandler.totalT():.1f} s')
        mheightLabel.grid(row=3, column=0)

        mheightLabel = tk.Label(stats, text=f'Velocity of impact: {self.physicshandler.finalV():.1f} m/s')
        mheightLabel.grid(row=4, column=0)

        mheightLabel = tk.Label(stats, text=f'Angle of impact: {self.physicshandler.finalTheta():.1f} degrees')
        mheightLabel.grid(row=5, column=0)


    def saveCSV(self):
        if self.physicshandler.data is None:
            self.userlabel['text'] = "No computed data exists"
        else:
            fname = filedialog.asksaveasfilename(initialdir = ".", title = "Select file",filetypes = (("CSV files","*.csv"),("all files","*.*")))
            self.physicshandler.save_csv(fname+".csv")
            self.userlabel['text'] = "File saved"

    def savePNG(self):
        if self.physicshandler.data is None:
            self.userlabel['text'] = "No computed data exists"
        else:
            fname = filedialog.asksaveasfilename(initialdir=".", title="Select file",
                                                 filetypes=(("PNG files", "*.png"), ("all files", "*.*")))
            self.mostrecentfig.savefig(fname+".png")
            self.userlabel['text'] = "File saved"

    def bye(self):
        self.quit()
        self.destroy()


if __name__ == "__main__":
    app = IdealGUI()
    app.mainloop()