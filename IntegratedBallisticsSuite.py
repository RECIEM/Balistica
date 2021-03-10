# Red Ciudadana de Estaciones Meteorologicas
#
# Copyright @ 2019
#
# Author: Santiago Nunez-Corrales <snunezcr@gmail.com>
import tkinter as tk
from tkinter import messagebox as mb
from tkinterhtml import HtmlFrame
from GUI.Ideal import IdealGUI
from GUI.AnalyticV import AnalyticVGUI
from GUI.NumericalV import NumericalVGUI
from GUI.NumericalV2 import NumericalV2GUI
from GUI.NumericalVWind import NumericalVWindGUI
from GUI.NumericalV2Wind import NumericalV2WindGUI
from GUI.NumericalCombined import NumericalV2WindExtGUI
from GUI.NumericalThermoFluid import NumericalV2WindThermoGUI


class IntegratedGUI(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("RECIEM - Integrated Ballistics Suite")
        self.geometry("1100x700")
        menubar = tk.Menu(master=self)
        self.config(menu=menubar)

        modelsmenu = tk.Menu(menubar, tearoff=0)
        modelsmenu.add_command(label='Ideal', command=self.add_ideal_model)
        modelsmenu.add_separator()
        modelsmenu.add_command(label='Analytic ~ velocity', command=self.add_anlyt_v)
        modelsmenu.add_command(label='Numerical ~ velocity', command=self.add_nmrcl_v)
        modelsmenu.add_command(label='Numerical ~ velocity squared', command=self.add_nmrcl_v2)
        modelsmenu.add_separator()
        modelsmenu.add_command(label='Numerical ~ velocity + wind', command=self.add_nmrcl_v_wind)
        modelsmenu.add_command(label='Numerical ~ velocity sq. + wind', command=self.add_nmrcl_v2_wind)
        modelsmenu.add_command(label='Numerical ~ velocity sq. + wind + density and sphericity',
                               command=self.add_nmrcl_v2_wind_ext)
        modelsmenu.add_command(label='Numerical ~ vel. sq + wind + atmospherics',
                               command=self.add_nmrcl_v2_wind_ext_thermo)
        menubar.add_cascade(label="Models", menu=modelsmenu)

        helpmenu = tk.Menu(menubar, tearoff=0)
        helpmenu.add_command(label='Usage instructions', command=self.show_instructions)
        helpmenu.add_command(label='About this program', command=self.show_about)
        menubar.add_cascade(label="Help", menu=helpmenu)

        self.most_recent_frame = None
        self.images = {}

    def onExit(self):
        self.quit()

    def remove_prior(self):
        if self.most_recent_frame is not None:
            self.most_recent_frame.pack_forget()
            self.most_recent_frame.destroy()
        else:
            pass

    def add_ideal_model(self):
        self.remove_prior()
        self.most_recent_frame = IdealGUI(self)
        self.most_recent_frame.pack()

    def add_anlyt_v(self):
        self.remove_prior()
        self.most_recent_frame = AnalyticVGUI(self)
        self.most_recent_frame.pack()

    def add_nmrcl_v(self):
        self.remove_prior()
        self.most_recent_frame = NumericalVGUI(self)
        self.most_recent_frame.pack()

    def add_nmrcl_v2(self):
        self.remove_prior()
        self.most_recent_frame = NumericalV2GUI(self)
        self.most_recent_frame.pack()

    def add_nmrcl_v_wind(self):
        self.remove_prior()
        self.most_recent_frame = NumericalVWindGUI(self)
        self.most_recent_frame.pack()

    def add_nmrcl_v2_wind(self):
        self.remove_prior()
        self.most_recent_frame = NumericalV2WindGUI(self)
        self.most_recent_frame.pack()

    def add_nmrcl_v2_wind_ext(self):
        self.remove_prior()
        self.most_recent_frame = NumericalV2WindExtGUI(self)
        self.most_recent_frame.pack()

    def add_nmrcl_v2_wind_ext_thermo(self):
        self.remove_prior()
        self.most_recent_frame = NumericalV2WindThermoGUI(self)
        self.most_recent_frame.pack()

    def view_plots(self):
        pass

    def compare_plots(self):
        pass

    def show_instructions(self):
        help_text = '''
        <html>
        <body>
        
        <h1> RECIEM - Integrated Ballistics Simulation Suite </h1>
        
        <p>Jose Brenes-Andre</p>
        <p>Copyright @ 2019</p>
        
        <h2> Generalities </h2>
        
        <h2> Ideal model </h2>
        
        <h2> Analytic model - proportional to velocity </h2>
        
        <h2> Numerical model - proportional to velocity </h2>
        
        <h2> Numerical model - proportional to velocity squared </h2>
        
        <h2> Numerical model - proportional to velocity + wind </h2>
        
        <h2> Numerical model - proportional to velocity squared + wind </h2>
        
        <h2> Numerical model - proportional to velocity squared + wind + sphericity </h2>
        
        <h2> Numerical model - proportional to velocity squared + wind + thermodynamics </h2>
        
        </body>
        </html>
        '''
        new = tk.Toplevel()
        new.wm_title('Help')
        new.wm_geometry('800x400')
        frame = HtmlFrame(new, horizontal_scrollbar="auto")
        frame.grid(sticky=tk.NSEW)
        frame.set_content(help_text)
        frame.pack()

    def show_about(self):
        message = '''
        RECIEM Costa Rica
        All rights reserved, 2020
        '''
        mb.showinfo("About us", message)

if __name__ == "__main__":
    app = IntegratedGUI()
    app.mainloop()
