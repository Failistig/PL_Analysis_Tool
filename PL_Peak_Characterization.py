import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tkinter import Tk, filedialog, Button, Label, Checkbutton, IntVar, Entry, Frame, Toplevel
from tkinter.simpledialog import askstring
from tkinter import StringVar
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from datetime import datetime


class PLAnalysisApp:
    def __init__(self, master):
        self.master = master
        self.master.title("PL Analysis Tool")
        self.master.geometry("1200x700")
        self.master.grid_columnconfigure(1, weight=1)  # Center column expands

        # Variables
        self.file_path = None
        self.data = None
        self.wavelength = None
        self.raw_counts = None
        self.time_labels = None
        self.relative_times = None

        self.plot_spectra_var = IntVar(value=1)
        self.plot_instability_var = IntVar(value=1)
        self.plot_segregation_var = IntVar(value=1)

        # Save options
        self.save_raw_var = IntVar(value=1)
        self.save_instability_var = IntVar(value=1)
        self.save_segregation_var = IntVar(value=1)
        self.save_combined_var = IntVar(value=1)
        self.save_names = {
            "Raw Data": StringVar(value="Raw_Data"),
            "Instability Gradient": StringVar(value="Instability_Gradient"),
            "Halide Segregation": StringVar(value="Halide_Segregation"),
            "Combined": StringVar(value="Combined_Plots"),
        }

        self.time_checkboxes = []
        self.setup_gui()

    def setup_gui(self):
        # Title
        Label(self.master, text="PL Analysis Tool", font=("Arial", 16)).grid(row=0, column=0, columnspan=3, pady=10)

        # File Load Section
        Button(self.master, text="Load File", width=15, command=self.load_file).grid(row=1, column=0, padx=10, pady=5, sticky="w")
        self.file_path_entry = Entry(self.master, width=80, state='disabled')
        self.file_path_entry.grid(row=1, column=1, columnspan=2, pady=5, sticky="w")

        # Buttons Section
        Button(self.master, text="Save Plots", width=15, command=self.save_plots_dialog).grid(row=2, column=0, padx=10, pady=5, sticky="w")
        Button(self.master, text="Plot All", width=15, command=self.plot_in_window).grid(row=2, column=1, pady=5, sticky="w")

        # Plot Selection
        Label(self.master, text="Select Plots:").grid(row=3, column=0, pady=5, sticky="w")
        Checkbutton(self.master, text="Spectra", variable=self.plot_spectra_var).grid(row=4, column=0, sticky="w")
        Checkbutton(self.master, text="Instability Gradient", variable=self.plot_instability_var).grid(row=5, column=0, sticky="w")
        Checkbutton(self.master, text="Halide Segregation", variable=self.plot_segregation_var).grid(row=6, column=0, sticky="w")

        # Info Buttons
        Button(self.master, text="Info: Spectra", width=15, command=lambda: self.show_info("Spectra")).grid(row=4, column=1, sticky="w")
        Button(self.master, text="Info: Instability", width=15, command=lambda: self.show_info("Instability")).grid(row=5, column=1, sticky="w")
        Button(self.master, text="Info: Segregation", width=15, command=lambda: self.show_info("Segregation")).grid(row=6, column=1, sticky="w")

        # Plot Frame
        self.plot_frame = Frame(self.master, width=900, height=400, bg="white")
        self.plot_frame.grid(row=7, column=0, columnspan=3, pady=10)

        # Raw Data Selection
        self.time_check_frame = Frame(self.master)
        self.time_check_frame.grid(row=8, column=0, columnspan=3, pady=5)
        Button(self.master, text="Update Raw Plot", width=15, command=self.update_raw_plot).grid(row=9, column=1, pady=5)

    def load_file(self):
        self.file_path = filedialog.askopenfilename(title="Select Data File", filetypes=[("Text Files", "*.txt")])
        if not self.file_path:
            return

        self.file_path_entry.config(state='normal')
        self.file_path_entry.delete(0, 'end')
        self.file_path_entry.insert(0, self.file_path)
        self.file_path_entry.config(state='disabled')

        self.data = pd.read_csv(self.file_path, sep="\t", skiprows=31)
        self.wavelength = self.data.iloc[:, 0].values
        self.raw_counts = self.data.iloc[:, 1:].values

        with open(self.file_path, 'r') as file:
            lines = file.readlines()
        raw_time_labels = lines[0].strip().split("\t")[1:]
        self.relative_times = self.calculate_relative_times(raw_time_labels)
        self.setup_time_checkboxes()

    def calculate_relative_times(self, raw_times):
        try:
            time_format = "%H:%M:%S"
            first_time = datetime.strptime(raw_times[0], time_format)
            relative_times = []
            for time_str in raw_times:
                current_time = datetime.strptime(time_str, time_format)
                delta_seconds = int((current_time - first_time).total_seconds())
                relative_times.append(f"{delta_seconds}s")
            return relative_times
        except Exception as e:
            print(f"Error parsing times: {e}")
            return None

    def setup_time_checkboxes(self):
        for widget in self.time_check_frame.winfo_children():
            widget.destroy()
        self.time_checkboxes = []

        for i, time_label in enumerate(self.relative_times):
            var = IntVar(value=1)
            cb = Checkbutton(self.time_check_frame, text=time_label, variable=var)
            cb.pack(side='left', padx=5)
            self.time_checkboxes.append((var, i))

    def save_plots_dialog(self):
        save_window = Toplevel(self.master)
        save_window.title("Save Options")

        Label(save_window, text="Select Plots to Save:").grid(row=0, column=0, columnspan=2, pady=10)

        # Checkboxes for save options
        Checkbutton(save_window, text="Raw Data", variable=self.save_raw_var).grid(row=1, column=0, sticky="w")
        Checkbutton(save_window, text="Instability Gradient", variable=self.save_instability_var).grid(row=2, column=0, sticky="w")
        Checkbutton(save_window, text="Halide Segregation", variable=self.save_segregation_var).grid(row=3, column=0, sticky="w")
        Checkbutton(save_window, text="Combined Plots", variable=self.save_combined_var).grid(row=4, column=0, sticky="w")

        # Entry fields for file names
        for i, (plot_name, var) in enumerate(self.save_names.items()):
            Label(save_window, text=f"{plot_name} Name:").grid(row=i + 1, column=1, sticky="e", padx=10)
            Entry(save_window, textvariable=var, width=20).grid(row=i + 1, column=2, sticky="w")

        Button(save_window, text="Save", command=lambda: self.save_plots(save_window)).grid(row=5, column=1, columnspan=2, pady=10)

    def save_plots(self, save_window):
        save_window.destroy()
        save_dir = filedialog.askdirectory(title="Select Save Directory")
        if not save_dir:
            return

        if self.save_raw_var.get():
            self.save_plot(self.plot_spectra, save_dir, self.save_names["Raw Data"].get())
        if self.save_instability_var.get():
            self.save_plot(self.plot_instability, save_dir, self.save_names["Instability Gradient"].get())
        if self.save_segregation_var.get():
            self.save_plot(self.plot_segregation, save_dir, self.save_names["Halide Segregation"].get())
        if self.save_combined_var.get():
            self.save_combined_plot(save_dir)

    def save_plot(self, plot_function, save_dir, plot_name):
        fig, ax = plt.subplots(figsize=(8, 5))
        plot_function(ax)
        fig.savefig(f"{save_dir}/{plot_name}.png")
        plt.close(fig)
        print(f"Saved {plot_name}")

    def save_combined_plot(self, save_dir):
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        self.plot_spectra(axes[0])
        self.plot_instability(axes[1])
        self.plot_segregation(axes[2])
        fig.tight_layout()
        fig.savefig(f"{save_dir}/{self.save_names['Combined'].get()}.png")
        plt.close(fig)
        print("Saved combined plots.")

    def show_info(self, plot_type):
        info = {
            "Spectra": (
                "Spectra:\n"
                "The spectra plot shows the photoluminescence (PL) signal intensity as a function "
                "of wavelength for different time intervals.\n\n"
                "The raw data is taken from the PL measurements, where each curve corresponds "
                "to a specific time interval.\n\n"
                "X-axis: Wavelength (nm)\nY-axis: Intensity (counts)."
            ),
            "Instability": (
                "Instability Gradient:\n"
                "The instability gradient is calculated using the gradient (rate of change) "
                "of the PL signal intensity over time for each wavelength.\n\n"
                "This highlights how unstable or dynamic the PL signal is across wavelengths "
                "over the measured time intervals.\n\n"
                "X-axis: Wavelength (nm)\nY-axis: Mean Absolute Gradient."
            ),
            "Segregation": (
                "Halide Segregation:\n"
                "Halide segregation measures the changes in PL signal intensity between adjacent "
                "time intervals. It reflects changes in material composition or structure over time.\n\n"
                "This is calculated as the absolute difference between consecutive PL spectra.\n\n"
                "X-axis: Wavelength (nm)\nY-axis: Mean Change in Counts."
            ),
        }

        top = Toplevel(self.master)
        top.title(f"{plot_type} Information")
        Label(top, text=info.get(plot_type, "No information available."), wraplength=500, justify="left").pack(padx=20, pady=20)

    def plot_in_window(self):
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        if self.plot_spectra_var.get():
            self.plot_spectra(axes[0])
        if self.plot_instability_var.get():
            self.plot_instability(axes[1])
        if self.plot_segregation_var.get():
            self.plot_segregation(axes[2])

        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        toolbar = NavigationToolbar2Tk(canvas, self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack()

    def plot_spectra(self, ax):
        selected_indices = [i for var, i in self.time_checkboxes if var.get()]
        colors = plt.cm.viridis(np.linspace(0, 1, len(selected_indices)))
        for idx, i in enumerate(selected_indices):
            ax.plot(self.wavelength, self.raw_counts[:, i], color=colors[idx], label=self.relative_times[i])
        ax.legend()
        ax.set_title("Spectra")

    def plot_instability(self, ax):
        gradient = np.gradient(self.raw_counts, axis=0)
        mean_gradient = np.mean(np.abs(gradient), axis=1)
        ax.plot(self.wavelength, mean_gradient, color="red")
        ax.set_title("Instability Gradient")

    def plot_segregation(self, ax):
        shape_change = np.abs(np.diff(self.raw_counts, axis=0))
        mean_segregation = np.mean(shape_change, axis=1)
        ax.plot(self.wavelength[:-1], mean_segregation, color="green")
        ax.set_title("Halide Segregation")

    def update_raw_plot(self):
        self.plot_in_window()


if __name__ == "__main__":
    root = Tk()
    app = PLAnalysisApp(root)
    root.mainloop()
