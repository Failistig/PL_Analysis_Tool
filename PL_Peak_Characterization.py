import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tkinter import Tk, filedialog, Button, Label, Checkbutton, IntVar, Entry, Frame, Toplevel, PhotoImage
from tkinter.simpledialog import askstring
from tkinter import StringVar, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from PIL import Image, ImageTk
from datetime import datetime

class PLAnalysisApp:
    def __init__(self, master):
        self.master = master
        self.master.title("PL Analysis Tool")
        self.master.geometry("1200x700")
        self.master.grid_columnconfigure(1, weight=1)  # Center column expands

        # Bind the window close event to a custom method
        self.master.protocol("WM_DELETE_WINDOW", self.on_close)

        # Variables
        self.file_path = None
        self.data = None
        self.wavelength = None
        self.raw_counts = None
        self.time_labels = None
        self.relative_times = None

        self.plot_spectra_var = IntVar(value=1)
        self.plot_absolute_gradient_var = IntVar(value=1)
        self.plot_relative_change_var = IntVar(value=1)
        self.plot_log_spectra_var = IntVar(value=1)
        self.plot_log_absolute_gradient_var = IntVar(value=1)
        self.plot_log_relative_change_var = IntVar(value=1)
        self.show_legend_var = IntVar(value=1)

        # Save options
        self.save_raw_var = IntVar(value=1)
        self.save_absolute_gradient_var = IntVar(value=1)
        self.save_relative_change_var = IntVar(value=1)
        self.save_log_spectra_var = IntVar(value=1)
        self.save_log_absolute_gradient_var = IntVar(value=1)
        self.save_log_relative_change_var = IntVar(value=1)
        self.save_combined_var = IntVar(value=1)
        self.save_names = {
            "Raw Data": StringVar(value="Raw_Data"),
            "Absolute Gradient": StringVar(value="Absolute_Gradient"),
            "Relative Change": StringVar(value="Relative_Change"),
            "Log Spectra": StringVar(value="Log_Spectra"),
            "Log Absolute Gradient": StringVar(value="Log_Absolute_Gradient"),
            "Log Relative Change": StringVar(value="Log_Relative_Change"),
            "Combined": StringVar(value="Combined_Plots"),
        }

        self.time_checkboxes = []
        self.logo_image = None  # Placeholder for the logo

        self.setup_gui()

    def on_close(self):
        """Ask for confirmation before closing the window."""
        response = messagebox.askyesno(
            "Confirm Exit",
            "Are you sure you want to close? Any unsaved plots will be lost."
        )
        if response:  # If the user presses "Yes"
            self.master.destroy()  # Close the window
        # If the user presses "No", do nothing and keep the window open.

    def setup_gui(self):
        # Title
        Label(self.master, text="PL Analysis Tool", font=("Arial", 16)).grid(row=0, column=0, columnspan=3, pady=10)

        # Load Logo in Top-Right Corner
        self.load_logo("./hzb_logo.jpg")

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
        Checkbutton(self.master, text="Absolute Gradient", variable=self.plot_absolute_gradient_var).grid(row=5, column=0, sticky="w")
        Checkbutton(self.master, text="Relative Change", variable=self.plot_relative_change_var).grid(row=6, column=0, sticky="w")
        Checkbutton(self.master, text="Log Spectra", variable=self.plot_log_spectra_var).grid(row=4, column=1, sticky="w")
        Checkbutton(self.master, text="Log Absolute Gradient", variable=self.plot_log_absolute_gradient_var).grid(row=5, column=1, sticky="w")
        Checkbutton(self.master, text="Log Relative Change", variable=self.plot_log_relative_change_var).grid(row=6, column=1, sticky="w")

        # Info Buttons
        Button(self.master, text="Info: Spectra", width=20, anchor="w", command=lambda: self.show_info("Spectra")).grid(row=4, column=2, sticky="w")
        Button(self.master, text="Info: Absolute Gradient", width=20, anchor="w", command=lambda: self.show_info("Absolute Gradient")).grid(row=5, column=2, sticky="w")
        Button(self.master, text="Info: Relative Change", width=20, anchor="w", command=lambda: self.show_info("Relative Change")).grid(row=6, column=2, sticky="w")
        # Legend Toggle
        Checkbutton(self.master, text="Show Legend", variable=self.show_legend_var).grid(row=7, column=0, sticky="w")

        # Raw Data Selection
        self.time_check_frame = Frame(self.master)
        self.time_check_frame.grid(row=7, column=0, columnspan=3, pady=5)
        Button(self.master, text="Update Raw Plot", width=15, command=self.update_raw_plot).grid(row=8, column=1, pady=5)

        # Plot Frame
        self.plot_frame = Frame(self.master, width=900, height=800, bg="white")
        self.plot_frame.grid(row=9, column=0, columnspan=3, pady=10)

    def load_logo(self, file_path):
        """Loads and places the logo image in the top-right corner."""
        try:
            image = Image.open(file_path)
            image = image.resize((300, 100), Image.Resampling.LANCZOS)  # Resize the logo to 2x size
            self.logo_image = ImageTk.PhotoImage(image)

            # Display the logo
            Label(self.master, image=self.logo_image).grid(row=0, column=2, sticky="ne", padx=10, pady=10)
        except Exception as e:
            print(f"Error loading logo: {e}")

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
        Checkbutton(save_window, text="Absolute Gradient", variable=self.save_absolute_gradient_var).grid(row=2, column=0, sticky="w")
        Checkbutton(save_window, text="Relative Change", variable=self.save_relative_change_var).grid(row=3, column=0, sticky="w")
        Checkbutton(save_window, text="Log Spectra", variable=self.save_log_spectra_var).grid(row=4, column=0, sticky="w")
        Checkbutton(save_window, text="Log Absolute Gradient", variable=self.save_log_absolute_gradient_var).grid(row=5, column=0, sticky="w")
        Checkbutton(save_window, text="Log Relative Change", variable=self.save_log_relative_change_var).grid(row=6, column=0, sticky="w")
        Checkbutton(save_window, text="Combined Plots", variable=self.save_combined_var).grid(row=7, column=0, sticky="w")

        # Entry fields for file names
        for i, (plot_name, var) in enumerate(self.save_names.items()):
            Label(save_window, text=f"{plot_name} Name:").grid(row=i + 1, column=1, sticky="e", padx=10)
            Entry(save_window, textvariable=var, width=20).grid(row=i + 1, column=2, sticky="w")

        Button(save_window, text="Save", command=lambda: self.save_plots(save_window)).grid(row=8, column=1, columnspan=2, pady=10)

    def save_plots(self, save_window):
        save_window.destroy()
        save_dir = filedialog.askdirectory(title="Select Save Directory")
        if not save_dir:
            return

        if self.save_raw_var.get():
            self.save_plot(self.plot_spectra, save_dir, self.save_names["Raw Data"].get())
        if self.save_absolute_gradient_var.get():
            self.save_plot(self.plot_absolute_gradient, save_dir, self.save_names["Absolute Gradient"].get())
        if self.save_relative_change_var.get():
            self.save_plot(self.plot_relative_change, save_dir, self.save_names["Relative Change"].get())
        if self.save_log_spectra_var.get():
            self.save_plot(self.plot_log_spectra, save_dir, self.save_names["Log Spectra"].get())
        if self.save_log_absolute_gradient_var.get():
            self.save_plot(self.plot_log_absolute_gradient, save_dir, self.save_names["Log Absolute Gradient"].get())
        if self.save_log_relative_change_var.get():
            self.save_plot(self.plot_log_relative_change, save_dir, self.save_names["Log Relative Change"].get())
        if self.save_combined_var.get():
            self.save_combined_plot(save_dir)

    def save_plot(self, plot_function, save_dir, plot_name):
        fig, ax = plt.subplots(figsize=(8, 5))
        plot_function(ax)
        fig.savefig(f"{save_dir}/{plot_name}.png")
        plt.close(fig)
        print(f"Saved {plot_name}")

    def save_combined_plot(self, save_dir):
        fig, axes = plt.subplots(2, 3, figsize=(14, 6), constrained_layout=True)
        self.plot_spectra(axes[0, 0])
        self.plot_absolute_gradient(axes[0, 1])
        self.plot_relative_change(axes[0, 2])
        self.plot_log_spectra(axes[1, 0])
        self.plot_log_absolute_gradient(axes[1, 1])
        self.plot_log_relative_change(axes[1, 2])
        fig.tight_layout()
        fig.savefig(f"{save_dir}/{self.save_names['Combined'].get()}.png")
        plt.close(fig)
        print("Saved combined plots.")

    def show_info(self, plot_type):
        info = {
            "Spectra": (
                "Spectra:\n"
                "The spectra plot shows the photoluminescence (PL) signal intensity as a function "
                "of wavelength for different time intervals. The logarithmic plot shows the same data on a logarithmic scale.\n\n"
                "The raw data is taken from the PL measurements, where each curve corresponds "
                "to a specific time interval.\n\n"
                "X-axis: Wavelength (nm)\nY-axis: Intensity (a.u.)"
            ),
            "Absolute Gradient": (
                "Absolute Intensity Gradient:\n"
                "This plot represents the absolute gradient calculated from the PL spectra. "
                "It shows the variation of PL intensity over time. The logarithmic plot shows the same data on a logarithmic scale.\n\n"
                "It is computed as the numerical gradient of the raw photoluminescence (PL) counts "
                "along the time axis (row-wise) and measures how the PL intensity changes with respect "
                "to time for each wavelength.\n\n"
                "Then the absolute value of the gradient is taken for each wavelength, and the mean "
                "across all time intervals is calculated. This results in a single average gradient value "
                "for each wavelength, representing how unstable or variable the PL intensity is over time "
                "at that wavelength.\n\n"
                "X-axis: Time (s)\nY-axis: Absolute Gradient (a.u.)"
            ),
            "Relative Change": (
                "Relative Intensity Change:\n"
                "This plot tracks relative changes in PL intensity. The logarithmic plot shows the same data on a logarithmic scale.\n\n"
                "It is computed as the absolute differences in PL intensity between consecutive time points "
                "and represents the magnitude of change in intensity over time at each wavelength.\n\n"
                "The change is normalized by the intensity at the previous time step to obtain a relative metric.\n\n"
                "The mean of the relative changes is then calculated across all time intervals for each wavelength "
                "providing an average metric for each wavelength.\n\n"
                "X-axis: Time (s)\nY-axis: Relative Change Intensity (a.u.)"
            ),
        }
        messagebox.showinfo(f"Info: {plot_type}", info.get(plot_type, "No information available."))

    def plot_in_window(self):
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        fig, axes = plt.subplots(2, 3, figsize=(14, 6), constrained_layout=True)

        # Clear all axes
        for ax in axes.flatten():
            ax.clear()

        # Plot in fixed positions
        if self.plot_spectra_var.get():
            self.plot_spectra(axes[0, 0])
        else:
            axes[0, 0].axis('off')

        if self.plot_absolute_gradient_var.get():
            self.plot_absolute_gradient(axes[0, 1])
        else:
            axes[0, 1].axis('off')

        if self.plot_relative_change_var.get():
            self.plot_relative_change(axes[0, 2])
        else:
            axes[0, 2].axis('off')

        if self.plot_log_spectra_var.get():
            self.plot_log_spectra(axes[1, 0])
        else:
            axes[1, 0].axis('off')

        if self.plot_log_absolute_gradient_var.get():
            self.plot_log_absolute_gradient(axes[1, 1])
        else:
            axes[1, 1].axis('off')

        if self.plot_log_relative_change_var.get():
            self.plot_log_relative_change(axes[1, 2])
        else:
            axes[1, 2].axis('off')

        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        toolbar = NavigationToolbar2Tk(canvas, self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack()

    def plot_spectra(self, ax):
        selected_indices = [i for var, i in self.time_checkboxes if var.get()]
        colors = plt.cm.viridis(np.linspace(0, 1, len(selected_indices)))
        for idx, i in enumerate(selected_indices):
            ax.plot(self.wavelength, self.raw_counts[:, i], color=colors[idx], label=self.relative_times[i])
        if self.show_legend_var.get():
            ax.legend()
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Counts")
        ax.set_title("Raw PL Spectra")

    def plot_absolute_gradient(self, ax):
        gradient = np.gradient(self.raw_counts, axis=0)
        mean_gradient = np.mean(np.abs(gradient), axis=1)
        ax.plot(self.wavelength, mean_gradient, color="blue", label="Absolute Gradient")
        ax.set_title("Absolute Gradient")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Gradient Value")
        if self.show_legend_var.get():
            ax.legend()

    def plot_relative_change(self, ax):
        shape_change = np.abs(np.diff(self.raw_counts, axis=1))
        epsilon = 1e-8
        quotient = shape_change / (np.abs(self.raw_counts[:, :-1]) + epsilon)
        average_quotient = np.mean(quotient, axis=1)
        ax.plot(self.wavelength[:-1], average_quotient[:len(self.wavelength) - 1], color="green",
                label="Relative Change")
        ax.set_title("Relative Change")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Relative Change Value")
        if self.show_legend_var.get():
            ax.legend()

    def plot_log_spectra(self, ax):
        selected_indices = [i for var, i in self.time_checkboxes if var.get()]
        colors = plt.cm.viridis(np.linspace(0, 1, len(selected_indices)))
        for idx, i in enumerate(selected_indices):
            ax.plot(self.wavelength, self.raw_counts[:, i], color=colors[idx], label=self.relative_times[i])
        if self.show_legend_var.get():
            ax.legend()
        ax.set_yscale('log')
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Counts")
        ax.set_title("Raw PL Spectra (log)")

    def plot_log_absolute_gradient(self, ax):
        gradient = np.gradient(self.raw_counts, axis=0)
        mean_gradient = np.mean(np.abs(gradient), axis=1)
        ax.plot(self.wavelength, mean_gradient, color="blue", label="Absolute Gradient")
        ax.set_title("Absolute Gradient")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Gradient Value (log)")
        ax.set_yscale('log')
        if self.show_legend_var.get():
            ax.legend()

    def plot_log_relative_change(self, ax):
        shape_change = np.abs(np.diff(self.raw_counts, axis=1))
        epsilon = 1e-8
        quotient = shape_change / (np.abs(self.raw_counts[:, :-1]) + epsilon)
        average_quotient = np.mean(quotient, axis=1)
        ax.plot(self.wavelength[:-1], average_quotient[:len(self.wavelength) - 1], color="green",
                label="Relative Change")
        ax.set_title("Relative Change")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Relative Change Value (log)")
        ax.set_yscale('log')
        if self.show_legend_var.get():
            ax.legend()

    def update_raw_plot(self):
        self.plot_in_window()

if __name__ == "__main__":
    root = Tk()
    app = PLAnalysisApp(root)
    root.mainloop()