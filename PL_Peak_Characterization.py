import tkinter as tk
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tkinter import Tk, filedialog, Button, Label, Checkbutton, IntVar, Entry, Frame, Toplevel, StringVar, messagebox, Text, Scrollbar, Canvas, OptionMenu
from tkinter.simpledialog import askstring
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from PIL import Image, ImageTk
from datetime import datetime
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator


#####################
# Fitting Functions #
#####################

def gaussian(x, A, mu, sigma):
    return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

def double_gaussian(x, A1, mu1, sigma1, f, delta, r):
    """
    Constrained double Gaussian model:
      Second peak parameters are linked to the first:
      A2 = A1 * f,  mu2 = mu1 + delta, sigma2 = sigma1 * r
    """
    A2 = A1 * f
    mu2 = mu1 + delta
    sigma2 = sigma1 * r
    return gaussian(x, A1, mu1, sigma1) + gaussian(x, A2, mu2, sigma2)

def sum_of_two_gaussians(x, A1, mu1, sigma1, A2, mu2, sigma2):
    """Sum of two independent Gaussian functions."""
    return gaussian(x, A1, mu1, sigma1) + gaussian(x, A2, mu2, sigma2)

def Multi_Gaussfit(x, y, max_gaussians=5):
    """
    Iteratively fits up to max_gaussians to the data and sums them.
    Returns the total fit and a list of parameters for each Gaussian.
    """
    residual = np.copy(y)
    total_fit = np.zeros_like(y)
    params_list = []
    for _ in range(max_gaussians):
        try:
            # Estimate initial guess: amplitude, peak, sigma
            A0 = np.nanmax(residual)
            mu0 = x[np.argmax(residual)]
            sigma0 = (x[-1] - x[0]) / 6
            p0 = [A0, mu0, sigma0]
            popt, _ = curve_fit(gaussian, x, residual, p0=p0, maxfev=10000)
            fit_component = gaussian(x, *popt)
            total_fit += fit_component
            params_list.append(popt)
            residual = residual - fit_component
            # Stop if residual is very low
            if np.nanmax(residual) < 0.05 * A0:
                break
        except Exception as e:
            print(f"Multi-Gauss fit iteration failed: {e}")
            break
    return total_fit, params_list

class PLAnalysisApp:
    def __init__(self, master):
        self.master = master
        self.master.title("PL Analysis Tool")
        self.master.geometry("1920x1080")

        # Create a canvas and a scrollable frame
        self.canvas = Canvas(master)
        self.scrollable_frame = Frame(self.canvas)
        self.v_scrollbar = Scrollbar(master, orient="vertical", command=self.canvas.yview)
        self.h_scrollbar = Scrollbar(master, orient="horizontal", command=self.canvas.xview)
        self.canvas.configure(yscrollcommand=self.v_scrollbar.set, xscrollcommand=self.h_scrollbar.set)
        self.v_scrollbar.pack(side="right", fill="y")
        self.h_scrollbar.pack(side="bottom", fill="x")
        self.canvas.pack(side="left", fill="both", expand=True)
        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.scrollable_frame.bind("<Configure>", self.on_frame_configure)

        self.master.protocol("WM_DELETE_WINDOW", self.on_close)

        # Variables
        self.file_path = None
        self.data = None
        self.wavelength = None
        self.raw_counts = None
        self.time_labels = None
        self.relative_times = None
        self.measurement_date = None

        # QFLS data
        self.qfls_data = []
        self.qfls_times = []

        # For quantification results (to be set by fitting)
        self.selected_fit_times = []      # times corresponding to each spectrum fitted
        self.selected_fit_metrics = []    # metric for each fitted spectrum

        self.plot_spectra_var = IntVar(value=1)
        self.plot_luminescence_flux_density_var = IntVar(value=1)
        self.plot_absolute_gradient_var = IntVar(value=1)
        self.plot_log_spectra_var = IntVar(value=1)
        self.plot_log_luminescence_flux_density_var = IntVar(value=1)
        self.plot_log_absolute_gradient_var = IntVar(value=1)
        self.plot_norm_spectra_var = IntVar(value=1)
        self.plot_norm_luminescence_flux_density_var = IntVar(value=1)
        self.plot_norm_absolute_gradient_var = IntVar(value=1)
        self.show_legend_var = IntVar(value=1)
        self.use_single_fit_var = IntVar(value=1)  # 1 = use single Gaussian fit
        self.use_double_fit_var = IntVar(value=0)  # 1 = use double Gaussian fit
        self.fit_lower_bound_double = StringVar()
        self.fit_upper_bound_double = StringVar()
        self.preview_time_option = StringVar()  # For dropdown preview time; will be set after file load
        self.normalize_fit_var = IntVar(value=0)

        # Quantification fit options variables
        self.fit_type_var = StringVar(value="single")  # "single" or "double"
        self.fit_lower_bound = StringVar()
        self.fit_upper_bound = StringVar()
        # For single Gaussian: initial guesses (as strings)
        self.single_a = StringVar()
        self.single_mu = StringVar()
        self.single_sigma = StringVar()
        # For double Gaussian: initial guesses
        self.double_a1 = StringVar()
        self.double_mu1 = StringVar()
        self.double_sigma1 = StringVar()
        self.double_a2 = StringVar()
        self.double_mu2 = StringVar()
        self.double_sigma2 = StringVar()
        self.use_single_fit_var = tk.IntVar(value=1)
        self.use_double_fit_var = tk.IntVar(value=1)
        self.normalize_fit_var = tk.IntVar(value=1)
        self.use_split_fit_var = IntVar(value=1)

        self.fit_lower_bound = StringVar(value="500")
        self.fit_upper_bound = StringVar(value="700")
        self.fit_lower_bound_double = StringVar(value="500")
        self.fit_upper_bound_double = StringVar(value="900")

        self.preview_time_option = StringVar(value="")  # will be set after loading file

        self.single_a = StringVar(value="")
        self.single_mu = StringVar(value="600")
        self.single_sigma = StringVar(value="5")

        self.double_a1 = StringVar(value="")
        self.double_mu1 = StringVar(value="600")
        self.double_sigma1 = StringVar(value="5")
        self.double_a2 = StringVar(value="")
        self.double_mu2 = StringVar(value="760")
        self.double_sigma2 = StringVar(value="3")

        # Save options
        self.save_raw_var = IntVar(value=1)
        self.save_luminescence_flux_density_var = IntVar(value=1)
        self.save_absolute_gradient_var = IntVar(value=1)
        self.save_log_spectra_var = IntVar(value=1)
        self.save_log_luminescence_flux_density_var = IntVar(value=1)
        self.save_log_absolute_gradient_var = IntVar(value=1)
        self.save_norm_spectra_var = IntVar(value=1)
        self.save_norm_luminescence_flux_density_var = IntVar(value=1)
        self.save_norm_absolute_gradient_var = IntVar(value=1)
        self.save_combined_var = IntVar(value=1)
        self.save_names = {
            "Raw Data": StringVar(value="Raw_Data"),
            "Luminescence Flux Density": StringVar(value="Luminescence_Flux_Density"),
            "Absolute Gradient": StringVar(value="Absolute_Gradient"),
            "Log Spectra": StringVar(value="Log_Spectra"),
            "Log Luminescence Flux Density": StringVar(value="Log_Luminescence_Flux_Density"),
            "Log Absolute Gradient": StringVar(value="Log_Absolute_Gradient"),
            "Normalized Spectra": StringVar(value="Norm_Spectra"),
            "Normalized Luminescence Flux Density": StringVar(value="Norm_Luminescence_Flux_Density"),
            "Normalized Absolute Gradient": StringVar(value="Norm_Absolute_Gradient"),
            "Combined": StringVar(value="Combined_Plots"),
        }

        self.time_checkboxes = []
        self.logo_image = None  # Placeholder for the logo

        self.setup_gui()

    def on_frame_configure(self, event):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def on_close(self):
        response = messagebox.askyesno("Confirm Exit", "Are you sure you want to close? Any unsaved plots will be lost.")
        if response:
            self.master.destroy()

    def setup_gui(self):
        # [Your existing setup_gui code, using self.scrollable_frame instead of self.master]
        Label(self.scrollable_frame, text="PL Analysis Tool", font=("Arial", 16)).grid(row=0, column=0, columnspan=3, pady=10)
        self.load_logo("./hzb_logo.jpg")
        Button(self.scrollable_frame, text="Load File", width=15, command=self.load_file).grid(row=1, column=0, padx=10, pady=5, sticky="w")
        self.file_path_entry = Entry(self.scrollable_frame, width=80, state='disabled')
        self.file_path_entry.grid(row=1, column=1, columnspan=2, sticky="w")
        Label(self.scrollable_frame, text="Timestamp:").grid(row=0, column=0, padx=10, pady=5, sticky="nw")
        self.date_entry = Entry(self.scrollable_frame, width=18, state='disabled')
        self.date_entry.grid(row=0, column=0, padx=10, pady=5, columnspan=1, sticky="w")
        Button(self.scrollable_frame, text="Save Plots", width=15, command=self.save_plots_dialog).grid(row=2, column=0, padx=10, pady=5, sticky="w")
        Button(self.scrollable_frame, text="Plot All", width=15, command=self.plot_in_window).grid(row=2, column=1, pady=5, sticky="w")
        Label(self.scrollable_frame, text="Select Plots:").grid(row=3, column=0, pady=5, sticky="w", padx=10)
        Checkbutton(self.scrollable_frame, text="Spectra", variable=self.plot_spectra_var).grid(row=4, column=0, sticky="w", padx=10)
        Checkbutton(self.scrollable_frame, text="Flux Density", variable=self.plot_luminescence_flux_density_var).grid(row=5, column=0, sticky="w", padx=10)
        Checkbutton(self.scrollable_frame, text="Absolute Gradient", variable=self.plot_absolute_gradient_var).grid(row=6, column=0, sticky="w", padx=10)
        Checkbutton(self.scrollable_frame, text="Log Spectra", variable=self.plot_log_spectra_var).grid(row=4, column=1, sticky="w")
        Checkbutton(self.scrollable_frame, text="Log Flux Density", variable=self.plot_log_luminescence_flux_density_var).grid(row=5, column=1, sticky="w")
        Checkbutton(self.scrollable_frame, text="Log Absolute Gradient", variable=self.plot_log_absolute_gradient_var).grid(row=6, column=1, sticky="w")
        Checkbutton(self.scrollable_frame, text="Normalized Spectra", variable=self.plot_norm_spectra_var).grid(row=7, column=0, sticky="w", padx=10)
        Checkbutton(self.scrollable_frame, text="Normalized Flux Density", variable=self.plot_norm_luminescence_flux_density_var).grid(row=7, column=1, sticky="w")
        Checkbutton(self.scrollable_frame, text="Normalized Absolute Gradient", variable=self.plot_norm_absolute_gradient_var).grid(row=8, column=0, sticky="w", padx=10)
        Button(self.scrollable_frame, text="Wavelength Filter", width=20, command=self.open_filter_wavelength_range_window).grid(row=5, column=2, sticky="w")
        Button(self.scrollable_frame, text="Intensity Filter", width=20, command=self.open_filter_by_intensity_range_window).grid(row=6, column=2, sticky="w")
        Button(self.scrollable_frame, text="Smoothing", width=20, command=self.open_filter_by_moving_average_window).grid(row=7, column=2, sticky="w")
        Button(self.scrollable_frame, text="Show Metadata", width=20, command=self.show_metadata).grid(row=8, column=2, sticky="w")
        Button(self.scrollable_frame, text="Reset Filter", width=20, command=self.reset_filters).grid(row=4, column=2, sticky="w")
        Button(self.scrollable_frame, text="QFLS", width=15, command=self.open_qfls_window).grid(row=4, column=3, padx=10, pady=5, sticky="w")
        Button(self.scrollable_frame, text="Quantify Segregation", width=15, command=self.open_quantification_window).grid(row=5, column=3, padx=10, pady=5, sticky="w")
        Button(self.scrollable_frame, text="Info: Spectra", width=20, anchor="w", command=lambda: self.show_info("Spectra")).grid(row=4, column=4, sticky="w")
        Button(self.scrollable_frame, text="Info: Flux Density", width=20, anchor="w", command=lambda: self.show_info("Flux Density")).grid(row=5, column=4, sticky="w")
        Button(self.scrollable_frame, text="Info: Absolute Gradient", width=20, anchor="w", command=lambda: self.show_info("Absolute Gradient")).grid(row=6, column=4, sticky="w")
        Checkbutton(self.scrollable_frame, text="Show Legend", variable=self.show_legend_var).grid(row=8, column=1, sticky="w")
        self.plot_frame = Frame(self.scrollable_frame, width=1200, height=800, bg="white")
        self.plot_frame.grid(row=9, column=0, columnspan=4, pady=10)
        self.time_check_frame = Frame(self.scrollable_frame)
        self.time_check_frame.grid(row=11, column=0, columnspan=4, pady=5)
        Button(self.scrollable_frame, text="Update Plots", width=15, command=self.update_raw_plot).grid(row=10, column=1)
        self.select_all_var = IntVar(value=1)
        Checkbutton(self.scrollable_frame, text="Toggle All Times", variable=self.select_all_var, command=self.toggle_all_time_checkboxes).grid(row=10, column=2)

    def load_logo(self, file_path):
        try:
            image = Image.open(file_path)
            image = image.resize((200, 100), Image.Resampling.LANCZOS)
            self.logo_image = ImageTk.PhotoImage(image)
            Label(self.scrollable_frame, image=self.logo_image).grid(row=0, column=3, sticky="ne")
        except Exception as e:
            print(f"Error loading logo: {e}")

    ##############################
    # Quantification GUI Methods #
    ##############################

    def open_quantification_window(self):
        """Opens a window to set initial guesses, wavelength bounds, and select fit types and preview time."""
        if not self.relative_times:
            messagebox.showerror("Error", "No relative time data available. Load a file first.")
            return

        if not self.preview_time_option.get():
            self.preview_time_option.set(self.relative_times[0])

        quant_window = Toplevel(self.master)
        quant_window.title("Quantify Halide Segregation")

        # --- Fit Type Selection ---
        Label(quant_window, text="Select Fit Types (tick one or more):", font=("Arial", 12)).grid(row=0, column=0,
                                                                                                  columnspan=2, padx=10,
                                                                                                  pady=5, sticky="w")
        Checkbutton(quant_window, text="Single Gaussian", variable=self.use_single_fit_var).grid(row=0, column=2,
                                                                                                 padx=5, pady=5,
                                                                                                 sticky="w")
        Checkbutton(quant_window, text="Double Gaussian", variable=self.use_double_fit_var).grid(row=0, column=3,
                                                                                                 padx=5, pady=5,
                                                                                                 sticky="w")
        # Remove Sum fit checkbutton and add Split Fit instead:
        Checkbutton(quant_window, text="Split Fit", variable=self.use_split_fit_var).grid(row=0, column=4,
                                                                                          padx=5, pady=5, sticky="w")

        # --- Wavelength Range for Single Gaussian Fitting ---
        Label(quant_window, text="Single Gaussian Wavelength Range (nm):", font=("Arial", 12)).grid(row=1, column=0,
                                                                                                    columnspan=2,
                                                                                                    padx=10, pady=5,
                                                                                                    sticky="w")
        Label(quant_window, text="Min:").grid(row=1, column=2, sticky="w")
        lower_entry = Entry(quant_window, textvariable=self.fit_lower_bound, width=10)
        lower_entry.grid(row=1, column=3, sticky="w", padx=5)
        Label(quant_window, text="Max:").grid(row=1, column=4, sticky="w")
        upper_entry = Entry(quant_window, textvariable=self.fit_upper_bound, width=10)
        upper_entry.grid(row=1, column=5, sticky="w", padx=5)

        # --- Wavelength Range for Double (and Split/Multi) Fitting ---
        Label(quant_window, text="Double Gaussian Wavelength Range (nm):", font=("Arial", 12)).grid(row=2, column=0,
                                                                                                    columnspan=2,
                                                                                                    padx=10, pady=5,
                                                                                                    sticky="w")
        Label(quant_window, text="Min:").grid(row=2, column=2, sticky="w")
        lower_double_entry = Entry(quant_window, textvariable=self.fit_lower_bound_double, width=10)
        lower_double_entry.grid(row=2, column=3, sticky="w", padx=5)
        Label(quant_window, text="Max:").grid(row=2, column=4, sticky="w")
        upper_double_entry = Entry(quant_window, textvariable=self.fit_upper_bound_double, width=10)
        upper_double_entry.grid(row=2, column=5, sticky="w", padx=5)

        # --- Preview Time Selection as Dropdown ---
        Label(quant_window, text="Select Preview Time:", font=("Arial", 12)).grid(row=3, column=0, padx=10, pady=5,
                                                                                  sticky="w")
        OptionMenu(quant_window, self.preview_time_option, *self.relative_times).grid(row=3, column=1, padx=10, pady=5,
                                                                                      sticky="w")

        # --- Initial Guesses for Single Peak ---
        Label(quant_window, text="Initial Guesses: Single Peak", font=("Arial", 12)).grid(row=4, column=0,
                                                                                          columnspan=2, padx=10,
                                                                                          pady=5, sticky="w")
        Label(quant_window, text="Amplitude:").grid(row=4, column=2, sticky="w")
        amp_entry = Entry(quant_window, textvariable=self.single_a, width=10)
        amp_entry.grid(row=4, column=3, sticky="w", padx=5)
        Label(quant_window, text="Peak (nm):").grid(row=4, column=4, sticky="w")
        mu_entry = Entry(quant_window, textvariable=self.single_mu, width=10)
        mu_entry.grid(row=4, column=5, sticky="w", padx=5)
        Label(quant_window, text="Sigma:").grid(row=4, column=6, sticky="w")
        sigma_entry = Entry(quant_window, textvariable=self.single_sigma, width=10)
        sigma_entry.grid(row=4, column=7, sticky="w", padx=5)

        # --- Initial Guesses for Double Peak (Used for Double, Split & Multi Fits) ---
        Label(quant_window, text="Initial Guesses: Double Peak", font=("Arial", 12)).grid(row=5, column=0,
                                                                                          columnspan=2, padx=10,
                                                                                          pady=5, sticky="w")
        Label(quant_window, text="Amp1:").grid(row=5, column=2, sticky="w")
        amp1_entry = Entry(quant_window, textvariable=self.double_a1, width=10)
        amp1_entry.grid(row=5, column=3, sticky="w", padx=5)
        Label(quant_window, text="Peak1 (nm):").grid(row=5, column=4, sticky="w")
        mu1_entry = Entry(quant_window, textvariable=self.double_mu1, width=10)
        mu1_entry.grid(row=5, column=5, sticky="w", padx=5)
        Label(quant_window, text="Sigma1:").grid(row=5, column=6, sticky="w")
        sigma1_entry = Entry(quant_window, textvariable=self.double_sigma1, width=10)
        sigma1_entry.grid(row=5, column=7, sticky="w", padx=5)

        Label(quant_window, text="Amp2:").grid(row=6, column=2, sticky="w")
        amp2_entry = Entry(quant_window, textvariable=self.double_a2, width=10)
        amp2_entry.grid(row=6, column=3, sticky="w", padx=5)
        Label(quant_window, text="Peak2 (nm):").grid(row=6, column=4, sticky="w")
        mu2_entry = Entry(quant_window, textvariable=self.double_mu2, width=10)
        mu2_entry.grid(row=6, column=5, sticky="w", padx=5)
        Label(quant_window, text="Sigma2:").grid(row=6, column=6, sticky="w")
        sigma2_entry = Entry(quant_window, textvariable=self.double_sigma2, width=10)
        sigma2_entry.grid(row=6, column=7, sticky="w", padx=5)

        # --- Normalization Option & Buttons ---
        Checkbutton(quant_window, text="Normalize Single Fit Area", variable=self.normalize_fit_var).grid(row=7,
                                                                                                          column=0,
                                                                                                          columnspan=2,
                                                                                                          padx=10,
                                                                                                          pady=5,
                                                                                                          sticky="w")
        Button(quant_window, text="Apply Fit", command=lambda: self.apply_peak_fit(quant_window)).grid(row=8, column=0,
                                                                                                       columnspan=4,
                                                                                                       pady=10)
        Button(quant_window, text="Preview Fit", command=self.preview_peak_fit).grid(row=8, column=4, columnspan=4,
                                                                                     pady=10)

    def plot_quantification_results_side_by_side(self, times_single, metrics_single, times_double, metrics_double,
                                                 times_sum, metrics_sum):
        """Opens a new window and plots all four quantification results side by side."""
        result_window = Toplevel(self.master)
        result_window.title("Peak Quantification Results (Comparison)")
        fig, axs = plt.subplots(1, 3, figsize=(18, 4))

        axs[0].plot(times_single, metrics_single, marker="o", linestyle="-", color="purple", label="Single Gaussian")
        axs[0].set_xlabel("Time [s]")
        axs[0].set_ylabel("Peak Area [a.u.]")
        axs[0].set_title("Single Gaussian Fit")
        axs[0].legend()
        axs[0].grid(True)

        axs[1].plot(times_double, metrics_double, marker="o", linestyle="-", color="blue", label="Double Gaussian")
        axs[1].set_xlabel("Time [s]")
        axs[1].set_ylabel("Area Ratio")
        axs[1].set_title("Double Gaussian Fit")
        axs[1].legend()
        axs[1].grid(True)

        axs[2].plot(times_sum, metrics_sum, marker="o", linestyle="-", color="green", label="Sum of Two Gaussians")
        axs[2].set_xlabel("Time [s]")
        axs[2].set_ylabel("Area Ratio")
        axs[2].set_title("Split Fit")
        axs[2].legend()
        axs[2].grid(True)

        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=result_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        NavigationToolbar2Tk(canvas, result_window).pack(side="top", fill="x")

    def preview_peak_fit(self):
        """
        Opens a window to preview the fitted curves for the selected time point.
        For each fitting method ticked in the quantification window (Single Gaussian,
        Double Gaussian, Split Fit, and Multi-Gaussian), a separate subplot is shown,
        allowing you to compare which fit works best.
        The full spectrum is shown as blue dots while the fit is computed only on the
        relevant wavelength region.
        (No ratio is displayed in the title.)
        """
        preview_time = self.preview_time_option.get()
        try:
            j = self.relative_times.index(preview_time)
        except ValueError:
            messagebox.showerror("Error", "Selected preview time not found in relative times.")
            return

        x_data = self.wavelength
        y_data = self.raw_counts[:, j] if self.contains_counts_flag and self.raw_counts is not None else self.luminescence_flux_density[:, j] if self.contains_flux_flag and self.luminescence_flux_density is not None else None
        selected_fits = []
        if self.use_single_fit_var.get():
            selected_fits.append("single")
        if self.use_double_fit_var.get():
            selected_fits.append("double")
        if self.use_split_fit_var.get():
            selected_fits.append("split")

        if not selected_fits:
            messagebox.showerror("Error", "No fitting method selected for preview.")
            return

        num_plots = len(selected_fits)
        preview_window = Toplevel(self.master)
        preview_window.title(f"Fit Preview for Time {preview_time}")
        fig, axs = plt.subplots(1, num_plots, figsize=(6 * num_plots, 4))
        if num_plots == 1:
            axs = [axs]
        plot_index = 0

        # --- Single Gaussian Preview ---
        if "single" in selected_fits:
            try:
                try:
                    lower_single = float(self.fit_lower_bound.get())
                    upper_single = float(self.fit_upper_bound.get())
                    indices_single = np.where((x_data >= lower_single) & (x_data <= upper_single))[0]
                except:
                    indices_single = np.arange(len(x_data))
                x_fit = x_data[indices_single]
                y_fit = y_data[indices_single]
                amp_guess = float(self.single_a.get()) if self.single_a.get() else np.nanmax(y_fit)
                mu_guess = float(self.single_mu.get()) if self.single_mu.get() else (x_fit[0] + x_fit[-1]) / 2
                sigma_guess = float(self.single_sigma.get()) if self.single_sigma.get() else (x_fit[-1] - x_fit[0]) / 4
                p0 = [amp_guess, mu_guess, sigma_guess]
                popt, _ = curve_fit(gaussian, x_fit, y_fit, p0=p0)
                fitted_curve = gaussian(x_fit, *popt)
                axs[plot_index].plot(x_data, y_data, "b.", label="Data")
                axs[plot_index].plot(x_fit, fitted_curve, "r-", label="Single Gaussian")
                axs[plot_index].set_title("Single Gaussian Fit")
                axs[plot_index].legend()
            except Exception as e:
                axs[plot_index].text(0.5, 0.5, f"Single fit failed:\n{e}", ha="center", va="center",
                                     transform=axs[plot_index].transAxes)
            plot_index += 1

        # --- Double Gaussian Preview (Constrained) ---
        if "double" in selected_fits:
            try:
                try:
                    lower_double = float(self.fit_lower_bound_double.get())
                    upper_double = float(self.fit_upper_bound_double.get())
                    indices_double = np.where((x_data >= lower_double) & (x_data <= upper_double))[0]
                except:
                    indices_double = np.arange(len(x_data))
                x_fit = x_data[indices_double]
                y_fit = y_data[indices_double]
                p0_double = [
                    float(self.double_a1.get()) if self.double_a1.get() else np.nanmax(y_fit) / 2,
                    float(self.double_mu1.get()) if self.double_mu1.get() else x_fit[np.argmax(y_fit)] - 5,
                    float(self.double_sigma1.get()) if self.double_sigma1.get() else (x_fit[-1] - x_fit[0]) / 8,
                    0.6, 20, 1.2
                ]
                bounds_constrained = ([0, x_fit[0], 0, 0, -np.inf, 0],
                                      [np.inf, x_fit[-1], np.inf, 1, np.inf, np.inf])
                popt, _ = curve_fit(double_gaussian, x_fit, y_fit, p0=p0_double,
                                    bounds=bounds_constrained, method='dogbox')
                fitted_curve = double_gaussian(x_fit, *popt)
                axs[plot_index].plot(x_data, y_data, "b.", label="Data")
                axs[plot_index].plot(x_fit, fitted_curve, "r-", label="Double Gaussian")
                axs[plot_index].set_title("Double Gaussian Fit")
                axs[plot_index].legend()
            except Exception as e:
                axs[plot_index].text(0.5, 0.5, f"Double fit failed:\n{e}", ha="center", va="center",
                                     transform=axs[plot_index].transAxes)
            plot_index += 1

        # --- Split Fit Preview ---
        if "split" in selected_fits:
            try:
                try:
                    lower_double = float(self.fit_lower_bound_double.get())
                    upper_double = float(self.fit_upper_bound_double.get())
                except:
                    lower_double = x_data[0]
                    upper_double = x_data[-1]
                indices_double = np.where((x_data >= lower_double) & (x_data <= upper_double))[0]
                x_fit_full = x_data[indices_double]
                y_fit_full = y_data[indices_double]
                mid = (lower_double + upper_double) / 2
                left_mask = x_fit_full < mid
                right_mask = x_fit_full >= mid
                if np.sum(left_mask) == 0 or np.sum(right_mask) == 0:
                    raise ValueError("Insufficient data on one side of the midpoint.")
                x_left = x_fit_full[left_mask]
                y_left = y_fit_full[left_mask]
                x_right = x_fit_full[right_mask]
                y_right = y_fit_full[right_mask]
                # Use the double peak initial guesses if provided.
                amp_left = float(self.double_a1.get()) if self.double_a1.get() else np.nanmax(y_left)
                mu_left = float(self.double_mu1.get()) if self.double_mu1.get() else (x_left[0] + x_left[-1]) / 2
                sigma_left = float(self.double_sigma1.get()) if self.double_sigma1.get() else (x_left[-1] - x_left[
                    0]) / 4
                p0_left = [amp_left, mu_left, sigma_left]
                popt_left, _ = curve_fit(gaussian, x_left, y_left, p0=p0_left)
                fitted_left = gaussian(x_left, *popt_left)
                amp_right = float(self.double_a2.get()) if self.double_a2.get() else np.nanmax(y_right)
                mu_right = float(self.double_mu2.get()) if self.double_mu2.get() else (x_right[0] + x_right[-1]) / 2
                sigma_right = float(self.double_sigma2.get()) if self.double_sigma2.get() else (x_right[-1] - x_right[
                    0]) / 4
                p0_right = [amp_right, mu_right, sigma_right]
                popt_right, _ = curve_fit(gaussian, x_right, y_right, p0=p0_right)
                fitted_right = gaussian(x_right, *popt_right)
                axs[plot_index].plot(x_data, y_data, "b.", label="Data")
                axs[plot_index].plot(x_left, fitted_left, "m-", label="Left Peak")
                axs[plot_index].plot(x_right, fitted_right, "c-", label="Right Peak")
                axs[plot_index].set_title("Split Fit")
                axs[plot_index].legend()
            except Exception as e:
                axs[plot_index].text(0.5, 0.5, f"Split fit failed:\n{e}", ha="center", va="center",
                                     transform=axs[plot_index].transAxes)
            plot_index += 1

        # --- Multi-Gaussian Preview ---
        if "multi" in selected_fits:
            try:
                total_fit, params_list = Multi_Gaussfit(x_data, y_data)
                if len(params_list) >= 2:
                    A1, mu1, sigma1 = params_list[0]
                    A2, mu2, sigma2 = params_list[1]
                    area1 = A1 * sigma1 * np.sqrt(2 * np.pi)
                    area2 = A2 * sigma2 * np.sqrt(2 * np.pi)
                    ratio = area2 / area1 if area1 != 0 else np.nan
                else:
                    ratio = np.nan
                axs[plot_index].plot(x_data, y_data, "b.", label="Data")
                axs[plot_index].plot(x_data, total_fit, "orange", label="Multi-Gaussian")
                axs[plot_index].set_title("Multi-Gaussian Fit")
                axs[plot_index].legend()
            except Exception as e:
                axs[plot_index].text(0.5, 0.5, f"Multi fit failed:\n{e}", ha="center", va="center",
                                     transform=axs[plot_index].transAxes)
            plot_index += 1

        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=preview_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        NavigationToolbar2Tk(canvas, preview_window).pack(side="top", fill="x")

    def apply_peak_fit(self, quant_window):
        """
        Applies the selected Gaussian fits to the spectral data within the specified wavelength ranges.
        Only uses time points where the checkboxes are ticked.
        Computes metrics for:
          - Single Gaussian: peak area (normalized if chosen)
          - Double Gaussian: ratio of the two peak areas (from the constrained double Gaussian)
          - Split Fit: ratio of areas from separate single Gaussian fits on left and right segments,
                       computed as (peak1 / peak2) where peak2 corresponds to the right segment.
        Then, opens a window with side-by-side plots for comparison.
        """
        try:
            lower_single = float(self.fit_lower_bound.get())
            upper_single = float(self.fit_upper_bound.get())
        except ValueError:
            messagebox.showerror("Error", "Please enter valid single Gaussian wavelength bounds.")
            return

        try:
            lower_double = float(self.fit_lower_bound_double.get())
            upper_double = float(self.fit_upper_bound_double.get())
        except ValueError:
            messagebox.showerror("Error", "Please enter valid double Gaussian wavelength bounds.")
            return

        indices_single = np.where((self.wavelength >= lower_single) & (self.wavelength <= upper_single))[0]
        indices_double = np.where((self.wavelength >= lower_double) & (self.wavelength <= upper_double))[0]

        if self.use_single_fit_var.get() == 1 and len(indices_single) == 0:
            messagebox.showerror("Error", "No data points in the single Gaussian wavelength range.")
            return
        if ((self.use_double_fit_var.get() or (
                hasattr(self, "use_split_fit_var") and self.use_split_fit_var.get()) or self.use_multi_fit_var.get())
                and len(indices_double) == 0):
            messagebox.showerror("Error", "No data points in the double Gaussian wavelength range.")
            return

        single_fit_metrics = []
        single_fit_times = []
        double_fit_metrics = []
        double_fit_times = []
        split_fit_metrics = []
        split_fit_times = []

        selected_time_indices = [i for var, i in self.time_checkboxes if var.get() == 1]
        if not selected_time_indices:
            messagebox.showerror("Error", "No time points selected for fitting.")
            return

        for j in selected_time_indices:
            t = float(self.relative_times[j].replace("s", ""))
            # --- Single Gaussian Fit ---
            if self.use_single_fit_var.get() == 1:
                x_single = self.wavelength[indices_single]
                y_single = self.raw_counts[indices_single, j] if self.contains_counts_flag and self.raw_counts is not None else self.luminescence_flux_density[indices_single, j] if self.contains_flux_flag and self.luminescence_flux_density is not None else None

                try:
                    amp_guess = float(self.single_a.get()) if self.single_a.get() else np.nanmax(y_single)
                    mu_guess = float(self.single_mu.get()) if self.single_mu.get() else (x_single[0] + x_single[-1]) / 2
                    sigma_guess = float(self.single_sigma.get()) if self.single_sigma.get() else (x_single[-1] -
                                                                                                  x_single[0]) / 4
                    p0 = [amp_guess, mu_guess, sigma_guess]
                    popt, _ = curve_fit(gaussian, x_single, y_single, p0=p0)
                    area = popt[0] * popt[2] * np.sqrt(2 * np.pi)
                    single_fit_metrics.append(area)
                    single_fit_times.append(t)
                except Exception as e:
                    print(f"Single fit failed for time index {j}: {e}")

            # --- Double Gaussian Fit (Constrained) ---
            if self.use_double_fit_var.get() == 1:
                x_double = self.wavelength[indices_double]
                y_double = self.raw_counts[
                    indices_double, j] if self.contains_counts_flag and self.raw_counts is not None else \
                    self.luminescence_flux_density[
                        indices_double, j] if self.contains_flux_flag and self.luminescence_flux_density is not None else None
                try:
                    p0_double = [
                        float(self.double_a1.get()) if self.double_a1.get() else np.nanmax(y_double) / 2,
                        float(self.double_mu1.get()) if self.double_mu1.get() else x_double[np.argmax(y_double)] - 5,
                        float(self.double_sigma1.get()) if self.double_sigma1.get() else (x_double[-1] - x_double[
                            0]) / 8,
                        0.6, 20, 1.2
                    ]
                    bounds_constrained = ([0, x_double[0], 0, 0, -np.inf, 0],
                                          [np.inf, x_double[-1], np.inf, 1, np.inf, np.inf])
                    popt, _ = curve_fit(double_gaussian, x_double, y_double, p0=p0_double,
                                        bounds=bounds_constrained, method='dogbox')
                    area1 = popt[0] * popt[2] * np.sqrt(2 * np.pi)
                    area2 = (popt[0] * popt[3]) * (popt[2] * popt[5]) * np.sqrt(2 * np.pi)
                    ratio = area1 / area2 if area2 != 0 else np.nan
                    double_fit_metrics.append(ratio)
                    double_fit_times.append(t)
                except Exception as e:
                    print(f"Double fit failed for time index {j}: {e}")

            # --- Split Fit (Separate Single Fits on Left & Right) ---
            if hasattr(self, "use_split_fit_var") and self.use_split_fit_var.get():
                try:
                    x_full = self.wavelength[indices_double]
                    y_full = self.raw_counts[indices_double, j] if self.contains_counts_flag and self.raw_counts is not None else \
                    self.luminescence_flux_density[indices_double, j] if self.contains_flux_flag and self.luminescence_flux_density is not None else None

                    mid = (lower_double + upper_double) / 2
                    left_mask = x_full < mid
                    right_mask = x_full >= mid
                    if np.sum(left_mask) == 0 or np.sum(right_mask) == 0:
                        raise ValueError("Insufficient data on one side of the midpoint.")
                    x_left = x_full[left_mask]
                    y_left = y_full[left_mask]
                    x_right = x_full[right_mask]
                    y_right = y_full[right_mask]
                    # Use double peak initial guesses if provided:
                    amp_left = float(self.double_a1.get()) if self.double_a1.get() else np.nanmax(y_left)
                    mu_left = float(self.double_mu1.get()) if self.double_mu1.get() else (x_left[0] + x_left[-1]) / 2
                    sigma_left = float(self.double_sigma1.get()) if self.double_sigma1.get() else (x_left[-1] - x_left[
                        0]) / 4
                    p0_left = [amp_left, mu_left, sigma_left]
                    popt_left, _ = curve_fit(gaussian, x_left, y_left, p0=p0_left)
                    amp_right = float(self.double_a2.get()) if self.double_a2.get() else np.nanmax(y_right)
                    mu_right = float(self.double_mu2.get()) if self.double_mu2.get() else (x_right[0] + x_right[-1]) / 2
                    sigma_right = float(self.double_sigma2.get()) if self.double_sigma2.get() else (x_right[-1] -
                                                                                                    x_right[0]) / 4
                    p0_right = [amp_right, mu_right, sigma_right]
                    popt_right, _ = curve_fit(gaussian, x_right, y_right, p0=p0_right)
                    area_left = popt_left[0] * popt_left[2] * np.sqrt(2 * np.pi)
                    area_right = popt_right[0] * popt_right[2] * np.sqrt(2 * np.pi)
                    # Always compute peak ratio as (peak2 / peak1)
                    ratio = area_left / area_right if area_right != 0 else np.nan
                    split_fit_metrics.append(ratio)
                    split_fit_times.append(t)
                except Exception as e:
                    print(f"Split fit failed for time index {j}: {e}")

        if self.normalize_fit_var.get() == 1:
            single_fit_metrics = single_fit_metrics / np.max(single_fit_metrics)

        quant_window.destroy()

        self.plot_quantification_results_side_by_side(single_fit_times, single_fit_metrics,
                                                      double_fit_times, double_fit_metrics,
                                                      split_fit_times, split_fit_metrics)

    def show_option_info(self, filter_type):
        info_messages = {
            "Filter by Wavelength Range": (
                "Filter by Wavelength Range:\n\n"
                "This filter allows you to specify a wavelength range (min and max in nm).\n"
                "All data points outside this range will be removed from the dataset.\n\n"
                "Example:\n- Min Wavelength: 400\n- Max Wavelength: 800"
            ),
            "Filter by Intensity Threshold": (
                "Filter by Intensity Threshold:\n\n"
                "This filter removes data points where the intensity is above or below the specified thresholds.\n"
                "You can use this to exclude low- or high-intensity data from the analysis.\n\n"
                "If one field is left blank, the filter will not apply to that bound.\n\n"
                "Example:\n Min=100, Max=1000:\n\n Retains intensities between 100 and 1000.\n"
            ),
            "Filter by Moving Average": (
                "Filter by Moving Average:\n\n"
                "This filter smooths the data using a moving average with a specified window size.\n"
                "A larger window size results in stronger smoothing.\n\n"
                "Example:\n- Window Size: 5"
            ),
        }
        messagebox.showinfo("Filter Info", info_messages.get(filter_type, "No information available."))

    def apply_filter_by_wavelength_range(self, min_wavelength, max_wavelength, filter_window):
        """Filters the data based on the provided wavelength range."""
        try:
            # Use full range if inputs are empty
            min_wavelength = float(min_wavelength) if min_wavelength.strip() else float('-inf')
            max_wavelength = float(max_wavelength) if max_wavelength.strip() else float('inf')

            print(f"Applying wavelength filter with Min: {min_wavelength}, Max: {max_wavelength}")

            # Mask for valid wavelengths
            wavelength_mask = (self.wavelength >= min_wavelength) & (self.wavelength <= max_wavelength)

            # Apply mask to wavelength data
            self.wavelength = self.wavelength[wavelength_mask]

            # Apply mask to raw counts if available
            if self.contains_counts_flag and self.raw_counts is not None:
                if self.raw_counts.ndim == 1:  # Single measurement
                    self.raw_counts = self.raw_counts[wavelength_mask]
                else:  # Continuous measurement
                    self.raw_counts = self.raw_counts[wavelength_mask, :]

            # Apply mask to luminescence flux density if available
            if self.contains_flux_flag and self.luminescence_flux_density is not None:
                if self.luminescence_flux_density.ndim == 1:  # Single measurement
                    self.luminescence_flux_density = self.luminescence_flux_density[wavelength_mask]
                else:  # Continuous measurement
                    self.luminescence_flux_density = self.luminescence_flux_density[wavelength_mask, :]

            # Update plots
            self.update_raw_plot()

            # Close filter window
            filter_window.destroy()

        except ValueError as e:
            print(f"Invalid wavelength input: {e}")
            messagebox.showerror("Error", "Please enter valid numerical values for wavelength range.")
        except Exception as e:
            print(f"An error occurred while applying the wavelength filter: {e}")
            messagebox.showerror("Error", f"An unexpected error occurred: {e}")

    def plot_qfls(self, ax):
        """Plots QFLS extracted from iVoc values with confidence = 1 based on selected time checkboxes."""
        ax.clear()  # Clear previous plot

        if not self.qfls_times or not self.qfls_data:
            ax.text(0.5, 0.5, "No data available", ha='center', va='center', fontsize=12)
            return

        # Get the selected times from the checkboxes
        selected_times = [self.relative_times[i] for var, i in self.time_checkboxes if var.get() == 1]

        # Debugging: Print selected times and qfls_times
        print(f"Selected times (checkboxes): {selected_times}")
        print(f"QFLS times: {self.qfls_times}")

        self.selected_qfls_data = []
        self.selected_qfls_times = []

        # Match selected times with qfls_times based on values (not index)
        for i, qfls_time in enumerate(self.qfls_times):
            if qfls_time in selected_times:
                self.selected_qfls_data.append(self.qfls_data[i])
                self.selected_qfls_times.append(qfls_time)

        # If no data is found for the selected times, show a message
        if not self.selected_qfls_data:
            ax.text(0.5, 0.5, "No matching QFLS data for selected times", ha='center', va='center', fontsize=12)
            print(f"No matching QFLS data for selected times: {selected_times}")
            return

        # Debugging: Print the selected data
        print(f"Selected QFLS data: {self.selected_qfls_data}")
        print(f"Selected QFLS times: {self.selected_qfls_times}")

        # Plot only the selected QFLS data
        ax.plot(self.selected_qfls_times, self.selected_qfls_data, marker="o", linestyle="-", color="blue", label="QFLS")

        # Use MaxNLocator to automatically adjust the number of ticks on x-axis
        ax.xaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))

        ax.set_xlabel("Time (s)")
        ax.set_ylabel("QFLS (eV)")
        ax.set_title("Quasi Fermi Level Splitting Over Time")
        ax.legend()
        ax.grid(True)

    def open_qfls_window(self):
        """Creates the window for displaying QFLS plot."""
        if not self.qfls_data or not self.qfls_times:
            messagebox.showwarning("No Data", "No QFLS data available to plot.")
            return

        # Create a new window
        qfls_window = tk.Toplevel(self.master)
        qfls_window.title("QFLS Plot")

        # Create a figure and canvas
        fig, ax = plt.subplots(figsize=(6, 4))
        self.plot_qfls(ax)

        # Embed the plot inside the Tkinter window
        canvas = FigureCanvasTkAgg(fig, master=qfls_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Button Frame
        button_frame = tk.Frame(qfls_window)
        button_frame.pack(fill=tk.X, padx=10, pady=5)

        def save_plot():
            """Saves the plot to a file."""
            file_path = filedialog.asksaveasfilename(defaultextension=".png",
                                                     filetypes=[("PNG Image", "*.png"), ("JPG Image", "*.jpg")])
            if file_path:
                fig.savefig(file_path, dpi=300)
                print(f"Plot saved as {file_path}")

        def save_data():
            """Saves the QFLS data to a text file."""
            if not hasattr(self, "selected_qfls_times") or not hasattr(self, "selected_qfls_data"):
                messagebox.showerror("Error", "No QFLS data selected to save!")
                return

            file_path = filedialog.asksaveasfilename(defaultextension=".txt",
                                                     filetypes=[("Text Files", "*.txt")])
            if not file_path:
                return  # User canceled the file save

            try:
                with open(file_path, "w") as file:
                    file.write("Time (s)\tQFLS (eV)\n")
                    for time, qfls in zip(self.selected_qfls_times, self.selected_qfls_data):
                        file.write(f"{time}\t{qfls}\n")
                messagebox.showinfo("Success", "QFLS data saved successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save data: {e}")

        # Buttons
        save_plot_btn = Button(button_frame, text="Save Plot", command=save_plot)
        save_plot_btn.pack(side=tk.LEFT, padx=5, pady=5)

        save_data_btn = Button(button_frame, text="Save Data", command=save_data)
        save_data_btn.pack(side=tk.RIGHT, padx=5, pady=5)

    def close_qfls_window(self):
        """Handle window closure to clean up and prevent multiple windows."""
        if hasattr(self, 'qfls_window'):
            self.qfls_window.destroy()
            del self.qfls_window

    def open_filter_wavelength_range_window(self):
        """Open a filter window to specify the wavelength range."""
        filter_window = Toplevel(self.master)
        filter_window.title("Filter by Wavelength Range")

        # Add labels and entry fields
        Label(filter_window, text="Min Wavelength (nm):", font=("Arial", 12)).pack(pady=5, anchor="w", padx=10)
        min_wavelength_entry = Entry(filter_window, width=15)
        min_wavelength_entry.pack(pady=5, anchor="w", padx=10)

        Label(filter_window, text="Max Wavelength (nm):", font=("Arial", 12)).pack(pady=5, anchor="w", padx=10)
        max_wavelength_entry = Entry(filter_window, width=15)
        max_wavelength_entry.pack(pady=5, anchor="w", padx=10)

        # Add Apply button
        Button(filter_window, text="Apply",
               command=lambda: self.apply_filter_by_wavelength_range(
                   min_wavelength_entry.get(),
                   max_wavelength_entry.get(),
                   filter_window  # Pass the window reference to close it
               )).pack(pady=10)

        Button(filter_window, text="Info",
               command=lambda: self.show_option_info("Filter by Wavelength Range")).pack(anchor="w")

    def apply_filter_by_intensity_range(self, min_intensity, max_intensity, filter_window):
        try:
            min_intensity = float(min_intensity) if min_intensity.strip() else float('-inf')
            max_intensity = float(max_intensity) if max_intensity.strip() else float('inf')
            print(f"Applying intensity filter with Min: {min_intensity}, Max: {max_intensity}")

            mask = None

            # Check and filter raw counts
            if self.contains_counts_flag and self.raw_counts is not None:
                if len(self.raw_counts.shape) == 2:  # Continuous measurement
                    mask_raw = (self.raw_counts >= min_intensity) & (self.raw_counts <= max_intensity)
                    mask = mask_raw.all(axis=1) if mask is None else mask & mask_raw.all(axis=1)
                    self.raw_counts = self.raw_counts[mask]
                else:  # Single measurement
                    mask_raw = (self.raw_counts >= min_intensity) & (self.raw_counts <= max_intensity)
                    mask = mask_raw if mask is None else mask & mask_raw
                    self.raw_counts = self.raw_counts[mask]

            # Check and filter luminescence flux density
            if self.contains_flux_flag and self.luminescence_flux_density is not None:
                if len(self.luminescence_flux_density.shape) == 2:  # Continuous measurement
                    mask_flux = (self.luminescence_flux_density >= min_intensity) & (
                                self.luminescence_flux_density <= max_intensity)
                    mask = mask_flux.all(axis=1) if mask is None else mask & mask_flux.all(axis=1)
                    self.luminescence_flux_density = self.luminescence_flux_density[mask]
                else:  # Single measurement
                    mask_flux = (self.luminescence_flux_density >= min_intensity) & (
                                self.luminescence_flux_density <= max_intensity)
                    mask = mask_flux if mask is None else mask & mask_flux
                    self.luminescence_flux_density = self.luminescence_flux_density[mask]

            # Apply the mask to the wavelength array
            if mask is not None:
                self.wavelength = self.wavelength[mask]

            # Close the filter window and update the plots
            filter_window.destroy()
            self.update_raw_plot()

        except Exception as e:
            print(f"An error occurred while applying the intensity filter: {e}")

    def open_filter_by_intensity_range_window(self):
        """Opens a window to set intensity filter range."""
        filter_window = Toplevel(self.master)
        filter_window.title("Filter by Intensity Range")

        Label(filter_window, text="Set Intensity Range:", font=("Arial", 12)).pack(pady=10)

        # Input fields for min and max intensity
        Label(filter_window, text="Min Intensity:").pack(anchor="w", padx=10)
        min_intensity_entry = Entry(filter_window, width=15)
        min_intensity_entry.pack(anchor="w", padx=10)

        Label(filter_window, text="Max Intensity:").pack(anchor="w", padx=10)
        max_intensity_entry = Entry(filter_window, width=15)
        max_intensity_entry.pack(anchor="w", padx=10)

        # Apply button
        Button(filter_window, text="Apply",
               command=lambda: self.apply_filter_by_intensity_range(
                   min_intensity_entry.get(),
                   max_intensity_entry.get(),
                   filter_window
               )).pack(pady=10)

        Button(filter_window, text="Info",
               command=lambda: self.show_option_info("Filter by Intensity Threshold")).pack(anchor="w")

    def apply_filter_by_moving_average(self, window_size, filter_window):
        if self.raw_counts is not None:
            self.raw_counts = np.apply_along_axis(
                lambda m: np.convolve(m, np.ones(window_size) / window_size, mode='valid'),
                axis=0, arr=self.raw_counts
            )
            self.wavelength = self.wavelength[:len(self.raw_counts)]
        if self.luminescence_flux_density is not None:
            self.luminescence_flux_density = np.apply_along_axis(
                lambda m: np.convolve(m, np.ones(window_size) / window_size, mode='valid'),
                axis=0, arr=self.luminescence_flux_density
            )
            self.wavelength = self.wavelength[:len(self.luminescence_flux_density)]

        filter_window.destroy()  # Close the filter window
        self.update_raw_plot()

    def open_filter_by_moving_average_window(self):
        filter_window = Toplevel(self.master)
        filter_window.title("Filter by Moving Average")

        Label(filter_window, text="Enter Window Size:", font=("Arial", 12)).pack(pady=10)
        window_size_entry = Entry(filter_window, width=10)
        window_size_entry.pack(padx=10)

        Button(filter_window, text="Apply",
               command=lambda: self.apply_filter_by_moving_average(
                   int(window_size_entry.get()), filter_window)
               ).pack(pady=10)

        Button(filter_window, text="Info",
               command=lambda: self.show_option_info("Filter by Moving Average")).pack(anchor="w")

    def reset_filters(self):
        if not self.file_path:
            messagebox.showinfo("No File Loaded", "Please load a file to reset filters.")
            return

        # Ensure self.data and self.headers are not None
        if self.data is None or self.headers is None:
            messagebox.showerror("Error", "Data or headers not loaded correctly.")
            return

        try:
            # Reinitialize attributes
            self.raw_counts = None
            self.luminescence_flux_density = None
            self.contains_counts_flag = False
            self.contains_flux_flag = False
            self.is_single_measurement_flag = None

            # Reload wavelength data
            self.wavelength = self.data.iloc[:, 0].values

            # Determine the available columns
            columns = self.headers.columns[1:].tolist()
            # Convert columns to lowercase for case-insensitive comparison
            columns_lower = [col.lower() for col in columns]
            print("Columns:", columns)  # Debug print

            # Check for flux density column
            if any('flux' in col and 'density' in col for col in columns_lower):
                matching_col_flux = next(col for col in columns if col.lower() in columns_lower)
                col_index_flux = columns.index(matching_col_flux) + 1  # +1 because wavelength is first column
                if self.is_single_measurement_flag:
                    self.luminescence_flux_density = self.data.iloc[:, col_index_flux:].values  # Store as 1D array
                else:
                    self.luminescence_flux_density = self.data.iloc[:, col_index_flux:].values  # Store as 2D array
                self.contains_flux_flag = True
                print(f"Luminescence flux density loaded from column: {matching_col_flux}")

            # Check for raw counts column
            if any('raw' in col and 'counts' in col for col in columns_lower):
                matching_col_raw = next(col for col in columns if col.lower() in columns_lower)
                if self.contains_flux_flag:
                    col_index_raw = columns.index(
                        matching_col_raw) + 2  # +2 because wavelength is first column and flux second
                else:
                    col_index_raw = columns.index(matching_col_raw) + 1  # +1 because wavelength is first column
                if self.is_single_measurement_flag:
                    self.raw_counts = self.data.iloc[:, col_index_raw:].values  # Store as 1D array
                else:
                    self.raw_counts = self.data.iloc[:, col_index_raw:].values  # Store as 2D array
                self.contains_counts_flag = True
                print(f"Raw counts loaded from column: {matching_col_raw}")

            # Read the raw time labels from the file
            with open(self.file_path, 'r') as file:
                lines = file.readlines()
                # For continuous measurements
                if "\t" in lines[0]:  # First line contains tab-separated timestamps
                    date = lines[0].strip().split("\t")[0]  # Extract date
                    time = lines[0].strip().split("\t")[1]  # Extract time
                    date_time_string = f"{date} {time}"
                    date_time = datetime.strptime(date_time_string, "%d.%m.%Y %H:%M:%S")
                else:  # For single measurements
                    date_time_string = lines[0].strip()  # Use the first line as a single date-time string
                    date_time = datetime.strptime(date_time_string, "%d.%m.%Y %H:%M:%S")

                # Update the GUI date entry
                self.date_entry.config(state='normal')
                self.date_entry.delete(0, 'end')
                self.date_entry.insert(0, date_time.strftime("%d.%m.%Y %H:%M:%S"))
                self.date_entry.config(state='disabled')

            if len(lines) > 0:
                raw_time_labels = lines[0].strip().split("\t")[1:]
                if len(raw_time_labels) > 1:
                    self.relative_times = self.calculate_relative_times(raw_time_labels)
                    self.is_single_measurement_flag = False
                else:
                    self.relative_times = []
                    self.is_single_measurement_flag = True

            self.setup_time_checkboxes()

            # Determine if the file contains single or continuous measurements
            self.is_single_measurement = (self.raw_counts is not None and self.raw_counts.ndim == 1) or (
                    self.luminescence_flux_density is not None and self.luminescence_flux_density.ndim == 1)

            # Add debug prints to check data shapes
            print(f"Raw counts shape: {self.raw_counts.shape if self.raw_counts is not None else 'None'}")
            print(
                f"Luminescence flux density shape: {self.luminescence_flux_density.shape if self.luminescence_flux_density is not None else 'None'}")
            print(f"contains_counts_flag: {self.contains_counts_flag}")
            print(f"contains_flux_flag: {self.contains_flux_flag}")

            self.update_raw_plot()

        except pd.errors.EmptyDataError:
            messagebox.showerror("Error", "Error loading file: No columns to parse from file")
            self.data = None
        except ValueError as e:
            messagebox.showerror("Error", "Error loading file: {}".format(e))
            self.data = None
        except IndexError:
            messagebox.showerror("Error", "Error parsing times: list index out of range")
            self.data = None
        except Exception as e:
            messagebox.showerror("Error", "An unexpected error occurred: {}".format(e))
            self.data = None

    def show_metadata(self):
        if not self.file_path:
            messagebox.showinfo("No File Loaded", "Please load a file to view its metadata.")
            return

        # Extract metadata
        metadata = []
        with open(self.file_path, 'r') as file:
            for line in file:
                if "Wavelength (nm)" in line:
                    break
                metadata.append(line.strip())

        # Create metadata window
        metadata_window = Toplevel(self.master)
        metadata_window.title("Metadata")
        metadata_window.geometry("600x400")

        # Add a scrollable text widget
        text_widget = Text(metadata_window, wrap="word")
        scrollbar = Scrollbar(metadata_window, command=text_widget.yview)
        text_widget.configure(yscrollcommand=scrollbar.set)

        # Populate the text widget with metadata
        for line in metadata:
            text_widget.insert("end", line + "\n")

        text_widget.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        text_widget.config(state="disabled")

    def is_single_measurement(self):
        self.is_single_measurement_flag = False

    def contains_counts(self):
        self.contains_counts_flag = False

    def contains_flux(self):
        self.contains_flux_flag = False

    def load_file(self):
        self.file_path = filedialog.askopenfilename(title="Select Data File", filetypes=[("Text Files", "*.txt")])
        if not self.file_path:
            return

        self.file_path_entry.config(state='normal')
        self.file_path_entry.delete(0, 'end')
        self.file_path_entry.insert(0, self.file_path)
        self.file_path_entry.config(state='disabled')

        # Initialize attributes
        self.raw_counts = None
        self.luminescence_flux_density = None
        self.contains_counts_flag = False
        self.contains_flux_flag = False
        self.is_single_measurement_flag = None
        self.qfls_data = []
        self.qfls_times = []

        # Scan the file to find the header line
        header_line = "Wavelength (nm)"
        skip_rows = 0
        encoding_used = 'utf-8'  # Default encoding

        # Read file as raw text first to extract metadata
        ivoc_values, ivoc_confidence = None, None

        try:
            with open(self.file_path, 'r', encoding='utf-8') as file:
                lines = file.readlines()
        except UnicodeDecodeError:
            with open(self.file_path, 'r', encoding='ISO-8859-1') as file:
                lines = file.readlines()
            encoding_used = 'ISO-8859-1'

        # Extract the first timestamp from metadata
        if len(lines) > 0:
            raw_time_labels = lines[0].strip().split("\t")[1:]  # Extract all timestamps from the first line

            if len(raw_time_labels) > 0:
                first_time_str = raw_time_labels[0]  # The first timestamp in the metadata

                # Determine the correct time format
                if " " in first_time_str:
                    time_format = "%d.%m.%Y %H:%M:%S"
                else:
                    time_format = "%H:%M:%S"

                try:
                    first_time = datetime.strptime(first_time_str, time_format)
                    # Update the timestamp entry box
                    self.date_entry.config(state='normal')
                    self.date_entry.delete(0, 'end')
                    self.date_entry.insert(0, first_time.strftime("%d.%m.%Y %H:%M:%S"))
                    self.date_entry.config(state='disabled')
                except ValueError:
                    print(f"Error: Unable to parse first timestamp '{first_time_str}'")

        # Extract iVoc and iVoc confidence
        for line in lines:
            if line.startswith("iVoc (V)"):
                ivoc_values = line.strip().split("\t")[1:]  # Skip the label
            elif line.startswith("iVoc confidence"):
                ivoc_confidence = line.strip().split("\t")[1:]

        # Convert values and filter by confidence
        if ivoc_values and ivoc_confidence:
            ivoc_values = [float(val) if val.replace('.', '', 1).isdigit() else np.nan for val in ivoc_values]
            ivoc_confidence = [int(float(val)) if val.replace('.', '', 1).isdigit() else 0 for val in ivoc_confidence]

            # Compute QFLS only for valid iVoc values where confidence is 1
            self.qfls_data = [ivoc_values[i] for i in range(len(ivoc_values)) if ivoc_confidence[i] == 1]

            # Use relative times for valid data
            if raw_time_labels:
                self.relative_times = self.calculate_relative_times(raw_time_labels)
                if self.relative_times is not None:
                    self.qfls_times = [self.relative_times[i] for i in range(len(ivoc_values)) if
                                       ivoc_confidence[i] == 1]
                else:
                    print("Error: Unable to calculate relative times.")
            else:
                print("Error: No raw time labels found.")

            print(f"Extracted {len(self.qfls_data)} QFLS values with confidence 1")
        else:
            print("Warning: No valid iVoc (QFLS) data found with confidence 1")

        # Now read the main spectral data
        def read_file_with_encoding(encoding):
            nonlocal skip_rows, encoding_used
            encoding_used = encoding
            with open(self.file_path, 'r', encoding=encoding) as file:
                for line in file:
                    if header_line in line:
                        break
                    skip_rows += 1

            self.data = pd.read_csv(self.file_path, sep="\t", skiprows=skip_rows + 1, encoding=encoding)
            self.headers = pd.read_csv(self.file_path, sep="\t", skiprows=skip_rows, encoding=encoding)

        try:
            read_file_with_encoding(encoding_used)
        except Exception as e:
            messagebox.showerror("Error", f"Error reading spectral data: {e}")
            return

        # Fix column headers if needed
        if "Unnamed" in self.headers.columns[0]:
            self.headers.columns = self.data.columns  # Try fixing headers

        self.wavelength = self.data.iloc[:, 0].values
        columns = self.headers.columns[1:].tolist()
        columns_lower = [col.lower() for col in columns]

        # Extract other data (Flux density, raw counts, etc.)
        if any('flux' in col and 'density' in col for col in columns_lower):
            matching_col_flux = next(col for col in columns if col.lower() in columns_lower)
            col_index_flux = columns.index(matching_col_flux) + 1
            self.luminescence_flux_density = self.data.iloc[:, col_index_flux:].values
            self.contains_flux_flag = True
            print(f"Luminescence flux density loaded from column: {matching_col_flux}")

        if any('raw' in col and 'counts' in col for col in columns_lower):
            matching_col_raw = next(col for col in columns if col.lower() in columns_lower)
            col_index_raw = columns.index(matching_col_raw) + (2 if self.contains_flux_flag else 1)
            self.raw_counts = self.data.iloc[:, col_index_raw:].values
            self.contains_counts_flag = True
            print(f"Raw counts loaded from column: {matching_col_raw}")

        # Read the raw time labels from the file
        with open(self.file_path, 'r') as file:
            lines = file.readlines()
            # For continuous measurements
            if "\t" in lines[0]:  # First line contains tab-separated timestamps
                date = lines[0].strip().split("\t")[0]  # Extract date
                time = lines[0].strip().split("\t")[1]  # Extract time
                date_time_string = f"{date} {time}"
                date_time = datetime.strptime(date_time_string, "%d.%m.%Y %H:%M:%S")
            else:  # For single measurements
                date_time_string = lines[0].strip()  # Use the first line as a single date-time string
                date_time = datetime.strptime(date_time_string, "%d.%m.%Y %H:%M:%S")

            # Update the GUI date entry
            self.date_entry.config(state='normal')
            self.date_entry.delete(0, 'end')
            self.date_entry.insert(0, date_time.strftime("%d.%m.%Y %H:%M:%S"))
            self.date_entry.config(state='disabled')

        # Check if relative times can be calculated
        if len(lines) > 0:
            raw_time_labels = lines[0].strip().split("\t")[1:]  # Extract all timestamps from the first line
            if len(raw_time_labels) > 1:
                self.relative_times = self.calculate_relative_times(raw_time_labels)
                self.is_single_measurement_flag = False
            else:
                self.relative_times = []
                self.is_single_measurement_flag = True

        self.setup_time_checkboxes()

        self.is_single_measurement = (
                (self.raw_counts is not None and self.raw_counts.ndim == 1) or
                (self.luminescence_flux_density is not None and self.luminescence_flux_density.ndim == 1)
        )

        # Debug prints
        print(f"Raw counts shape: {self.raw_counts.shape if self.raw_counts is not None else 'None'}")
        print(
            f"Luminescence flux density shape: {self.luminescence_flux_density.shape if self.luminescence_flux_density is not None else 'None'}")
        print(f"contains_counts_flag: {self.contains_counts_flag}")
        print(f"contains_flux_flag: {self.contains_flux_flag}")

    def calculate_relative_times(self, raw_times):
        """
        Parses a list of timestamp strings and computes relative times (in seconds)
        based on consecutive differences. The first time is set to 0. For each subsequent
        timestamp, if the gap from the previous timestamp is negative or more than three
        times the median gap (computed from all positive consecutive differences),
        the gap is replaced with the median gap.
        Returns a list of strings (e.g. "0s", "2s", "4s", …).
        """
        if len(raw_times) <= 1:
            return None
        try:
            # Determine the time format.
            if " " in raw_times[0]:
                time_format = "%d.%m.%Y %H:%M:%S"
            else:
                time_format = "%H:%M:%S"

            # Parse all timestamps.
            times = [datetime.strptime(ts, time_format) for ts in raw_times]

            # Compute differences (in seconds) between consecutive timestamps.
            diffs = []
            for i in range(1, len(times)):
                d = (times[i] - times[i - 1]).total_seconds()
                diffs.append(d)

            # Compute the typical (median) gap using only positive differences.
            valid_diffs = [d for d in diffs if d > 0]
            median_gap = np.median(valid_diffs) if valid_diffs else 2

            # Build the cumulative relative times, starting at 0.
            cum_times = [0]
            for d in diffs:
                # If the gap is negative or unusually large, replace it with median_gap.
                if d <= 0 or d > 3 * median_gap:
                    cum_times.append(cum_times[-1] + median_gap)
                else:
                    cum_times.append(cum_times[-1] + d)

            # Return the cumulative times as strings with "s" appended.
            return [f"{int(t)}s" for t in cum_times]
        except Exception as e:
            print("Error parsing times:", e)
            return None

    def setup_time_checkboxes(self):
        # Clear existing widgets
        for widget in self.time_check_frame.winfo_children():
            widget.destroy()
        self.time_checkboxes = []

        # If it's a single measurement or no relative times, don't create checkboxes
        if self.is_single_measurement_flag or not self.relative_times:
            print("Single measurement or no relative times - skipping checkbox setup")
            # Set a default "selected" state for the single measurement
            self.time_checkboxes = [(IntVar(value=1), 0)]
            return

        # Define the maximum number of checkboxes per row
        max_per_row = 15

        # Dynamically create checkboxes in a grid layout
        for i, time_label in enumerate(self.relative_times):
            var = IntVar(value=1)  # Default value is checked
            cb = Checkbutton(self.time_check_frame, text=time_label, variable=var)

            # Calculate grid row and column
            row = i // max_per_row
            col = i % max_per_row

            # Place the checkbox in the grid
            cb.grid(row=row, column=col, padx=5, pady=2, sticky="w")

            # Store the variable and index
            self.time_checkboxes.append((var, i))

    def toggle_all_time_checkboxes(self):
        """
        Ticks/unticks all time checkboxes based on the 'Select All' state.
        """
        for var, _ in self.time_checkboxes:
            var.set(self.select_all_var.get())  # Set all checkboxes to match the Select All state

    def save_plots_dialog(self):
        save_window = Toplevel(self.master)
        save_window.title("Save Options")
        Label(save_window, text="Select Plots to Save:").grid(row=0, column=0, columnspan=2, pady=10)
        # Checkboxes for save options
        Checkbutton(save_window, text="Raw Counts", variable=self.save_raw_var).grid(row=1, column=0, sticky="w")
        Checkbutton(save_window, text="Flux Density", variable=self.save_luminescence_flux_density_var).grid(row=2, column=0, sticky="w")
        Checkbutton(save_window, text="Absolute Gradient", variable=self.save_absolute_gradient_var).grid(row=3, column=0, sticky="w")
        Checkbutton(save_window, text="Log Raw Counts", variable=self.save_log_spectra_var).grid(row=4, column=0, sticky="w")
        Checkbutton(save_window, text="Log Flux Density", variable=self.save_log_luminescence_flux_density_var).grid(row=5,column=0,sticky="w")
        Checkbutton(save_window, text="Log Absolute Gradient", variable=self.save_log_absolute_gradient_var).grid(row=6, column=0, sticky="w")
        Checkbutton(save_window, text="Normalized Raw Counts", variable=self.save_norm_spectra_var).grid(row=7, column=0, sticky="w")
        Checkbutton(save_window, text="Normalized Flux Density", variable=self.save_norm_luminescence_flux_density_var).grid(row=8,column=0,sticky="w")
        Checkbutton(save_window, text="Normalized Absolute Gradient", variable=self.save_norm_absolute_gradient_var).grid(row=9, column=0, sticky="w")
        Checkbutton(save_window, text="Combined Plots", variable=self.save_combined_var).grid(row=10, column=0, sticky="w")

        # Entry fields for file names
        for i, (plot_name, var) in enumerate(self.save_names.items()):
            Label(save_window, text=f"{plot_name}:").grid(row=i + 1, column=1, sticky="e", padx=10)
            Entry(save_window, textvariable=var, width=30).grid(row=i + 1, column=2, sticky="w")

        Button(save_window, text="Save", command=lambda: self.save_plots(save_window)).grid(row=11, column=1, columnspan=2, pady=10)

    def save_plots(self, save_window):
        save_window.destroy()
        save_dir = filedialog.askdirectory(title="Select Save Directory")
        if not save_dir:
            return

        if self.save_raw_var.get() and self.contains_counts_flag:
            self.save_plot(self.plot_spectra, save_dir, self.save_names["Raw Data"].get())
        if self.save_luminescence_flux_density_var.get() and self.contains_flux_flag:
            self.save_plot(self.plot_luminescence_flux_density, save_dir, self.save_names["Luminescence Flux Density"].get())
        if self.save_absolute_gradient_var.get() and not self.is_single_measurement_flag:
            self.save_plot(self.plot_absolute_gradient, save_dir, self.save_names["Absolute Gradient"].get())
        if self.save_log_spectra_var.get() and self.contains_counts_flag:
            self.save_plot(self.plot_log_spectra, save_dir, self.save_names["Log Spectra"].get())
        if self.save_log_luminescence_flux_density_var.get() and self.contains_flux_flag:
            self.save_plot(self.plot_log_luminescence_flux_density, save_dir, self.save_names["Log Luminescence Flux Density"].get())
        if self.save_log_absolute_gradient_var.get() and not self.is_single_measurement_flag:
            self.save_plot(self.plot_log_absolute_gradient, save_dir, self.save_names["Log Absolute Gradient"].get())
        if self.save_norm_spectra_var.get() and self.contains_counts_flag:
            self.save_plot(self.plot_norm_spectra, save_dir, self.save_names["Normalized Spectra"].get())
        if self.save_norm_luminescence_flux_density_var.get() and self.contains_flux_flag:
            self.save_plot(self.plot_norm_luminescence_flux_density, save_dir, self.save_names["Normalized Luminescence Flux Density"].get())
        if self.save_norm_absolute_gradient_var.get() and not self.is_single_measurement_flag:
            self.save_plot(self.plot_norm_absolute_gradient, save_dir, self.save_names["Normalized Absolute Gradient"].get())
        if self.save_combined_var.get():
            self.save_combined_plot(save_dir)

    def save_plot(self, plot_function, save_dir, plot_name):
        fig, ax = plt.subplots(figsize=(8, 5))
        plot_function(ax)
        fig.savefig(f"{save_dir}/{plot_name}.png")
        plt.close(fig)
        print(f"Saved {plot_name}")

    def save_combined_plot(self, save_dir):
        fig, axes = plt.subplots(3, 3, figsize=(18, 6), constrained_layout=True)
        self.plot_spectra(axes[0, 0])
        self.plot_luminescence_flux_density(axes[0, 1])
        self.plot_absolute_gradient(axes[0, 2])
        self.plot_log_spectra(axes[1, 0])
        self.plot_log_luminescence_flux_density(axes[1, 1])
        self.plot_log_absolute_gradient(axes[1, 2])
        self.plot_norm_spectra(axes[2, 0])
        self.plot_norm_luminescence_flux_density(axes[2, 1])
        self.plot_norm_absolute_gradient(axes[2, 2])

        fig.tight_layout()
        fig.savefig(f"{save_dir}/{self.save_names['Combined'].get()}.png")
        plt.close(fig)
        print("Saved combined plots.")

    def show_info(self, plot_type):
        info = {
            "Spectra": (
                "Spectra:\n"
                "The spectra plot shows the photoluminescence (PL) signal intensity as a function "
                "of wavelength for one or multiple time intervals. The logarithmic plot shows the same data on a logarithmic scale.\n\n"
                "The raw data is taken from the PL measurements. \n\n"
                "X-axis: Wavelength (nm)\nY-axis: Intensity (a.u.)"
            ),
            "Flux Density": (
                "Flux Density:\n"
                "The luminescence flux density plot shows the lumiescence flux density in photons per second per square centimeter as a function "
                "of wavelength for one or multiple time intervals. The logarithmic plot shows the same data on a logarithmic scale.\n\n"
                "The luminescence flux density data is taken from the PL measurements.\n\n"
                "X-axis: Wavelength (nm)\nY-axis: Luminescence flux density (photons/s cm^2)"
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
        }
        messagebox.showinfo(f"Info: {plot_type}", info.get(plot_type, "No information available."))

    def plot_in_window(self):
        # Check if a file has been loaded with valid data
        if len(self.wavelength) < 2:
            messagebox.showerror("Error", "Please load a file first.")
            return

        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        # Enlarge plot window for better visualization
        self.plot_frame.config(width=1800, height=1200)

        fig, axes = plt.subplots(3, 3, figsize=(15, 10), constrained_layout=True)

        # Clear all axes and turn off initially
        for ax in axes.flatten():
            ax.clear()
            ax.axis('off')  # Turn off the axis to make the plot area completely blank

        # Plot in fixed positions
        if self.plot_spectra_var.get():
            if self.raw_counts is not None and self.raw_counts.size > 0 and self.contains_counts_flag:
                axes[0, 0].axis('on')  # Turn on axis for this plot
                self.plot_spectra(axes[0, 0])
            else:
                print("Skipping plot_spectra")  # Debug print
                axes[0, 0].axis('off')

        if self.plot_luminescence_flux_density_var.get():
            if self.luminescence_flux_density is not None and self.luminescence_flux_density.size > 0 and self.contains_flux_flag:
                axes[0, 1].axis('on')  # Turn on axis for this plot
                self.plot_luminescence_flux_density(axes[0, 1])
            else:
                print("Skipping plot_luminescence_flux_density")  # Debug print
                axes[0, 1].axis('off')

        if self.plot_absolute_gradient_var.get():
            if not self.is_single_measurement_flag:
                axes[0, 2].axis('on')  # Turn on axis for this plot
                self.plot_absolute_gradient(axes[0, 2])
            else:
                print("Skipping plot_absolute_gradient")  # Debug print
                axes[0, 2].axis('off')

        if self.plot_log_spectra_var.get():
            if self.raw_counts is not None and self.raw_counts.size > 0 and self.contains_counts_flag:
                axes[1, 0].axis('on')  # Turn on axis for this plot
                self.plot_log_spectra(axes[1, 0])
            else:
                print("Skipping plot_log_spectra")  # Debug print
                axes[1, 0].axis('off')

        if self.plot_log_luminescence_flux_density_var.get():
            if self.luminescence_flux_density is not None and self.luminescence_flux_density.size > 0 and self.contains_flux_flag:
                axes[1, 1].axis('on')  # Turn on axis for this plot
                self.plot_log_luminescence_flux_density(axes[1, 1])
            else:
                print("Skipping plot_log_luminescence_flux_density")  # Debug print
                axes[1, 1].axis('off')

        if self.plot_log_absolute_gradient_var.get():
            if not self.is_single_measurement_flag:
                axes[1, 2].axis('on')  # Turn on axis for this plot
                self.plot_log_absolute_gradient(axes[1, 2])
            else:
                print("Skipping plot_log_absolute_gradient")  # Debug print
                axes[1, 2].axis('off')

        if self.plot_norm_spectra_var.get():
            if self.raw_counts is not None and self.raw_counts.size > 0 and self.contains_counts_flag:
                axes[2, 0].axis('on')  # Turn on axis for this plot
                self.plot_norm_spectra(axes[2, 0])
            else:
                print("Skipping plot_normalized_spectra")  # Debug print
                axes[2, 0].axis('off')

        if self.plot_norm_luminescence_flux_density_var.get():
            if self.luminescence_flux_density is not None and self.luminescence_flux_density.size > 0 and self.contains_flux_flag:
                axes[2, 1].axis('on')  # Turn on axis for this plot
                self.plot_norm_luminescence_flux_density(axes[2, 1])
            else:
                print("Skipping plot_normalized_luminescence_flux_density")  # Debug print
                axes[2, 1].axis('off')

        if self.plot_norm_absolute_gradient_var.get():
            if not self.is_single_measurement_flag:
                axes[2, 2].axis('on')  # Turn on axis for this plot
                self.plot_norm_absolute_gradient(axes[2, 2])
            else:
                print("Skipping plot_normalized_absolute_gradient")  # Debug print
                axes[2, 2].axis('off')

        canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        toolbar = NavigationToolbar2Tk(canvas, self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack()

    def plot_spectra(self, ax):
        colors = plt.cm.viridis(np.linspace(0, 1, len(self.time_checkboxes)))

        if self.raw_counts is None or self.wavelength is None or self.relative_times is None or self.contains_counts_flag == False:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Raw PL Spectra")
            return

        selected_indices = [i for var, i in self.time_checkboxes if var.get()]
        colors = plt.cm.viridis(np.linspace(0, 1, len(selected_indices)))
        if self.is_single_measurement_flag:
            ax.plot(self.wavelength, self.raw_counts[:, 0], color="blue", label="Data")
        else:
            for idx, i in enumerate(selected_indices):
                ax.plot(self.wavelength, self.raw_counts[:, i], color=colors[idx], label=f"{self.relative_times[i]}")

        if self.show_legend_var.get():
            ax.legend()
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Counts")
        ax.set_title("Raw PL Spectra")

    def plot_luminescence_flux_density(self, ax):
        colors = plt.cm.viridis(np.linspace(0, 1, len(self.time_checkboxes)))

        if self.luminescence_flux_density is None or self.wavelength is None or self.relative_times is None or self.contains_flux_flag == False:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Luminescence Flux Density")
            return

        selected_indices = [i for var, i in self.time_checkboxes if var.get()]
        colors = plt.cm.viridis(np.linspace(0, 1, len(selected_indices)))
        if self.is_single_measurement_flag:
            ax.plot(self.wavelength, self.luminescence_flux_density[:, 0], color="green", label="Data")
        else:
            for idx, i in enumerate(selected_indices):
                ax.plot(self.wavelength, self.luminescence_flux_density[:, i], color=colors[idx], label=f"{self.relative_times[i]}")

        if self.show_legend_var.get():
            ax.legend()
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Flux Density (photons/s/cm²/nm))")
        ax.set_title("Luminescence Flux Density")

    def plot_log_luminescence_flux_density(self, ax):
        colors = plt.cm.viridis(np.linspace(0, 1, len(self.time_checkboxes)))

        if self.luminescence_flux_density is None or self.wavelength is None or self.relative_times is None or self.contains_flux_flag == False:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Luminescence Flux Density (log)")
            return

        selected_indices = [i for var, i in self.time_checkboxes if var.get()]
        colors = plt.cm.viridis(np.linspace(0, 1, len(selected_indices)))
        if self.is_single_measurement_flag:
            ax.plot(self.wavelength, self.luminescence_flux_density[:, 0], color="green", label="Data")
        else:
            for idx, i in enumerate(selected_indices):
                ax.plot(self.wavelength, self.luminescence_flux_density[:, i], color=colors[idx], label=f"{self.relative_times[i]}")

        if self.show_legend_var.get():
            ax.legend()
        ax.set_yscale('log')
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Flux Density (photons/s/cm²/nm))")
        ax.set_title("Luminescence Flux Density (log)")

    def plot_norm_luminescence_flux_density(self, ax):
        colors = plt.cm.viridis(np.linspace(0, 1, len(self.time_checkboxes)))

        if self.luminescence_flux_density is None or self.wavelength is None or self.relative_times is None or self.contains_flux_flag == False:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Normalized Luminescence Flux Density")
            return

        selected_indices = [i for var, i in self.time_checkboxes if var.get()]
        colors = plt.cm.viridis(np.linspace(0, 1, len(selected_indices)))
        if self.is_single_measurement_flag:
            ax.plot(self.wavelength, self.luminescence_flux_density[:, 0]/self.luminescence_flux_density[:, 0].max(), color="green", label="Data")
        else:
            for idx, i in enumerate(selected_indices):
                ax.plot(self.wavelength, self.luminescence_flux_density[:, i]/self.luminescence_flux_density.max(), color=colors[idx], label=f"{self.relative_times[i]}")

        if self.show_legend_var.get():
            ax.legend()
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Normalized Flux Density (a. u.)")
        ax.set_title("Normalized Luminescence Flux Density")

    def plot_log_spectra(self, ax):
        colors = plt.cm.viridis(np.linspace(0, 1, len(self.time_checkboxes)))

        if self.raw_counts is None or self.wavelength is None or self.relative_times is None or self.contains_counts_flag == False:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Raw PL Spectra (log)")
            return

        selected_indices = [i for var, i in self.time_checkboxes if var.get()]
        colors = plt.cm.viridis(np.linspace(0, 1, len(selected_indices)))
        if self.is_single_measurement_flag:
            ax.plot(self.wavelength, self.raw_counts[:, 0], color="blue", label="Data")
        else:
            for idx, i in enumerate(selected_indices):
                ax.plot(self.wavelength, self.raw_counts[:, i], color=colors[idx], label=f"{self.relative_times[i]}")

        if self.show_legend_var.get():
            ax.legend()
        ax.set_yscale('log')
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Counts")
        ax.set_title("Raw PL Spectra (log)")

    def plot_norm_spectra(self, ax):
        colors = plt.cm.viridis(np.linspace(0, 1, len(self.time_checkboxes)))
        if self.raw_counts is None or self.wavelength is None or self.relative_times is None or self.contains_counts_flag == False:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Normalized Raw PL Spectra")
            return

        selected_indices = [i for var, i in self.time_checkboxes if var.get()]
        colors = plt.cm.viridis(np.linspace(0, 1, len(selected_indices)))
        if self.is_single_measurement_flag:
            ax.plot(self.wavelength, self.raw_counts[:, 0]/self.raw_counts[:, 0].max(), color="blue", label="Data")
        else:
            for idx, i in enumerate(selected_indices):
                ax.plot(self.wavelength, self.raw_counts[:, i]/self.raw_counts.max(), color=colors[idx], label=f"{self.relative_times[i]}")

        if self.show_legend_var.get():
            ax.legend()
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Normalized Intensity (a. u.)")
        ax.set_title("Normalized Raw PL Spectra")

    def plot_absolute_gradient(self, ax):
        if self.is_single_measurement_flag:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Absolute Gradient")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Gradient Value")
        else:
            if self.contains_counts_flag:
                data = self.raw_counts
            elif self.contains_flux_flag:
                data = self.luminescence_flux_density
            else:
                ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
                ax.set_title("Absolute Gradient")
                ax.set_xlabel("Wavelength (nm)")
                ax.set_ylabel("Gradient Value")
                return

            gradient = np.gradient(data, axis=0)
            mean_gradient = np.mean(np.abs(gradient), axis=1)
            ax.plot(self.wavelength, mean_gradient, color="blue", label="Absolute Gradient")
            ax.set_title("Absolute Gradient")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Gradient Value")
            if self.show_legend_var.get():
                ax.legend()

    def plot_log_absolute_gradient(self, ax):
        if self.is_single_measurement_flag:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Log Absolute Gradient")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Gradient Value (log)")
            return

        if self.contains_counts_flag:
            data = self.raw_counts
        elif self.contains_flux_flag:
            data = self.luminescence_flux_density
        else:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Log Absolute Gradient")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Gradient Value (log)")
            return

        gradient = np.gradient(data, axis=0)
        mean_gradient = np.mean(np.abs(gradient), axis=1)
        ax.plot(self.wavelength, mean_gradient, color="blue", label="Absolute Gradient")
        ax.set_title("Log Absolute Gradient")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Gradient Value (log)")
        ax.set_yscale('log')
        if self.show_legend_var.get():
            ax.legend()

    def plot_norm_absolute_gradient(self, ax):
        if self.is_single_measurement_flag:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Normalized Absolute Gradient")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Normalized Gradient (a. u.)")
            return

        if self.contains_counts_flag:
            data = self.raw_counts
        elif self.contains_flux_flag:
            data = self.luminescence_flux_density
        else:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Normalized Absolute Gradient")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Normalized Gradient (a. u.)")
            return

        gradient = np.gradient(data, axis=0)
        mean_gradient = np.mean(np.abs(gradient), axis=1)
        ax.plot(self.wavelength, mean_gradient/mean_gradient.max(), color="blue", label="Absolute Gradient")
        ax.set_title("Normalzied Absolute Gradient")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Gradient Value")
        if self.show_legend_var.get():
            ax.legend()

    def update_raw_plot(self):
        self.plot_in_window()

if __name__ == "__main__":
    root = Tk()
    app = PLAnalysisApp(root)
    root.mainloop()

