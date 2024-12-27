import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tkinter import Tk, filedialog, Button, Label, Checkbutton, IntVar, Entry, Frame, Toplevel, PhotoImage, StringVar, messagebox, HORIZONTAL, VERTICAL, RIGHT, BOTTOM, Y, X, LEFT, BOTH, Text
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from PIL import Image, ImageTk
from datetime import datetime
from tkinter import Scrollbar, Canvas

global matching_col_raw
global col_index_raw
global matching_col_flux
global col_index_flux

class PLAnalysisApp:
    def __init__(self, master):
        self.master = master
        self.master.title("PL Analysis Tool")
        self.master.geometry("1920x1080")

        # Create a canvas and a frame to hold the widgets
        self.canvas = Canvas(master)
        self.scrollable_frame = Frame(self.canvas)

        # Add scrollbars
        self.v_scrollbar = Scrollbar(master, orient="vertical", command=self.canvas.yview)
        self.h_scrollbar = Scrollbar(master, orient="horizontal", command=self.canvas.xview)

        self.canvas.configure(yscrollcommand=self.v_scrollbar.set)
        self.canvas.configure(xscrollcommand=self.h_scrollbar.set)

        self.v_scrollbar.pack(side="right", fill="y")
        self.h_scrollbar.pack(side="bottom", fill="x")
        self.canvas.pack(side="left", fill="both", expand=True)

        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        self.scrollable_frame.bind("<Configure>", self.on_frame_configure)

        # Bind the window close event to a custom method
        self.master.protocol("WM_DELETE_WINDOW", self.on_close)

        # Variables
        self.file_path = None
        self.data = None
        self.wavelength = None
        self.raw_counts = None
        self.time_labels = None
        self.relative_times = None
        self.measurement_date = None

        self.plot_spectra_var = IntVar(value=1)
        self.plot_luminescence_flux_density_var = IntVar(value=1)
        self.plot_absolute_gradient_var = IntVar(value=1)
        self.plot_relative_change_var = IntVar(value=1)
        self.plot_log_spectra_var = IntVar(value=1)
        self.plot_log_luminescence_flux_density_var = IntVar(value=1)
        self.plot_log_absolute_gradient_var = IntVar(value=1)
        self.plot_log_relative_change_var = IntVar(value=1)
        self.show_legend_var = IntVar(value=1)

        # Save options
        self.save_raw_var = IntVar(value=1)
        self.save_luminescence_flux_density_var = IntVar(value=1)
        self.save_absolute_gradient_var = IntVar(value=1)
        self.save_relative_change_var = IntVar(value=1)
        self.save_log_spectra_var = IntVar(value=1)
        self.save_log_luminescence_flux_density_var = IntVar(value=1)
        self.save_log_absolute_gradient_var = IntVar(value=1)
        self.save_log_relative_change_var = IntVar(value=1)
        self.save_combined_var = IntVar(value=1)
        self.save_names = {
            "Raw Data": StringVar(value="Raw_Data"),
            "Luminescence Flux Density": StringVar(value="Luminescence_Flux_Density"),
            "Absolute Gradient": StringVar(value="Absolute_Gradient"),
            "Relative Change": StringVar(value="Relative_Change"),
            "Log Spectra": StringVar(value="Log_Spectra"),
            "Log Luminescence Flux Density": StringVar(value="Log_Luminescence_Flux_Density"),
            "Log Absolute Gradient": StringVar(value="Log_Absolute_Gradient"),
            "Log Relative Change": StringVar(value="Log_Relative_Change"),
            "Combined": StringVar(value="Combined_Plots"),
        }

        self.time_checkboxes = []
        self.logo_image = None  # Placeholder for the logo

        self.setup_gui()

    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(int(-1*(event.delta/120)), "units")

    def _on_shiftmousewheel(self, event):
        self.canvas.xview_scroll(int(-1*(event.delta/120)), "units")

    def on_frame_configure(self, event):
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

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
        Label(self.scrollable_frame, text="PL Analysis Tool", font=("Arial", 16)).grid(row=0, column=0, columnspan=3,
                                                                                       pady=10)

        # Load Logo in Top-Right Corner
        self.load_logo("./hzb_logo.jpg")

        # File Load Section
        Button(self.scrollable_frame, text="Load File", width=15, command=self.load_file).grid(row=1, column=0, padx=10,pady=5, sticky="w")

        self.file_path_entry = Entry(self.scrollable_frame, width=80, state='disabled')
        self.file_path_entry.grid(row=1, column=1, columnspan=2, sticky="w")

        Label(self.scrollable_frame, text="Timestamp:").grid(row=0, column=0, padx=10, pady=5, sticky="nw")
        self.date_entry = Entry(self.scrollable_frame, width=18, state='disabled')
        self.date_entry.grid(row=0, column=0, padx=10, pady=5, columnspan=1, sticky="w")

        # Buttons Section
        Button(self.scrollable_frame, text="Save Plots", width=15, command=self.save_plots_dialog).grid(row=2, column=0,
                                                                                                        padx=10, pady=5,
                                                                                                        sticky="w")
        Button(self.scrollable_frame, text="Plot All", width=15, command=self.plot_in_window).grid(row=2, column=1,
                                                                                                   pady=5, sticky="w")

        # Plot Selection
        Label(self.scrollable_frame, text="Select Plots:").grid(row=3, column=0, pady=5, sticky="w", padx=10)

        # First column of checkboxes
        Checkbutton(self.scrollable_frame, text="Spectra", variable=self.plot_spectra_var).grid(row=4, column=0,
                                                                                                sticky="w", padx=10)
        Checkbutton(self.scrollable_frame, text="Flux Density", variable=self.plot_luminescence_flux_density_var).grid(
            row=5, column=0, sticky="w", padx=10)
        Checkbutton(self.scrollable_frame, text="Absolute Gradient", variable=self.plot_absolute_gradient_var).grid(
            row=6, column=0, sticky="w", padx=10)
        Checkbutton(self.scrollable_frame, text="Relative Change", variable=self.plot_relative_change_var).grid(row=7,
                                                                                                                column=0,
                                                                                                                sticky="w",
                                                                                                                padx=10)

        # Second column of checkboxes
        Checkbutton(self.scrollable_frame, text="Log Spectra", variable=self.plot_log_spectra_var).grid(row=4, column=1,
                                                                                                        sticky="w")
        Checkbutton(self.scrollable_frame, text="Log Flux Density",
                    variable=self.plot_log_luminescence_flux_density_var).grid(row=5, column=1, sticky="w")
        Checkbutton(self.scrollable_frame, text="Log Absolute Gradient",
                    variable=self.plot_log_absolute_gradient_var).grid(row=6, column=1, sticky="w")
        Checkbutton(self.scrollable_frame, text="Log Relative Change", variable=self.plot_log_relative_change_var).grid(
            row=7, column=1, sticky="w")

        # Control buttons
        Button(self.scrollable_frame, text="Trim Data", width=20, command=self.open_trim_options_window).grid(row=4,
                                                                                                              column=2,
                                                                                                              sticky="w")
        Button(self.scrollable_frame, text="Filter Data", width=20, command=self.open_filter_options_window).grid(row=5,
                                                                                                                  column=2,
                                                                                                                  sticky="w")
        Button(self.scrollable_frame, text="Show Metadata", width=20, command=self.show_metadata).grid(row=6, column=2,
                                                                                                       sticky="w")

        # Info Buttons
        Button(self.scrollable_frame, text="Info: Spectra", width=20, anchor="w",
               command=lambda: self.show_info("Spectra")).grid(row=4, column=3, sticky="w")
        Button(self.scrollable_frame, text="Info: Flux Density", width=20, anchor="w",
               command=lambda: self.show_info("Flux Density")).grid(row=5, column=3, sticky="w")
        Button(self.scrollable_frame, text="Info: Absolute Gradient", width=20, anchor="w",
               command=lambda: self.show_info("Absolute Gradient")).grid(row=6, column=3, sticky="w")
        Button(self.scrollable_frame, text="Info: Relative Change", width=20, anchor="w",
               command=lambda: self.show_info("Relative Change")).grid(row=7, column=3, sticky="w")

        # Legend Toggle
        Checkbutton(self.scrollable_frame, text="Show Legend", variable=self.show_legend_var).grid(row=8, column=0,
                                                                                                   sticky="w", padx=10)

        # Plot Frame
        self.plot_frame = Frame(self.scrollable_frame, width=1200, height=800, bg="white")
        self.plot_frame.grid(row=9, column=0, columnspan=4, pady=10)

        # Raw Data Selection
        self.time_check_frame = Frame(self.scrollable_frame)
        self.time_check_frame.grid(row=11, column=0, columnspan=4, pady=5)

        Button(self.scrollable_frame, text="Update Plots", width=15, command=self.update_raw_plot).grid(row=10,
                                                                                                        column=1)
        # Select All Checkbox
        self.select_all_var = IntVar(value=1)  # Default: All selected
        Checkbutton(self.scrollable_frame,text="Toggle All Times",variable=self.select_all_var,command=self.toggle_all_time_checkboxes).grid(row=10 , column=2)


    def load_logo(self, file_path):
        """Loads and places the logo image in the top-right corner."""
        try:
            image = Image.open(file_path)
            image = image.resize((200, 100), Image.Resampling.LANCZOS)  # Resize the logo to 2x size
            self.logo_image = ImageTk.PhotoImage(image)

            # Display the logo
            Label(self.scrollable_frame, image=self.logo_image).grid(row=0, column=3, sticky="ne")
        except Exception as e:
            print(f"Error loading logo: {e}")

    def update_date_label(self):
        """Update the date label with the measurement date."""
        if self.measurement_date:
            self.date_label.config(text=f"Measurement Date: {self.measurement_date}")

    def open_trim_options_window(self):
        trim_window = Toplevel(self.master)
        trim_window.title("Trim Data Options")

        Label(trim_window, text="Trim Data Options", font=("Arial", 14)).pack(pady=10)

        # Add trim options here
        trim_option1 = IntVar()
        Checkbutton(trim_window, text="Trim Option 1", variable=trim_option1).pack(anchor="w")
        Button(trim_window, text="Info", command=lambda: self.show_option_info("Trim Option 1 Info")).pack(anchor="w")

        trim_option2 = IntVar()
        Checkbutton(trim_window, text="Trim Option 2", variable=trim_option2).pack(anchor="w")
        Button(trim_window, text="Info", command=lambda: self.show_option_info("Trim Option 2 Info")).pack(anchor="w")

        # Add more trim options as needed

        Button(trim_window, text="Apply", command=lambda: self.apply_trim_options([trim_option1, trim_option2])).pack(
            pady=10)

    def open_filter_options_window(self):
        filter_window = Toplevel(self.master)
        filter_window.title("Filter Data Options")

        Label(filter_window, text="Filter Data Options", font=("Arial", 14)).pack(pady=10)

        # Add filter options here
        filter_option1 = IntVar()
        Checkbutton(filter_window, text="Filter Option 1", variable=filter_option1).pack(anchor="w")
        Button(filter_window, text="Info", command=lambda: self.show_option_info("Filter Option 1 Info")).pack(
            anchor="w")

        filter_option2 = IntVar()
        Checkbutton(filter_window, text="Filter Option 2", variable=filter_option2).pack(anchor="w")
        Button(filter_window, text="Info", command=lambda: self.show_option_info("Filter Option 2 Info")).pack(
            anchor="w")

        # Add more filter options as needed

        Button(filter_window, text="Apply",
               command=lambda: self.apply_filter_options([filter_option1, filter_option2])).pack(pady=10)

    def show_option_info(self, info_text):
        messagebox.showinfo("Option Info", info_text)

    def apply_trim_options(self, options):
        # Apply the selected trim options to the data
        for option in options:
            if option.get():
                print(f"Applying {option}...")  # Replace with actual logic
        self.update_raw_plot()

    def apply_filter_options(self, options):
        # Apply the selected filter options to the data
        for option in options:
            if option.get():
                print(f"Applying {option}...")  # Replace with actual logic
        self.update_raw_plot()

    def show_metadata(self):
        """Display metadata in a new scrollable window."""
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

        # Scan the file to find the header line
        header_line = "Wavelength (nm)"
        skip_rows = 0

        def read_file_with_encoding(encoding):
            nonlocal skip_rows
            with open(self.file_path, 'r', encoding=encoding) as file:
                for line in file:
                    if header_line in line:
                        break
                    skip_rows += 1

            self.data = pd.read_csv(self.file_path, sep="\t", skiprows=skip_rows + 1, encoding=encoding)
            self.headers = pd.read_csv(self.file_path, sep="\t", skiprows=skip_rows, encoding=encoding)

        try:
            # Try reading the file with UTF-8 encoding first
            read_file_with_encoding('utf-8')
        except UnicodeDecodeError:
            # If UTF-8 fails, try reading the file with ISO-8859-1 encoding
            skip_rows = 0
            read_file_with_encoding('ISO-8859-1')

        try:
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
                if self.contains_flux_flag == True:
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
            print( f"Luminescence flux density shape: {self.luminescence_flux_density.shape if self.luminescence_flux_density is not None else 'None'}")
            print(f"contains_counts_flag: {self.contains_counts_flag}")
            print(f"contains_flux_flag: {self.contains_flux_flag}")

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

    def calculate_relative_times(self, raw_times):
        # Check if it's a single measurement first
        if len(raw_times) <= 1:
            return None

        try:
            # Check if time contains date
            if " " in raw_times[0]:
                time_format = "%d.%m.%Y %H:%M:%S"
            else:
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
        max_per_row = 12

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
        Checkbutton(save_window, text="Relative Change", variable=self.save_relative_change_var).grid(row=4, column=0, sticky="w")
        Checkbutton(save_window, text="Log Raw Counts", variable=self.save_log_spectra_var).grid(row=5, column=0, sticky="w")
        Checkbutton(save_window, text="Log Flux Density", variable=self.save_log_luminescence_flux_density_var).grid(row=6,column=0,sticky="w")
        Checkbutton(save_window, text="Log Absolute Gradient", variable=self.save_log_absolute_gradient_var).grid(row=7, column=0, sticky="w")
        Checkbutton(save_window, text="Log Relative Change", variable=self.save_log_relative_change_var).grid(row=8, column=0, sticky="w")
        Checkbutton(save_window, text="Combined Plots", variable=self.save_combined_var).grid(row=9, column=0, sticky="w")

        # Entry fields for file names
        for i, (plot_name, var) in enumerate(self.save_names.items()):
            Label(save_window, text=f"{plot_name}:").grid(row=i + 1, column=1, sticky="e", padx=10)
            Entry(save_window, textvariable=var, width=30).grid(row=i + 1, column=2, sticky="w")

        Button(save_window, text="Save", command=lambda: self.save_plots(save_window)).grid(row=10, column=1, columnspan=2, pady=10)

    def save_plots(self, save_window):
        # Check if the required data is available before attempting to save the plots
        if self.plot_spectra_var.get() and (self.raw_counts is None or self.raw_counts.size == 0):
            messagebox.showerror("Error", "No raw data available. Please deselect 'Raw PL Spectra / Log Raw PL Spectra'.")
            return
        elif self.plot_luminescence_flux_density_var.get() and (
                self.luminescence_flux_density is None or self.luminescence_flux_density.size == 0):
            messagebox.showerror("Error", "No flux density data available. Please deselect 'Luminescence Flux Density / Log Luminescence Flux Density'.")
            return
        elif self.is_single_measurement_flag_var.get() == True:
            messagebox.showerror("Error", "Not a continous measurement, no Gradient or relative Change data available. Please deselect all the related plots.")
            return

        else:
            return

        save_window.destroy()
        save_dir = filedialog.askdirectory(title="Select Save Directory")
        if not save_dir:
            return

        if self.save_raw_var.get():
            self.save_plot(self.plot_spectra, save_dir, self.save_names["Raw Data"].get())
        if self.save_luminescence_flux_density_var.get():
            self.save_plot(self.plot_luminescence_flux_density, save_dir, self.save_names["Flux Density"].get())
        if self.save_absolute_gradient_var.get():
            self.save_plot(self.plot_absolute_gradient, save_dir, self.save_names["Absolute Gradient"].get())
        if self.save_relative_change_var.get():
            self.save_plot(self.plot_relative_change, save_dir, self.save_names["Relative Change"].get())
        if self.save_log_spectra_var.get():
            self.save_plot(self.plot_log_spectra, save_dir, self.save_names["Log Spectra"].get())
        if self.save_log_luminescence_flux_density_var.get():
            self.save_plot(self.plot_log_luminescence_flux_density, save_dir, self.save_names["Log Flux Density"].get())
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
        fig, axes = plt.subplots(2, 4, figsize=(14, 6), constrained_layout=True)
        self.plot_spectra(axes[0, 0])
        self.plot_luminescence_flux_density(axes[0, 1])
        self.plot_absolute_gradient(axes[0, 2])
        self.plot_relative_change(axes[0, 3])
        self.plot_log_spectra(axes[1, 0])
        self.plot_log_luminescence_flux_density(axes[1, 1])
        self.plot_log_absolute_gradient(axes[1, 2])
        self.plot_log_relative_change(axes[1, 3])
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
        # Check if a file has been loaded with valid data
        if len(self.wavelength) < 2:
            messagebox.showerror("Error", "Please load a file first.")
            return

        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        # Enlarge plot window for better visualization
        self.plot_frame.config(width=1800, height=800)

        fig, axes = plt.subplots(2, 4, figsize=(15, 6), constrained_layout=True)

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

        if self.plot_relative_change_var.get():
            if not self.is_single_measurement_flag:
                axes[0, 3].axis('on')  # Turn on axis for this plot
                self.plot_relative_change(axes[0, 3])
            else:
                print("Skipping plot_relative_change")  # Debug print
                axes[0, 3].axis('off')

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

        if self.plot_log_relative_change_var.get():
            if not self.is_single_measurement_flag:
                axes[1, 3].axis('on')  # Turn on axis for this plot
                self.plot_log_relative_change(axes[1, 3])
            else:
                print("Skipping plot_log_relative_change")  # Debug print
                axes[1, 3].axis('off')

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
        for idx, i in enumerate(selected_indices):
            ax.plot(self.wavelength, self.raw_counts[:, i], color=colors[idx], label=f"Time {i}")

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
        for idx, i in enumerate(selected_indices):
            ax.plot(self.wavelength, self.luminescence_flux_density[:, i], color=colors[idx], label=f"Time {i}")

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
        for idx, i in enumerate(selected_indices):
            ax.plot(self.wavelength, self.luminescence_flux_density[:, i], color=colors[idx], label=f"Time {i}")

        if self.show_legend_var.get():
            ax.legend()
        ax.set_yscale('log')
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Flux Density (photons/s/cm²/nm))")
        ax.set_title("Luminescence Flux Density (log)")

    def plot_log_spectra(self, ax):
        colors = plt.cm.viridis(np.linspace(0, 1, len(self.time_checkboxes)))

        if self.raw_counts is None or self.wavelength is None or self.relative_times is None or self.contains_counts_flag == False:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Raw PL Spectra (log)")
            return

        selected_indices = [i for var, i in self.time_checkboxes if var.get()]
        colors = plt.cm.viridis(np.linspace(0, 1, len(selected_indices)))
        for idx, i in enumerate(selected_indices):
            ax.plot(self.wavelength, self.raw_counts[:, i], color=colors[idx], label=f"Time {i}")

        if self.show_legend_var.get():
            ax.legend()
        ax.set_yscale('log')
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Counts")
        ax.set_title("Raw PL Spectra (log)")

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

    def plot_relative_change(self, ax):
        if self.is_single_measurement_flag:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Relative Change")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Relative Change Value")
            return

        if self.contains_counts_flag:
            data = self.raw_counts
        elif self.contains_flux_flag:
            data = self.luminescence_flux_density
        else:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Relative Change")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Relative Change Value")
            return

        shape_change = np.abs(np.diff(data, axis=1))
        epsilon = 1e-8
        quotient = shape_change / (np.abs(data[:, :-1]) + epsilon)
        average_quotient = np.mean(quotient, axis=1)
        ax.plot(self.wavelength[:-1], average_quotient[:len(self.wavelength) - 1], color="green",
                label="Relative Change")
        ax.set_title("Relative Change")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Relative Change Value")
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

    def plot_log_relative_change(self, ax):
        if self.is_single_measurement_flag:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Log Relative Change")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Relative Change Value (log)")
            return

        if self.contains_counts_flag:
            data = self.raw_counts
        elif self.contains_flux_flag:
            data = self.luminescence_flux_density
        else:
            ax.text(0.5, 0.5, "No Data to Plot", ha='center', va='center', transform=ax.transAxes)
            ax.set_title("Log Relative Change")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Relative Change Value (log)")
            return

        shape_change = np.abs(np.diff(data, axis=1))
        epsilon = 1e-8
        quotient = shape_change / (np.abs(data[:, :-1]) + epsilon)
        average_quotient = np.mean(quotient, axis=1)
        ax.plot(self.wavelength[:-1], average_quotient[:len(self.wavelength) - 1], color="green",
                label="Relative Change")
        ax.set_title("Log Relative Change")
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

