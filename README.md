Photoluminescence (PL) Analysis Tool

PL Analysis Tool is a Python-based GUI application for analyzing photoluminescence measurements. This tool supports both continuous and non-continuous measurement data, providing comprehensive visualization of spectra, instability gradients, and relative changes. It offers interactive data selection, plot customization, and metadata display.


Features:

Load and filter photoluminescence measurement data from text files.
Visualize raw spectra, luminescence flux data, instability gradients, and relative changes.
Supports both linear and logarithmic plots.
Save plots in PNG format.


Requirements:

To run the PL Analysis Tool, you need to have Python 3.8 or later installed. Additionally, the following Python libraries are required:

pandas

matplotlib

numpy

Pillow (PIL)

tkinter (usually comes with Python installation)


Dependecies:

Install Python 3.8 or later. You can download Python from the official website: Python.org

Install required libraries

Open a terminal or command prompt and run the following commands to install the necessary libraries:

pip install pandas matplotlib numpy Pillow

Note: tkinter is usually included with your Python installation. If it's not, you may need to install it separately depending on your operating system.


Installation:

Clone the repository:

open bash

git clone https://github.com/Failistig/PL_Analysis_Tool.git

cd PL_Analysis_Tool

Run the application: python PL_Peak_Characterization.py


Using the GUI:

Click on "Load File" to select a text file containing photoluminescence measurement data.
Select the types of plots you want to visualize (Spectra, Instability Gradient, Relative Change, and their logarithmic counterparts).
Click on "Plot All" to generate the plots.
Use the "Save Plots" button to save the current plots as PNG files.
.
.
.
.
.
.
