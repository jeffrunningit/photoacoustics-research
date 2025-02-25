#!/usr/bin/env python3
import os
import subprocess
import sys

py_files = ["./Figures/transient_grating_plot.py",
            "./Data/peak fluence.py",
            "./Data/grating spectrum/spectrometer plot.py",
            "./Data/grating spectrum/18012024/white light spectroscopy plots.py",
            "./Data/grating 750nm measurement/grating 750nm probe measurement plots.py",
            "./Data/grating 750nm measurement/750nm -1 to +1 data/plot_script.py",
            "./Data/grating 750nm measurement/750nm -3 to +3 data/plot_script.py",
            "./Data/AFM/crosssection_plot.py",
            "./Data/grating whitelight measurement/Measurements/whitelight pumpprobe 2dplot.py",
            "./Data/Resonance angle calculation data/Al resonance angle calculation.py"]

def display_menu(py_files):
    print("\nPython files available:")
    for i, file in enumerate(py_files, start=1):
        print(f"{i}. {file}")
    print("0. Exit")

def get_user_choice(max_choice):
    while True:
        try:
            choice = int(input("Select a file to run (by number): "))
            if 0 <= choice <= max_choice:
                return choice
            else:
                print(f"Please enter a number between 0 and {max_choice}.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def main():
    while True:
        display_menu(py_files)
        choice = get_user_choice(len(py_files))
        if choice == 0:
            print("Exiting.")
            break

        file_to_run = py_files[choice - 1]
        print(f"\nRunning: {file_to_run}\n")
        try:
            # Run the selected Python file in a subprocess.
            subprocess.run([sys.executable, os.path.basename(file_to_run)], check=True, cwd=os.path.dirname(file_to_run))
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running {file_to_run}: {e}")
        except KeyboardInterrupt:
            print("\nInterrupted by user.")
            break

        print("\nFinished running the script.\n")

if __name__ == "__main__":
    main()