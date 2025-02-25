#!/usr/bin/env python3
import os
import subprocess
import sys

def find_py_files(root='.'):
    """Recursively find all Python (.py) files in the directory tree."""
    py_files = []
    for dirpath, _, filenames in os.walk(root):
        for filename in filenames:
            if filename.endswith('.py'):
                # Skip this file (the launcher) so it doesn't list itself.
                if os.path.abspath(os.path.join(dirpath, filename)) == os.path.abspath(__file__):
                    continue
                file_path = os.path.join(dirpath, filename)
                py_files.append(file_path)
    return py_files

def display_menu(py_files):
    print("\nPython files found:")
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
    py_files = find_py_files()
    if not py_files:
        print("No Python files found in the directory.")
        return

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
            subprocess.run([sys.executable, file_to_run], check=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running {file_to_run}: {e}")
        except KeyboardInterrupt:
            print("\nInterrupted by user.")
            break

        print("\nFinished running the script.\n")

if __name__ == "__main__":
    main()