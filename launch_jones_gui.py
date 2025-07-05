#!/usr/bin/env python3
"""
Jones Unified GUI Launcher

Simple launcher for the unified Jones calculator GUI.
This provides an easy way to start the graphical interface.

Usage: python launch_jones_gui.py
"""

import sys
import tkinter as tk
from tkinter import messagebox

def check_dependencies():
    """Check if required dependencies are available"""
    try:
        import numpy
        import matplotlib
        from tmm.gui.jones_unified_gui import JonesUnifiedGUI
        return True
    except ImportError as e:
        messagebox.showerror("Missing Dependencies", 
                           f"Required dependencies not found:\n{str(e)}\n\n"
                           "Please install the TMM package and its dependencies.")
        return False

def main():
    """Launch the Jones unified GUI"""
    print("Launching Jones Unified Calculator GUI...")
    
    # Create root window
    root = tk.Tk()
    
    # Check dependencies
    if not check_dependencies():
        root.destroy()
        return
    
    try:
        # Import and create GUI
        from tmm.gui.jones_unified_gui import JonesUnifiedGUI
        app = JonesUnifiedGUI(root)
        
        print("GUI started successfully!")
        print("Use the GUI to:")
        print("• Configure layer structures")
        print("• Set calculation parameters")
        print("• Choose grid resolution (coarse/fine)")
        print("• Run 3D E-kx-ky analysis")
        print("• Visualize results")
        
        # Start GUI event loop
        root.mainloop()
        
    except Exception as e:
        messagebox.showerror("GUI Error", 
                           f"Failed to start GUI:\n{str(e)}")
        print(f"Error: {str(e)}")
    
    print("GUI closed.")

if __name__ == "__main__":
    main() 