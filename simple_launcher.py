#!/usr/bin/env python3
"""
Simple Jones Calculator Launcher

A simple working launcher for the Jones calculator functionality.
"""

import sys

def main():
    """Simple launcher menu"""
    print("="*60)
    print("TMM JONES CALCULATOR LAUNCHER")
    print("="*60)
    print("Available options:")
    print("1. Run tutorial (simple_jones_tutorial.py)")
    print("2. Run fine grid demo")
    print("3. Run custom calculation")
    print("="*60)
    
    try:
        choice = input("Enter your choice (1-3): ").strip()
        
        if choice == "1":
            print("Launching tutorial...")
            import subprocess
            subprocess.run([sys.executable, "simple_jones_tutorial.py"])
            
        elif choice == "2":
            print("Launching fine grid demo...")
            import subprocess
            subprocess.run([sys.executable, "fine_grid_demo.py"])
            
        elif choice == "3":
            print("For custom calculations, use:")
            print("  from tmm.calculations import Jones3DAnalyzer, CalculationParams")
            print("  See simple_jones_tutorial.py for examples")
            
        else:
            print("Invalid choice. Please run again.")
            
    except KeyboardInterrupt:
        print("\nExiting...")
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main() 