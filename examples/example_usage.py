#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example Usage Script for LapTimeSim

This script demonstrates how to use the LapTimeSim modules to:
1. Process a vehicle specification file
2. Process a track data file
3. Run a lap time simulation
4. Analyze the results

Run this script from the examples/ directory.
"""

import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


def example_basic_simulation():
    """
    Basic example: Run a simulation with pre-processed model files.
    """
    print("=" * 60)
    print("Example 1: Basic Simulation")
    print("=" * 60)
    
    # This would normally import from LapSimV2
    # For demonstration, we'll show the workflow
    
    print("""
    To run a basic simulation:
    
    1. Ensure you have model files generated:
       - Model_<TrackName>_TrackSimV2.txt
       - Model_<CarName>_CarSimV2.txt
    
    2. In LapSimV2.py, set the file paths:
       nom_circuit, pays, tr = ouvertureFileTrack("path/to/track_model.txt")
       nom_voiture, veh = ouvertureFileCar("path/to/car_model.txt")
    
    3. Run the simulation:
       python src/LapSimV2.py
    
    The simulation will output:
    - Lap time
    - Speed profile plot
    - Track map with velocity coloring
    - Engine/gear data plot
    """)


def example_create_car_model():
    """
    Example: Create a vehicle model from specification file.
    """
    print("\n" + "=" * 60)
    print("Example 2: Creating a Vehicle Model")
    print("=" * 60)
    
    print("""
    Vehicle specification file format (tab-separated):
    
    ```
    My Racing Car
    Total Mass	1100	kg
    Front Mass Distribution	42	%
    Wheelbase	2420	mm
    Lift Coefficient CL	-2.8	-
    Drag Coefficient CD	-0.8	-
    CL Scale Multiplier	1	-
    CD Scale Multiplier	1	-
    Front Aero Distribution	50	%
    Frontal Area	1	m2
    Air Density	1.225	kg/m3
    Grip Factor Multiplier	1.5	-
    Tyre Radius	225	mm
    Rolling Resistance	-0.001	-
    Longitudinal Friction Coefficient	1.9	-
    Longitudinal Friction Load Rating	350	kg
    Longitudinal Friction Sensitivity	0.0001	1/N
    Lateral Friction Coefficient	1.4	-
    Lateral Friction Load Rating	350	kg
    Lateral Friction Sensitivity	0.0001	1/N
    Front Cornering Stiffness	600	N/deg
    Rear Cornering Stiffness	750	N/deg
    Power Factor Multiplier	1.25	-
    Drive Type	RWD	-
    Gear Shift Time	0.01	s
    Primary Gear Efficiency	1	-
    Final Gear Efficiency	0.92	-
    Gearbox Efficiency	0.98	-
    Primary Gear Reduction	1	-
    Final Gear Reduction	2.2	-
    1st Gear Ratio	2.9	-
    2nd Gear Ratio	2.28	-
    ...
    MARKER_info/data
    1000	150
    2000	200
    3000	260
    ...
    ```
    
    The section after MARKER_info/data is the engine torque curve:
    - Column 1: Engine RPM
    - Column 2: Torque in N·m
    
    Run: python src/CarModelSimV2.py
    Output: Model_<CarName>_CarSimV2.txt
    """)


def example_create_track_model():
    """
    Example: Create a track model from corner data.
    """
    print("\n" + "=" * 60)
    print("Example 3: Creating a Track Model")
    print("=" * 60)
    
    print("""
    Track data file format (tab-separated):
    
    ```
    France	Circuit de Monaco
    Straight	600	0
    Right	150	50
    Straight	200	0
    Left	180	75
    Right	100	30
    ...
    ```
    
    Header line: Country<tab>TrackName
    
    Each subsequent line:
    - Column 1: Corner type (Straight, Left, Right)
    - Column 2: Length/Arc length in meters
    - Column 3: Radius in meters (0 for straights)
    
    Notes:
    - Straights: Length is distance, radius is 0
    - Corners: Length is arc length (not chord), radius is turn radius
    - The sum of all lengths should equal total track length
    - Track should form a closed loop (returns to start)
    
    Run: python src/CircuitModelingSimV2.py
    Output: Model_<TrackName>_TrackSimV2.txt
    """)


def example_parameter_sensitivity():
    """
    Example: How to study parameter sensitivity.
    """
    print("\n" + "=" * 60)
    print("Example 4: Parameter Sensitivity Analysis")
    print("=" * 60)
    
    print("""
    To study how lap time varies with vehicle parameters:
    
    ```python
    import numpy as np
    
    # Example: Vary downforce coefficient
    Cl_values = np.linspace(-2.0, -4.0, 5)
    lap_times = []
    
    for Cl in Cl_values:
        # Modify vehicle parameter
        veh['Cl'] = Cl
        
        # Run simulation (you'd need to extract the simulation logic)
        lap_time = run_simulation(veh, tr)
        lap_times.append(lap_time)
    
    # Plot results
    plt.plot(Cl_values, lap_times)
    plt.xlabel('Downforce Coefficient')
    plt.ylabel('Lap Time [s]')
    plt.title('Sensitivity to Downforce')
    plt.show()
    ```
    
    Key parameters to study:
    - Downforce (Cl): Higher = faster in corners, slower on straights
    - Mass: Lower = faster everywhere
    - Power: More = faster on straights
    - Tire grip (mu): Higher = faster everywhere
    - Gear ratios: Optimization depends on track
    """)


def example_compare_cars():
    """
    Example: Compare different cars on the same track.
    """
    print("\n" + "=" * 60)
    print("Example 5: Comparing Multiple Vehicles")
    print("=" * 60)
    
    print("""
    To compare different cars on the same track:
    
    ```python
    cars = [
        "Model_Alpine_A110_CarSimV2.txt",
        "Model_Porsche_911_CarSimV2.txt",
        "Model_F1_Car_CarSimV2.txt"
    ]
    
    track = "Model_Monza_TrackSimV2.txt"
    
    results = {}
    for car_file in cars:
        car_name, veh = ouvertureFileCar(car_file)
        _, _, tr = ouvertureFileTrack(track)
        
        lap_time = run_simulation(veh, tr)
        results[car_name] = lap_time
    
    # Print comparison
    for car, time in sorted(results.items(), key=lambda x: x[1]):
        print(f"{car}: {time:.3f}s")
    ```
    
    This helps understand which car characteristics matter most
    for a given circuit (power vs. downforce vs. weight, etc.)
    """)


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("LapTimeSim - Example Usage Guide")
    print("=" * 60)
    
    example_basic_simulation()
    example_create_car_model()
    example_create_track_model()
    example_parameter_sensitivity()
    example_compare_cars()
    
    print("\n" + "=" * 60)
    print("End of Examples")
    print("=" * 60)
    print("\nFor more information, see the documentation in docs/")
