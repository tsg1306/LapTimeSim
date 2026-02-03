#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CarModelSimV2.py - Vehicle Parameter Processing Module

This module processes raw vehicle specification data and generates a comprehensive
vehicle model file containing all parameters needed for lap time simulation.

The module handles:
    - Mass and weight distribution
    - Aerodynamic coefficients (lift/downforce, drag)
    - Tire friction characteristics with load sensitivity
    - Powertrain modeling (engine curve, gearbox ratios)
    - Gear shift optimization
    - Force calculations (aerodynamic, rolling resistance, traction)
    - GGV (G-G-Velocity) diagram generation for traction envelope visualization

Author: Original development for TIPE project
Date: October 2023
Python Version: 3.7+

Dependencies:
    - numpy: Numerical computations
    - matplotlib: Visualization
    - scipy: Interpolation functions

Usage:
    python CarModelSimV2.py
    
    Modify the input file path in the main section to process different vehicles.
    Output is saved as 'Model_<CarName>_CarSimV2.txt'
"""

from math import sqrt, pi, cos, sin
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def ouvertureFileCar(fichier):
    """
    Read and parse vehicle specification file.
    
    The file format consists of:
    - Line 1: Vehicle name
    - Lines 2-N: Tab-separated parameters (name, value, unit)
    - After 'MARKER_info/data': Engine torque curve (RPM, Torque)
    
    Parameters
    ----------
    fichier : str
        Path to the vehicle specification file
        
    Returns
    -------
    tuple
        nom : str - Vehicle name
        info : list - List of parameter values (as strings)
        data : list - Engine curve data as [[RPM, Torque], ...]
        
    Example
    -------
    >>> nom, info, data = ouvertureFileCar("AlpineA110V2.txt")
    >>> print(nom)
    'Alpine A110 (cup)'
    """
    file = open(fichier, "r")
    nom = file.readline()[:-1]  # Remove newline character
    info = []
    data = []
    VU = False  # Flag to track if we've passed the data marker
    
    for ligne in file:
        if "MARKER_info/data" in ligne:
            VU = True
        elif not VU:
            # Extract parameter value (second column)
            info.append(ligne.split("\t")[1])
        elif VU and not "MARKER_info/data" in ligne:
            # Parse engine curve data
            data.append([float(ligne.split("\t")[0]), float(ligne.split("\t")[1][:-1])])
    
    return nom, info, data


def TraceCourbe(X, Y):
    """
    Simple 2D line plot utility.
    
    Parameters
    ----------
    X : array-like
        X-axis values
    Y : array-like
        Y-axis values
    """
    plt.plot(X, Y)
    plt.show()


def TraceCourbeCAR():
    """
    Generate comprehensive visualization of vehicle characteristics.
    
    Creates three plots:
    1. Gear shift diagram: Engine RPM and gear number vs vehicle speed
    2. GGV (traction circle) 3D surface plot
    3. Transmission curves: All forces vs vehicle speed
    
    Uses global variables for plotting (defined in main execution block).
    """
    # =========================================================================
    # Plot 1: Gear Shift Diagram
    # Shows engine RPM (blue) and current gear (red) against vehicle speed
    # =========================================================================
    fig, ax = plt.subplots()
    ax.plot(vehicule_speed * 3.6, engine_speed, "b")
    ax.set_xlabel("Vehicle Speed [km/h]", fontsize=14)
    ax.set_ylabel("Engine Speed [rpm]", fontsize=14)
    
    # Secondary Y-axis for gear number
    ax2 = ax.twinx()
    ax2.plot(vehicule_speed * 3.6, gear, "r")
    plt.title("Gear Shift Diagram")
    ax2.set_ylabel("Gear", fontsize=14)
    
    # Create legend for both axes
    lines = [ax.get_lines()[0], ax2.get_lines()[0]]
    plt.legend(lines, ["RPM", "Gear"])
    
    # =========================================================================
    # Plot 2: GGV (Traction Circle) 3D Surface
    # Shows the friction envelope as a function of velocity
    # X-axis: Longitudinal acceleration, Y-axis: Lateral acceleration, Z-axis: Velocity
    # =========================================================================
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(GGVx, GGVy, GGVz, cmap=cm.coolwarm, linewidth=0)
    plt.title("GGV Diagram (Traction Envelope)")
    ax.set_xlabel('Longitudinal Accel [m/s²]')
    ax.set_ylabel('Lateral Accel [m/s²]')
    ax.set_zlabel('Velocity [m/s]')
    plt.tight_layout()
    
    # =========================================================================
    # Plot 3: Transmission Force Curves
    # Shows all relevant forces as a function of vehicle speed
    # =========================================================================
    fig2, ax3 = plt.subplots()
    plt.title("Transmission & Force Curves")
    
    # Engine traction force (with power factor)
    ax3.plot(vehicule_speed * 3.6, scalList(fx_engine, factor_power))
    ax3.set_xlabel("Vehicle Speed [km/h]", fontsize=14)
    ax3.set_ylabel("Force [N]", fontsize=14)
    
    # Final traction (minimum of engine force and tire grip)
    ax3.plot(vehicule_speed * 3.6, 
             minListTAT(fx_tyre, scalList(fx_engine, factor_power)),
             color="black", linewidth=6)
    
    # Aerodynamic drag (negative because it opposes motion)
    ax3.plot(vehicule_speed * 3.6, scalList(fx_aero, -1), linewidth=3)
    
    # Rolling resistance
    ax3.plot(vehicule_speed * 3.6, scalList(fx_roll, -1))
    
    # Maximum tire traction force
    ax3.plot(vehicule_speed * 3.6, fx_tyre, linewidth=3)
    
    # Individual gear curves
    for i in range(nog):
        ax3.plot((vehicule_speed * 3.6)[1:], transpo(fx)[i])
    
    # Legend
    lines = [ax3.get_lines()[0], ax3.get_lines()[1], ax3.get_lines()[2], 
             ax3.get_lines()[3], ax3.get_lines()[4]]
    plt.legend(lines, ["Engine Traction", "Final Traction", "Aero Drag", 
                       "Rolling Resistance", "Max Tire Traction"])
    
    plt.show()


# =============================================================================
# MATRIX AND LIST OPERATIONS
# =============================================================================

def minMatrice(M):
    """
    Find minimum value in a 2D matrix.
    
    Parameters
    ----------
    M : list of lists
        2D matrix
        
    Returns
    -------
    float
        Minimum value in the matrix
    """
    res = M[0][0]
    for elt in M:
        if min(elt) < res:
            res = min(elt)
    return res


def maxMatrice(M):
    """
    Find maximum value in a 2D matrix.
    
    Parameters
    ----------
    M : list of lists
        2D matrix
        
    Returns
    -------
    float
        Maximum value in the matrix
    """
    res = M[0][0]
    for elt in M:
        if max(elt) > res:
            res = max(elt)
    return res


def indexGearMax(L):
    """
    Find index of maximum value in list (1-indexed for gear number).
    
    Parameters
    ----------
    L : list
        List of values (typically force values for each gear)
        
    Returns
    -------
    int
        Index of maximum value + 1 (gear number)
    """
    ind = 0
    m = L[0]
    for i in range(len(L)):
        if m < L[i]:
            ind = i
            m = L[i]
    return ind + 1


def transpo(M):
    """
    Transpose a 2D matrix.
    
    Parameters
    ----------
    M : list of lists
        Input matrix of shape (n, m)
        
    Returns
    -------
    list of lists
        Transposed matrix of shape (m, n)
    """
    res = [[0 for i in range(len(M))] for j in range(len(M[0]))]
    for i in range(len(M)):
        for j in range(len(M[0])):
            res[j][i] = M[i][j]
    return res


def scalList(L, scal):
    """
    Multiply each element of a list by a scalar.
    
    Parameters
    ----------
    L : list
        Input list
    scal : float
        Scalar multiplier
        
    Returns
    -------
    list
        Scaled list
    """
    res = [0 for i in range(len(L))]
    for i in range(len(L)):
        res[i] = L[i] * scal
    return res


def scalMatrice(M, scal):
    """
    Multiply each element of a matrix by a scalar.
    
    Parameters
    ----------
    M : list of lists
        Input matrix
    scal : float
        Scalar multiplier
        
    Returns
    -------
    list of lists
        Scaled matrix
    """
    res = [[0 for i in range(len(M[0]))] for j in range(len(M))]
    for i in range(len(M)):
        for j in range(len(M[0])):
            res[i][j] = M[i][j] * scal
    return res


def addList(A, B):
    """
    Element-wise addition of two lists.
    
    Parameters
    ----------
    A : list
        First list
    B : list
        Second list (must be same length as A)
        
    Returns
    -------
    list
        Element-wise sum
    """
    res = [0 for i in range(len(A))]
    for i in range(len(A)):
        res[i] = A[i] + B[i]
    return res


def addMatrice(M, B):
    """
    Element-wise addition of two matrices.
    
    Parameters
    ----------
    M : list of lists
        First matrix
    B : list of lists
        Second matrix (must be same shape as M)
        
    Returns
    -------
    list of lists
        Element-wise sum
    """
    res = [[0 for i in range(len(M[0]))] for j in range(len(M))]
    for i in range(len(M)):
        for j in range(len(M[0])):
            res[i][j] = M[i][j] + B[i][j]
    return res


def multListTAT(M, B):
    """
    Element-wise (term-by-term) multiplication of two lists.
    
    Parameters
    ----------
    M : list
        First list
    B : list
        Second list (must be same length as M)
        
    Returns
    -------
    list
        Element-wise product (Hadamard product)
    """
    res = [0 for i in range(len(M))]
    for i in range(len(M)):
        res[i] = M[i] * B[i]
    return res


def multMatriceTAT(M, B):
    """
    Element-wise (term-by-term) multiplication of two matrices.
    
    Parameters
    ----------
    M : list of lists
        First matrix
    B : list of lists
        Second matrix (must be same shape as M)
        
    Returns
    -------
    list of lists
        Element-wise product (Hadamard product)
    """
    res = [[0 for i in range(len(M[0]))] for j in range(len(M))]
    for i in range(len(M)):
        for j in range(len(M[0])):
            res[i][j] = M[i][j] * B[i][j]
    return res


def minListTAT(M, B):
    """
    Element-wise minimum of two lists.
    
    Parameters
    ----------
    M : list
        First list
    B : list
        Second list (must be same length as M)
        
    Returns
    -------
    list
        Element-wise minimum
    """
    res = [0 for i in range(len(M))]
    for i in range(len(M)):
        res[i] = min(M[i], B[i])
    return res


def numRapTORap(gear, ratio_gearbox):
    """
    Convert gear numbers to actual gear ratios.
    
    Parameters
    ----------
    gear : list
        List of gear numbers (1-indexed)
    ratio_gearbox : list
        List of gear ratios for each gear
        
    Returns
    -------
    list
        Gear ratios corresponding to each gear number
    """
    res = [0 for i in range(len(gear))]
    for i in range(len(gear)):
        res[i] = ratio_gearbox[gear[i] - 1]
    return res


def invList(L):
    """
    Compute element-wise reciprocal of a list.
    
    Parameters
    ----------
    L : list
        Input list (no zeros allowed)
        
    Returns
    -------
    list
        Element-wise reciprocal (1/x for each x)
    """
    res = [0 for i in range(len(L))]
    for i in range(len(L)):
        res[i] = 1 / L[i]
    return res


def deltaList(L):
    """
    Compute forward differences of a list.
    
    Parameters
    ----------
    L : list
        Input list
        
    Returns
    -------
    list
        List of differences: [L[1]-L[0], L[2]-L[1], ...]
    """
    res = [0 for i in range(len(L) - 1)]
    for i in range(len(L) - 1):
        res[i] = L[i + 1] - L[i]
    return res


def shiftANDrev(L):
    """
    Extract gear shift and arrival points from gear change indicator list.
    
    Identifies points where gear changes occur:
    - shift: Points where an upshift is initiated
    - rev: Points where the new gear is engaged
    
    Parameters
    ----------
    L : list
        Gear change indicator list (from deltaList of gear numbers)
        
    Returns
    -------
    tuple
        shift : list - RPM values at shift initiation
        rev : list - RPM values after shift completion
    """
    shift, rev = [], []
    for i in range(1, len(L) - 1):
        if L[i] != 0 and L[i + 1] != 0:
            shift.append(L[i])
        elif L[i] != 0 and L[i - 1] != 0:
            rev.append(L[i])
    return shift, rev


def absList(L):
    """
    Compute element-wise absolute value of a list.
    
    Parameters
    ----------
    L : list
        Input list
        
    Returns
    -------
    list
        List of absolute values
    """
    res = []
    for val in L:
        res.append(abs(val))
    return res


def fit(L, K):
    """
    Create a list of constant value with same length as input.
    
    Parameters
    ----------
    L : list
        Reference list (for length)
    K : float
        Constant value to fill
        
    Returns
    -------
    list
        List of length len(L) filled with K
    """
    return [K for i in range(len(L))]


def array2matrice(A):
    """
    Convert numpy 2D array to Python list of lists.
    
    Parameters
    ----------
    A : numpy.ndarray
        2D numpy array
        
    Returns
    -------
    list of lists
        Equivalent Python nested list
    """
    res = [[0 for i in range(len(A[0]))] for j in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A[0])):
            res[i][j] = A[i][j]
    return res


def array2list(A):
    """
    Convert numpy 1D array to Python list.
    
    Parameters
    ----------
    A : numpy.ndarray
        1D numpy array
        
    Returns
    -------
    list
        Equivalent Python list
    """
    res = [0 for i in range(len(A))]
    for i in range(len(A)):
        res[i] = A[i]
    return res


def array2cube(A):
    """
    Convert numpy 3D array to Python nested list (3 levels).
    
    Parameters
    ----------
    A : numpy.ndarray
        3D numpy array
        
    Returns
    -------
    list of lists of lists
        Equivalent Python nested list
    """
    res = [[[0 for i in range(len(A[0][0]))] for j in range(len(A[0]))] for p in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A[0])):
            for p in range(len(A[0][0])):
                res[i][j][p] = A[i][j][p]
    return res


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    
    # =========================================================================
    # SECTION 1: READ VEHICLE DATA
    # =========================================================================
    print("=" * 60)
    print("LapTimeSim - Vehicle Model Generator")
    print("=" * 60)
    
    # Read vehicle specification file
    # Modify this path to process different vehicles
    nom, info, data = ouvertureFileCar("../data/cars/essaieCARV2.txt")
    
    print(f"\nProcessing vehicle: {nom}")
    
    # =========================================================================
    # SECTION 2: PARSE VEHICLE PARAMETERS
    # =========================================================================
    
    # Initialize export dictionary (will contain all computed parameters)
    i = 0  # Parameter index
    export = {}
    
    # -------------------------------------------------------------------------
    # Mass Properties
    # -------------------------------------------------------------------------
    M = float(info[i])              # Total mass [kg]
    export["M"] = M
    i += 1
    
    df = float(info[i]) / 100       # Front mass distribution [fraction]
    export["df"] = df               # e.g., 0.42 means 42% front, 58% rear
    i += 1
    
    # -------------------------------------------------------------------------
    # Wheelbase
    # -------------------------------------------------------------------------
    L = float(info[i]) / 1000       # Wheelbase [m] (converted from mm)
    export["L"] = L
    i += 1
    
    # -------------------------------------------------------------------------
    # Aerodynamic Parameters
    # -------------------------------------------------------------------------
    Cl = float(info[i])             # Lift coefficient (negative = downforce)
    export["Cl"] = Cl
    i += 1
    
    Cd = float(info[i])             # Drag coefficient (always negative)
    export["Cd"] = Cd
    i += 1
    
    factor_Cl = float(info[i])      # Lift coefficient multiplier (for tuning)
    export["factor_Cl"] = factor_Cl
    i += 1
    
    factor_Cd = float(info[i])      # Drag coefficient multiplier (for tuning)
    export["factor_Cd"] = factor_Cd
    i += 1
    
    da = float(info[i]) / 100       # Front aerodynamic distribution [fraction]
    export["da"] = da               # Affects front/rear downforce split
    i += 1
    
    A = float(info[i])              # Frontal area [m²]
    export["A"] = A
    i += 1
    
    rho = float(info[i])            # Air density [kg/m³] (typically 1.225)
    export["rho"] = rho
    i += 1
    
    # -------------------------------------------------------------------------
    # Tire Parameters
    # -------------------------------------------------------------------------
    factor_grip = float(info[i])    # Overall grip multiplier
    export["factor_grip"] = factor_grip
    i += 1
    
    tyre_radius = float(info[i]) / 1000  # Tire radius [m] (converted from mm)
    export["tyre_radius"] = tyre_radius
    i += 1
    
    Cr = float(info[i])             # Rolling resistance coefficient
    export["Cr"] = Cr
    i += 1
    
    # Longitudinal (braking/acceleration) tire characteristics
    mu_x = float(info[i])           # Base longitudinal friction coefficient
    export["mu_x"] = mu_x
    i += 1
    
    mu_x_M = float(info[i])         # Rated load for longitudinal friction [kg]
    export["mu_x_M"] = mu_x_M
    i += 1
    
    sens_x = float(info[i])         # Load sensitivity [1/N]
    export["sens_x"] = sens_x       # Friction decreases as load increases
    i += 1
    
    # Lateral (cornering) tire characteristics
    mu_y = float(info[i])           # Base lateral friction coefficient
    export["mu_y"] = mu_y
    i += 1
    
    mu_y_M = float(info[i])         # Rated load for lateral friction [kg]
    export["mu_y_M"] = mu_y_M
    i += 1
    
    sens_y = float(info[i])         # Lateral load sensitivity [1/N]
    export["sens_y"] = sens_y
    i += 1
    
    # Cornering stiffness (for steering response, not used in basic sim)
    CF = float(info[i])             # Front cornering stiffness [N/deg]
    export["CF"] = CF
    i += 1
    
    CR = float(info[i])             # Rear cornering stiffness [N/deg]
    export["CR"] = CR
    i += 1
    
    # -------------------------------------------------------------------------
    # Powertrain Parameters
    # -------------------------------------------------------------------------
    factor_power = float(info[i])   # Power multiplier (for tuning)
    export["factor_power"] = factor_power
    i += 1
    
    # Transmission configuration
    drive = info[i]                 # Drive type: 'RWD', 'FWD', or 'AWD'
    export["drive"] = drive
    i += 1
    
    shift_time = float(info[i])     # Gear shift time [s]
    export["shift_time"] = shift_time
    i += 1
    
    # Drivetrain efficiencies
    n_primary = float(info[i])      # Primary gear efficiency
    export["n_primary"] = n_primary
    i += 1
    
    n_final = float(info[i])        # Final drive efficiency
    export["n_final"] = n_final
    i += 1
    
    n_gearbox = float(info[i])      # Gearbox efficiency
    export["n_gearbox"] = n_gearbox
    i += 1
    
    # Gear ratios
    ratio_primary = float(info[i])  # Primary gear ratio
    export["ratio_primary"] = ratio_primary
    i += 1
    
    ratio_final = float(info[i])    # Final drive ratio
    export["ratio_final"] = ratio_final
    i += 1
    
    # Individual gear ratios
    ratio_gearbox = []
    while i < len(info):
        ratio_gearbox.append(float(info[i]))
        i += 1
    export["ratio_gearbox"] = ratio_gearbox
    
    nog = len(ratio_gearbox)        # Number of gears
    export["nog"] = nog
    
    print("✓ Parameter parsing complete")
    
    # =========================================================================
    # SECTION 3: TRANSMISSION MODEL
    # =========================================================================
    # This section converts engine torque curve to wheel force vs vehicle speed
    # for each gear, then determines optimal gear selection
    
    # Parse engine torque curve
    en_speed_curve = []             # Engine speed [rpm]
    en_torque_curve = []            # Engine torque [N·m]
    en_power_curve = []             # Engine power [W]
    
    for ligne in data:
        en_speed_curve.append(ligne[0])
        en_torque_curve.append(ligne[1])
        # Power = Torque × Angular velocity = Torque × RPM × (2π/60)
        en_power_curve.append(ligne[0] * ligne[1] * (2 * pi / 60))
    
    # Initialize arrays for each gear
    # Dimensions: [engine_speed_points, number_of_gears]
    wheel_speed_gear = [[0 for i in range(nog)] for j in range(len(en_speed_curve))]
    vehicule_speed_gear = [[0 for i in range(nog)] for j in range(len(en_speed_curve))]
    wheel_torque_gear = [[0 for i in range(nog)] for j in range(len(en_speed_curve))]
    
    # Calculate wheel speed, vehicle speed, and wheel torque for each gear
    for j in range(nog):
        for i in range(len(en_speed_curve)):
            # Wheel speed = Engine speed / (total gear ratio)
            wheel_speed_gear[i][j] = en_speed_curve[i] / (ratio_final * ratio_primary * ratio_gearbox[j])
            
            # Vehicle speed = Wheel angular velocity × tire radius
            # Wheel angular velocity = wheel_rpm × 2π / 60
            vehicule_speed_gear[i][j] = (wheel_speed_gear[i][j] * 2 * pi * tyre_radius) / 60
            
            # Wheel torque = Engine torque × gear ratio × efficiency
            wheel_torque_gear[i][j] = (en_torque_curve[i] * ratio_final * ratio_primary * 
                                        ratio_gearbox[j] * n_final * n_primary * n_gearbox)
    
    # Store computed values
    export["vehicule_speed_gear"] = vehicule_speed_gear
    export["wheel_speed_gear"] = wheel_speed_gear
    export["wheel_torque_gear"] = wheel_torque_gear
    
    # Determine min and max vehicle speeds
    v_min = minMatrice(vehicule_speed_gear)
    export["v_min"] = v_min
    v_max = maxMatrice(vehicule_speed_gear)
    export["v_max"] = v_max
    
    # Create uniform velocity discretization
    dv = 0.5 / 3.6                  # Velocity step [m/s] (0.5 km/h)
    vehicule_speed = np.linspace(v_min, v_max, int((v_max - v_min) / dv))
    
    # Initialize arrays for optimal gear selection
    gear = [0 for i in range(len(vehicule_speed))]           # Optimal gear number
    fx_engine = [0 for i in range(len(vehicule_speed))]      # Engine traction force
    fx = [[0 for j in range(nog)] for i in range(len(vehicule_speed))]  # Force per gear
    
    # Optimize gear selection at each velocity
    # Strategy: Choose gear that maximizes traction force at each speed
    for i in range(len(vehicule_speed)):
        for j in range(nog):
            # Interpolate wheel torque at current speed for this gear
            temp = interp1d(transpo(vehicule_speed_gear)[j],
                           scalList(transpo(wheel_torque_gear)[j], 1 / tyre_radius),
                           'linear')
            
            # Check if current speed is within gear's operating range
            if (vehicule_speed[i] < transpo(vehicule_speed_gear)[j][0] or 
                vehicule_speed[i] > transpo(vehicule_speed_gear)[j][-1]):
                fx[i][j] = 0
            else:
                fx[i][j] = round(float(np.array2string(temp(vehicule_speed[i]))), 2)
        
        # Select gear with maximum traction force
        fx_engine[i] = max(fx[i])
        gear[i] = indexGearMax(fx[i])
    
    # Add initial point at v=0 for interpolation
    vehicule_speed = np.append(np.array([0]), vehicule_speed)
    export["vehicule_speed"] = array2list(vehicule_speed)
    
    gear = [gear[0]] + gear
    export["gear"] = gear
    
    fx_engine = [fx_engine[0]] + fx_engine
    export["fx_engine"] = fx_engine
    
    # Calculate engine speed profile
    # Engine RPM = Vehicle speed × gear ratio × final ratio × primary ratio / (tire circumference)
    engine_speed = scalList(
        multListTAT(numRapTORap(gear, ratio_gearbox), vehicule_speed),
        (ratio_final * ratio_primary * 60) / (tyre_radius * 2 * pi)
    )
    export["engine_speed"] = engine_speed
    
    # Calculate wheel torque profile
    wheel_torque = scalList(fx_engine, tyre_radius)
    export["wheel_torque"] = wheel_torque
    
    # Calculate engine torque profile
    engine_torque = scalList(
        multListTAT(wheel_torque, invList(numRapTORap(gear, ratio_gearbox))),
        1 / (ratio_final * ratio_primary * n_final * n_primary * n_gearbox)
    )
    export["engine_torque"] = engine_torque
    
    # Calculate engine power profile
    engine_power = scalList(multListTAT(engine_torque, engine_speed), 2 * pi / 60)
    export["engine_power"] = engine_power
    
    print("✓ Transmission model complete")
    
    # =========================================================================
    # SECTION 4: GEAR SHIFT MODEL
    # =========================================================================
    # Identify gear change points and calculate RPM drops
    
    # Detect gear changes (delta != 0)
    gear_change = addList(deltaList(gear) + [0], [0] + deltaList(gear))
    export["gear_change"] = gear_change
    
    # Calculate engine RPM at gear changes
    engine_speed_gear_change = multListTAT(gear_change, engine_speed)
    export["engine_speed_gear_change"] = engine_speed_gear_change
    
    # Identify shift initiation and completion points
    shift_points, arrive_points = shiftANDrev(engine_speed_gear_change)
    export["shift_points"] = shift_points
    export["arrive_points"] = arrive_points
    
    # Calculate RPM drop during each shift
    rev_drops = addList(shift_points, scalList(arrive_points, -1))
    export["rev_drops"] = rev_drops
    
    print("✓ Gear shift model complete")
    
    # =========================================================================
    # SECTION 5: FORCE MODEL
    # =========================================================================
    # Calculate all forces acting on the vehicle as functions of velocity
    
    g = 9.81  # Gravitational acceleration [m/s²]
    
    # Drive type factors
    # These determine which wheels receive power and how weight is distributed
    if drive == 'RWD':  # Rear-wheel drive
        factor_drive = (1 - df)     # Fraction of weight on driven wheels
        factor_aero = (1 - da)      # Fraction of aero load on driven wheels
        driven_wheels = 2
    elif drive == 'FWD':  # Front-wheel drive
        factor_drive = df
        factor_aero = da
        driven_wheels = 2
    elif drive == 'AWD':  # All-wheel drive
        factor_drive = 1
        factor_aero = 1
        driven_wheels = 4
    
    export["factor_drive"] = factor_drive
    export["factor_aero"] = factor_aero
    export["driven_wheels"] = driven_wheels
    
    # -------------------------------------------------------------------------
    # Vertical (Z-axis) Forces
    # -------------------------------------------------------------------------
    
    # Weight force (constant, pointing down)
    fz_mass = -M * g
    export["fz_mass"] = fz_mass
    
    # Aerodynamic downforce (varies with v²)
    # F_aero = 0.5 × ρ × Cl × A × v²
    fz_aero = scalList(multListTAT(vehicule_speed, vehicule_speed),
                       0.5 * rho * factor_Cl * Cl * A)
    export["fz_aero"] = fz_aero
    
    # Total vertical force
    fz_total = addList(fz_aero, fit(fz_aero, fz_mass))
    export["fz_total"] = fz_total
    
    # Normal force on driven wheels
    fz_tyre = scalList(
        addList(fit(fz_aero, fz_mass * factor_drive),
                scalList(fz_aero, factor_aero)),
        1 / driven_wheels
    )
    export["fz_tyre"] = fz_tyre
    
    # -------------------------------------------------------------------------
    # Longitudinal (X-axis) Forces
    # -------------------------------------------------------------------------
    
    # Aerodynamic drag (varies with v²)
    fx_aero = scalList(multListTAT(vehicule_speed, vehicule_speed),
                       0.5 * rho * factor_Cd * Cd * A)
    export["fx_aero"] = fx_aero
    
    # Rolling resistance (proportional to normal force)
    fx_roll = scalList(absList(fz_total), Cr)
    export["fx_roll"] = fx_roll
    
    # Maximum tire traction force (with load sensitivity)
    # μ = μ_0 + sensitivity × (N_rated - N_actual)
    fx_tyre = scalList(
        multListTAT(
            absList(fz_tyre),
            addList(
                fit(fz_tyre, mu_x),
                scalList(
                    addList(fit(fz_tyre, mu_x_M * g),
                           scalList(absList(fz_tyre), -1)),
                    sens_x
                )
            )
        ),
        driven_wheels
    )
    export["fx_tyre"] = fx_tyre
    
    print("✓ Force model complete")
    
    # =========================================================================
    # SECTION 6: GGV DIAGRAM (TRACTION ENVELOPE)
    # =========================================================================
    # Generate the 3D traction circle showing available acceleration
    # as a function of velocity
    
    # Lateral tire coefficients
    dmy = factor_grip * sens_y      # Load sensitivity (lateral)
    export["dmy"] = dmy
    muy = factor_grip * mu_y        # Base friction (lateral)
    export["muy"] = muy
    Ny = mu_y_M * g                 # Rated normal force (lateral)
    export["Ny"] = Ny
    
    # Longitudinal tire coefficients
    dmx = factor_grip * sens_x      # Load sensitivity (longitudinal)
    export["dmx"] = dmx
    mux = factor_grip * mu_x        # Base friction (longitudinal)
    export["mux"] = mux
    Nx = mu_x_M * g                 # Rated normal force (longitudinal)
    export["Nx"] = Nx
    
    # Total weight force on 4 tires
    Wz = M * g
    export["Wz"] = Wz
    
    # Velocity vector for GGV
    dv = 2  # Velocity step [m/s]
    v = np.arange(0, v_max, dv)
    if v[-1] != v_max:
        v = np.append(v, v_max)
    
    # GGV array: [velocity, angle, [ax, ay, v]]
    N = 45  # Number of angular divisions (0° to 180°)
    GGV = np.zeros((len(v), 2 * N - 1, 3))
    GGVx = np.zeros((len(v), 2 * N - 1))  # Longitudinal acceleration
    GGVy = np.zeros((len(v), 2 * N - 1))  # Lateral acceleration
    GGVz = np.zeros((len(v), 2 * N - 1))  # Velocity
    
    for i in range(len(v)):
        # Aerodynamic forces at current velocity
        Aero_Df = 0.5 * rho * factor_Cl * Cl * A * v[i]**2  # Downforce
        Aero_Dr = 0.5 * rho * factor_Cd * Cd * A * v[i]**2  # Drag
        
        # Rolling resistance
        Roll_Dr = Cr * abs(-Aero_Df + Wz)
        
        # Normal force on driven wheels
        Wd = (factor_drive * Wz + (-factor_aero * Aero_Df)) / driven_wheels
        
        # Acceleration due to drag (always negative/retarding)
        ax_drag = (Aero_Dr + Roll_Dr) / M
        
        # Maximum lateral acceleration
        # Uses load-sensitive friction model
        ay_max = (1 / M) * (muy + dmy * (Ny - (Wz - Aero_Df) / 4)) * (Wz - Aero_Df)
        
        # Maximum longitudinal acceleration (traction-limited)
        ax_tyre_max_acc = (1 / M) * (mux + dmx * (Nx - Wd)) * Wd * driven_wheels
        
        # Maximum deceleration (braking, all 4 wheels)
        ax_tyre_max_dec = -(1 / M) * (mux + dmx * (Nx - (Wz - Aero_Df) / 4)) * (Wz - Aero_Df)
        
        # Power-limited acceleration
        temp = interp1d(vehicule_speed, scalList(fx_engine, factor_power))
        ax_power_limit = (temp(v[i]) / M) * np.ones((N, 1))
        
        # Generate friction ellipse
        # ay varies from -ay_max to +ay_max
        ay = ay_max * np.cos(np.linspace(0, 180, N) * pi / 180)
        
        # Acceleration half of ellipse (combined grip)
        ax_tyre_acc = ax_tyre_max_acc * np.sqrt(1 - (ay / ay_max)**2)
        ax_acc = (np.minimum(ax_tyre_acc, ax_power_limit) + ax_drag)[0]
        
        # Deceleration half of ellipse (combined grip)
        ax_dec = ax_tyre_max_dec * np.sqrt(1 - (ay / ay_max)**2) + ax_drag
        
        # Combine acceleration and deceleration curves
        D1 = np.append(ax_acc, ax_dec[1:])  # Full range of ax
        D2 = np.append(ay, np.flipud(ay[1:]))  # Full range of ay (symmetric)
        D3 = [v[i] for k in range(2 * N - 1)]  # Constant velocity
        
        # Store in GGV arrays
        for j in range(2 * N - 1):
            GGV[i][j][0], GGVx[i][j] = D1[j], D1[j]
            GGV[i][j][1], GGVy[i][j] = D2[j], D2[j]
            GGV[i][j][2], GGVz[i][j] = D3[j], D3[j]
    
    # Store GGV data (converted to native Python lists for serialization)
    export["GGV"] = array2cube(GGV)
    export["GGVx"] = array2matrice(GGVx)
    export["GGVy"] = array2matrice(GGVy)
    export["GGVz"] = array2matrice(GGVz)
    
    print("✓ GGV diagram complete")
    
    # =========================================================================
    # SECTION 7: VISUALIZATION
    # =========================================================================
    TraceCourbeCAR()
    print("✓ Visualization complete")
    
    # =========================================================================
    # SECTION 8: SAVE MODEL
    # =========================================================================
    print("\n" + "=" * 60)
    print("Vehicle Modeling Complete!")
    print("=" * 60)
    
    # Save model file
    output_filename = f"Model_{nom}_CarSimV2.txt"
    save = open(output_filename, 'w')
    save.write(nom + '\n')
    save.write(str(export))
    save.close()
    
    print(f"\n✓ Model exported as: {output_filename}")
    print(f"\nVehicle Summary:")
    print(f"  Name: {nom}")
    print(f"  Mass: {M} kg")
    print(f"  Power: ~{max(en_power_curve)/1000:.0f} kW")
    print(f"  Top Speed: {v_max*3.6:.1f} km/h")
    print(f"  Gears: {nog}")
    print(f"  Drive: {drive}")
