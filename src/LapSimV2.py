#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LapSimV2.py - Quasi-Static Lap Time Simulation Engine

This is the core simulation module that calculates optimal lap times using a
quasi-static approach. It combines vehicle and track models to determine the
theoretical minimum lap time.

Algorithm Overview (Backward-Forward Method):
    1. Calculate maximum cornering speed at each point (lateral grip limit)
    2. Identify apex points (local velocity minima - braking zones)
    3. From each apex, integrate forward (acceleration) and backward (braking)
    4. At each point, take the minimum velocity from all integration passes

Key Physics:
    - Tire friction circle (combined grip model)
    - Load-sensitive friction coefficients
    - Aerodynamic downforce and drag
    - Engine power curve with gear selection
    - Rolling resistance

Author: Original development for TIPE project
Date: October 2023
Python Version: 3.7+
"""

from math import sqrt, pi, cos, sin, inf, atan
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib import cm
import ast
from scipy.signal import find_peaks


# =============================================================================
# VISUALIZATION
# =============================================================================

def TraceCourbeSIM():
    """
    Generate comprehensive visualization of simulation results.
    
    Creates four plots showing speed profile, engine/gear data,
    yaw rate, and a track map colored by velocity.
    """
    # Speed vs Time
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(temps, V * 3.6, "black", linewidth=1)
    ax.set_xlabel("Time [s]", fontsize=14)
    ax.set_ylabel("Speed [km/h]", fontsize=14)
    ax.set_title("Speed Profile")
    ax.grid(True, alpha=0.3)
    ax.fill_between(temps, 0, (V * 3.6).flatten(), alpha=0.3)
    
    # Engine RPM and Gear
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(temps, engine_speed, "b", linewidth=1)
    ax.set_xlabel("Time [s]", fontsize=14)
    ax.set_ylabel("Engine Speed [rpm]", fontsize=14, color='b')
    ax2 = ax.twinx()
    ax2.plot(temps, gear, "r", linewidth=1)
    ax2.set_ylabel("Gear", fontsize=14, color='r')
    plt.title("Engine RPM and Gear Selection")
    ax.grid(True, alpha=0.3)
    
    # Yaw Rate
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(temps, yaw_rate * 180 / pi, "black", linewidth=1)
    ax.set_xlabel("Time [s]", fontsize=14)
    ax.set_ylabel("Yaw Rate [deg/s]", fontsize=14)
    ax.set_title("Vehicle Yaw Rate")
    ax.grid(True, alpha=0.3)
    
    # Track Map with Velocity Coloring
    fig, ax = plt.subplots(figsize=(12, 10))
    scatter = plt.scatter(x=tr['X'], y=tr['Y'],
                         c=np.append(V, V[-1]) * 3.6, cmap='turbo', s=8)
    plt.colorbar(scatter, label="Speed [km/h]")
    plt.plot((tr['X'][-1] + tr['X'][0]) / 2, (tr['Y'][-1] + tr['Y'][0]) / 2,
             marker="<", color="black", markersize=15)
    plt.title("Circuit Map - Speed Distribution")
    lap_time_str = f"{int(temps[-1]//60)}:{temps[-1]%60:06.3f}"
    plt.annotate(f'Lap Time: {lap_time_str}', xycoords='figure fraction',
                 xy=(0.02, 0.02), fontsize=12, backgroundcolor='white')
    ax.set_xlabel("X [m]", fontsize=14)
    ax.set_ylabel("Y [m]", fontsize=14)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def transpo(M):
    """Transpose a 2D matrix."""
    res = [[0 for i in range(len(M))] for j in range(len(M[0]))]
    for i in range(len(M)):
        for j in range(len(M[0])):
            res[j][i] = M[i][j]
    return res


def listeInd(L, ind):
    """Extract elements at specified indices."""
    return [L[i] for i in ind]


def triParRapport(liste, liste2):
    """Sort two lists together by first list's values."""
    L, L2 = liste, liste2
    N = len(L)
    for n in range(1, N):
        cle, cle2 = L[n], L2[n]
        j = n - 1
        while j >= 0 and L[j] > cle:
            L[j + 1], L2[j + 1] = L[j], L2[j]
            j -= 1
        L[j + 1], L2[j + 1] = cle, cle2
    return L, L2


def ouvertureFileTrack(fichier):
    """Read pre-processed track model file."""
    file = open(fichier, "r")
    head = (file.readline()).split("\t")
    tr = ast.literal_eval(file.readline())
    return head[0], head[1][:-1], tr


def ouvertureFileCar(fichier):
    """Read pre-processed vehicle model file."""
    file = open(fichier, "r")
    head = file.readline()
    veh = ast.literal_eval(file.readline())
    return head[:-1], veh


def array2list(A):
    """Convert numpy array to Python list."""
    return [float(A[i]) for i in range(len(A))]


def scalList(L, scal):
    """Multiply list by scalar."""
    return [L[i] * scal for i in range(len(L))]


def fast_nearest_interp(xi, x, y):
    """
    Fast nearest-neighbor interpolation for monotonically increasing x.
    
    Parameters
    ----------
    xi : float
        Query point
    x : array
        Known x values (must be monotonically increasing)
    y : array
        Known y values
        
    Returns
    -------
    float
        Interpolated y value
    """
    spacing = np.diff(x) / 2
    x = x + np.hstack([spacing, spacing[-1]])
    y = np.hstack([y, y[-1]])
    return y[np.searchsorted(x, xi)]


def other_points(i, i_max):
    """Get list of all indices except i."""
    res = [j for j in range(i_max + 1)]
    return res[:i] + res[i + 1:]


def minDroiteGauche(A):
    """
    Find minimum value, its index, and mode (acceleration/deceleration).
    
    Parameters
    ----------
    A : array
        2D array where A[i][k] contains velocity for apex i, mode k
        
    Returns
    -------
    tuple
        (min_value, apex_index, mode)
        mode: 0 = acceleration, 1 = deceleration
    """
    mode, ind, m = -1, -1, inf
    for i in range(len(A)):
        for k in range(2):
            if A[i][k] < m:
                m, mode, ind = A[i][k], k, i
    return float(m), ind, mode


def next_point(j, j_max, mode):
    """
    Calculate next track point index for simulation iteration.
    
    Handles wrap-around at track ends for both acceleration (forward)
    and deceleration (backward) integration.
    
    Parameters
    ----------
    j : int
        Current point index
    j_max : int
        Maximum point index
    mode : int
        1 = acceleration (forward), -1 = deceleration (backward)
        
    Returns
    -------
    tuple
        (next_index, current_index_updated)
    """
    if mode == 1:  # Acceleration (forward integration)
        if j == j_max - 1:
            j, j_next = j_max, 1
        elif j == j_max:
            j, j_next = 1, 2
        else:
            j, j_next = j + 1, j + 2
    else:  # Deceleration (backward integration)
        if j == 2:
            j, j_next = 1, j_max
        elif j == 1:
            j, j_next = j_max, j_max - 1
        else:
            j, j_next = j - 1, j - 2
    return j_next, j


# =============================================================================
# PHYSICS MODELS
# =============================================================================

def vehicle_model_lat(veh, tr, p):
    """
    Calculate maximum cornering velocity based on lateral grip.
    
    Solves for the maximum speed at which the car can negotiate a corner
    considering tire grip, aerodynamic downforce, and drag forces.
    
    Physics:
        Required centripetal force: F_c = m × v² × κ (where κ = 1/R)
        Available lateral force: F_y = μ(N) × N
        Normal force: N = m×g - F_downforce(v²)
    
    Parameters
    ----------
    veh : dict
        Vehicle model parameters
    tr : dict
        Track model parameters  
    p : int
        Track point index
        
    Returns
    -------
    float or False
        Maximum velocity [m/s] or False if no solution
    """
    g = 9.81
    r = tr["tr"][p][3]  # Curvature (1/R)
    factor_grip = veh["factor_grip"]
    
    # Vehicle parameters
    factor_drive = veh["factor_drive"]
    M = veh["M"]
    Wz = M * g
    
    # Straight: limited by top speed only
    if r == 0:
        return veh["v_max"]
    
    # Solve quadratic for corner speed
    # Equation derived from: m×v²×κ = μ(v)×N(v)
    D = -0.5 * veh["rho"] * veh["factor_Cl"] * veh["Cl"] * veh["A"]
    
    dmy = factor_grip * veh["sens_y"]
    muy = factor_drive * veh["mu_y"]
    Ny = veh["mu_y_M"] * g
    
    sign = int(r / abs(r))
    a = -sign * dmy / 4 * D**2
    b = sign * (muy * D + dmy * Ny * D - 2 * (dmy / 4) * Wz * D) - M * r
    c = sign * (muy * Wz + dmy * Ny * Wz - (dmy / 4) * Wz**2)
    
    # Solve quadratic
    if a == 0:
        v = sqrt(-c / b)
    elif b**2 - 4 * a * c >= 0:
        if (-b + sqrt(b**2 - 4 * a * c)) / (2 * a) >= 0:
            v = sqrt((-b + sqrt(b**2 - 4 * a * c)) / (2 * a))
        elif (-b - sqrt(b**2 - 4 * a * c)) / (2 * a) >= 0:
            v = sqrt((-b - sqrt(b**2 - 4 * a * c)) / (2 * a))
        else:
            return False
    else:
        return False
    
    v = min(v, veh["v_max"])
    
    # Iterative adjustment for drag compensation
    ajust = True
    while ajust:
        Aero_Df = 0.5 * veh["rho"] * veh['factor_Cl'] * veh['Cl'] * veh['A'] * v**2
        Aero_Dr = 0.5 * veh["rho"] * veh['factor_Cd'] * veh['Cd'] * veh['A'] * v**2
        Roll_Dr = veh["Cr"] * (-Aero_Df + Wz)
        
        Wd = (factor_drive * Wz - veh["factor_aero"] * Aero_Df) / veh["driven_wheels"]
        ax_drag = (Aero_Dr + Roll_Dr) / M
        
        ay_max = sign / M * (muy + dmy * (Ny - (Wz - Aero_Df) / 4)) * (Wz - Aero_Df)
        ay_needed = v**2 * r
        
        dmx = factor_grip * veh["sens_x"]
        mux = factor_drive * veh["mu_x"]
        Nx = veh["mu_x_M"] * g
        
        if ax_drag <= 0:  # Need throttle to maintain speed
            ax_tyre_max_acc = 1 / M * (mux + dmx * (Nx - Wd)) * Wd * veh["driven_wheels"]
            temp = interp1d(veh["vehicule_speed"], veh["factor_power"] * np.array(veh["fx_engine"]))
            ay = ay_max * sqrt(1 - (ax_drag / ax_tyre_max_acc)**2)
        else:  # Need braking
            ax_tyre_max_dec = -1 / M * (mux + dmx * (Nx - (Wz - Aero_Df) / 4)) * (Wz - Aero_Df)
            ay = ay_max * sqrt(1 - (ax_drag / ax_tyre_max_dec)**2)
        
        if ay / ay_needed < 1:  # Not enough grip
            v = sqrt(ay / r) - 1e-3
        else:
            ajust = False
    
    return v


def vehicle_model_comb(veh, tr, v, v_max_next, j, mode):
    """
    Calculate next velocity using combined (longitudinal + lateral) grip model.
    
    This is the core integration step that accounts for the friction ellipse:
    when cornering, available longitudinal acceleration is reduced.
    
    The friction ellipse equation:
        (a_x / a_x_max)² + (a_y / a_y_max)² = 1
    
    Parameters
    ----------
    veh : dict
        Vehicle model parameters
    tr : dict
        Track model parameters
    v : float
        Current velocity [m/s]
    v_max_next : float
        Maximum allowed velocity at next point [m/s]
    j : int
        Current track point index
    mode : int
        1 = acceleration, -1 = deceleration
        
    Returns
    -------
    tuple
        (v_next, overshoot)
        v_next : float - Velocity at next point [m/s]
        overshoot : bool - True if velocity limit exceeded
    """
    overshoot = False
    dx = tr['tr'][j][1]  # Segment length
    r = tr['tr'][j][3]   # Curvature
    factor_grip = veh['factor_grip']
    g = 9.81
    
    # Mode-specific parameters
    if mode == 1:  # Acceleration
        factor_drive = veh['factor_drive']
        factor_aero = veh['factor_aero']
        driven_wheels = veh['driven_wheels']
    else:  # Braking (all 4 wheels)
        factor_drive = 1
        factor_aero = 1
        driven_wheels = 4
    
    M = veh['M']
    Wz = M * g
    
    # Aerodynamic forces
    Aero_Df = 0.5 * veh['rho'] * veh['factor_Cl'] * veh['Cl'] * veh['A'] * v**2
    Aero_Dr = 0.5 * veh['rho'] * veh['factor_Cd'] * veh['Cd'] * veh['A'] * v**2
    Roll_Dr = veh['Cr'] * (-Aero_Df + Wz)
    Wd = (factor_drive * Wz - factor_aero * Aero_Df) / driven_wheels
    
    # Required acceleration to reach v_max_next
    ax_max = mode * (v_max_next**2 - v**2) / (2 * dx)
    ax_drag = (Aero_Dr + Roll_Dr) / M
    ax_needed = ax_max - ax_drag
    
    # Current lateral acceleration
    ay = v**2 * r
    
    # Tire friction coefficients
    dmy = factor_grip * veh['sens_y']
    muy = factor_grip * veh['mu_y']
    Ny = veh['mu_y_M'] * g
    dmx = factor_grip * veh['sens_x']
    mux = factor_grip * veh['mu_x']
    Nx = veh['mu_x_M'] * g
    
    # Friction ellipse calculation
    if ay != 0:  # In a corner
        ay_max = 1 / M * (int(ay / abs(ay)) * (muy + dmy * (Ny - (Wz - Aero_Df) / 4)) * (Wz - Aero_Df))
        ellipse_multi = 0 if abs(ay / ay_max) > 1 else sqrt(1 - (ay / ay_max)**2)
    else:  # Straight
        ellipse_multi = 1
    
    # Calculate achievable acceleration
    if ax_needed >= 0:  # Need acceleration
        ax_tyre_max = 1 / M * (mux + dmx * (Nx - Wd)) * Wd * driven_wheels
        ax_tyre = ax_tyre_max * ellipse_multi
        
        temp = interp1d(veh['vehicule_speed'], veh['factor_power'] * np.array(veh['fx_engine']))
        if v < veh['vehicule_speed'][0] or v > veh['vehicule_speed'][-1]:
            ax_power_limit = 0
        else:
            ax_power_limit = 1 / M * temp(v)
        
        scale = min(ax_tyre, ax_needed) / ax_power_limit if ax_power_limit > 0 else 0
        tps = max(min(1, scale), 0)
        ax_com = tps * ax_power_limit
    else:  # Need braking
        ax_tyre_max = -1 / M * (mux + dmx * (Nx - (Wz - Aero_Df) / 4)) * (Wz - Aero_Df)
        ax_tyre = ax_tyre_max * ellipse_multi
        tps = 0
        ax_com = -min(-ax_tyre, -ax_needed)
    
    # Calculate next velocity: v² = v₀² + 2×a×dx
    ax = ax_com + ax_drag
    v_next = sqrt(v**2 + 2 * mode * ax * dx)
    
    # Full throttle correction at top speed
    if tps > 0 and v / veh['v_max'] >= 0.999:
        tps = 1
    
    # Check for overshoot
    if v_next / v_max_next > 1:
        if v > v_max_next:
            overshoot = True
            v_next = inf
            return v_next, overshoot
        else:
            v_next = (v + v_max_next) / 2
    
    return v_next, overshoot


# =============================================================================
# MAIN SIMULATION
# =============================================================================

if __name__ == "__main__":
    
    print("=" * 60)
    print("LapTimeSim - Quasi-Static Lap Time Simulator")
    print("=" * 60)
    
    # =========================================================================
    # LOAD MODELS
    # =========================================================================
    
    # Load track model (modify path as needed)
    nom_circuit, pays, tr = ouvertureFileTrack("../data/tracks/Model_Monza_TrackSimV2.txt")
    print(f"\nTrack: {nom_circuit} ({pays})")
    print(f"  Length: {tr['Lcircuit']} m")
    
    # Load vehicle model (modify path as needed)
    nom_voiture, veh = ouvertureFileCar("../data/cars/Model_Alpine_A110_CarSimV2.txt")
    print(f"\nVehicle: {nom_voiture}")
    print(f"  Mass: {veh['M']} kg")
    print(f"  Top Speed: {veh['v_max']*3.6:.1f} km/h")
    
    # =========================================================================
    # PHASE 1: CALCULATE MAXIMUM CORNERING SPEEDS
    # =========================================================================
    print("\n[1/3] Calculating maximum cornering speeds...")
    
    # Calculate lateral-limited velocity at each point
    v_max = np.zeros((len(tr["tr"]), 1))
    for i in range(len(tr["tr"])):
        v_max[i] = vehicle_model_lat(veh, tr, i)
    
    # Find apex points (local velocity minima = braking zones)
    # These are points where v_max is lower than vehicle top speed
    peak_ind = find_peaks(array2list(veh["v_max"] - v_max))[0]
    
    # Store apex data
    apex = np.zeros((len(peak_ind), 1))
    v_apex = np.zeros((len(peak_ind), 1))
    for i in range(len(peak_ind)):
        apex[i] = transpo(tr["tr"])[2][peak_ind[i]]  # Distance
        v_apex[i] = v_max[peak_ind[i]]               # Speed at apex
    
    print(f"  Found {len(peak_ind)} apex points")
    
    # =========================================================================
    # PHASE 2: FORWARD/BACKWARD INTEGRATION
    # =========================================================================
    print("[2/3] Running simulation...")
    
    N = len(apex)        # Number of apex points
    pt = len(tr["tr"])   # Number of track points
    flag = np.zeros((pt, 2))
    
    # 3D velocity array: [track_point, apex_index, mode]
    # mode: 0 = acceleration pass, 1 = deceleration pass
    v = np.ones((pt, N, 2)) * inf
    
    # For each apex, integrate both forward and backward
    for i in range(N):
        for k in range(2):  # k=0: acceleration, k=1: deceleration
            if k == 0:
                mode = 1    # Forward integration
                k_rest = 1
            else:
                mode = -1   # Backward integration
                k_rest = 0
            
            i_rest = other_points(i, N - 1)
            if len(i_rest) == 0:
                i_rest = [i]
            
            # Start at apex
            j = peak_ind[i]
            v[j][i][k] = v_apex[i]
            flag[j][k] = True
            
            # Initialize next point
            bin_j, j_next = next_point(j, pt - 1, mode)
            v[j_next][i][k] = v[j][i][k]
            j_next, j = next_point(j, pt - 1, mode)
            
            # Integrate until overshoot or lap complete
            while True:
                v[j_next][i][k], overshoot = vehicle_model_comb(
                    veh, tr, v[j][i][k], v_max[j_next], j, mode
                )
                
                if overshoot:
                    break
                
                flag[j][k] = True
                j_next, j = next_point(j, pt - 1, mode)
                
                if j == peak_ind[i]:  # Back to start
                    break
    
    # Handle boundary condition
    v[0] = v[1]
    
    print("  Simulation complete")
    
    # =========================================================================
    # PHASE 3: POST-PROCESSING
    # =========================================================================
    print("[3/3] Processing results...")
    
    # Take minimum velocity at each point across all passes
    V = np.zeros((pt, 1))
    for i in range(pt):
        V[i], ind, k = minDroiteGauche(v[i])
    
    # Calculate lap time
    temps = [0]
    for i in range(1, pt):
        temps.append(float(temps[-1] + tr['tr'][i][1] / V[i]))
    lap_time = temps[-1]
    
    # Calculate forces
    M = veh["M"]
    g = 9.81
    Fz_mass = -M * g
    Fz_aero = 0.5 * veh['rho'] * veh['factor_Cl'] * veh['Cl'] * veh['A'] * V**2
    Fz_total = Fz_aero + Fz_mass
    Fx_aero = 0.5 * veh['rho'] * veh['factor_Cd'] * veh['Cd'] * veh['A'] * V**2
    Fx_roll = veh['Cr'] * np.abs(Fz_total)
    
    # Calculate yaw rate
    yaw_rate = np.zeros((pt, 1))
    for i in range(pt):
        yaw_rate[i] = V[i] * tr['tr'][i][3]
    
    # Calculate engine data
    temp = interp1d(veh["vehicule_speed"], veh['engine_speed'])
    engine_speed = temp(V)
    
    gear = np.zeros((pt, 1))
    for i in range(pt):
        gear[i] = fast_nearest_interp(V[i], veh['vehicule_speed'], veh['gear'])
    
    # =========================================================================
    # RESULTS
    # =========================================================================
    print("\n" + "=" * 60)
    print("SIMULATION RESULTS")
    print("=" * 60)
    
    minutes = int(lap_time // 60)
    seconds = lap_time % 60
    print(f"\n  Lap Time: {minutes}:{seconds:06.3f}")
    print(f"  Avg Speed: {tr['Lcircuit'] / lap_time * 3.6:.1f} km/h")
    print(f"  Max Speed: {np.max(V) * 3.6:.1f} km/h")
    print(f"  Min Speed: {np.min(V) * 3.6:.1f} km/h")
    
    # =========================================================================
    # VISUALIZATION
    # =========================================================================
    TraceCourbeSIM()
