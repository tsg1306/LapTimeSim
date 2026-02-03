#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CircuitModelingSimV2.py - Track Geometry Processing Module

This module processes raw track/circuit data and generates a discretized track model
suitable for lap time simulation. It converts corner-by-corner track descriptions
into a fine mesh of points with associated curvature data.

The module handles:
    - Track data parsing (straights and corners with radii)
    - Track discretization into uniform segments
    - 2D track geometry reconstruction from curvature data
    - Closure error correction (ensuring the track forms a closed loop)
    - Track visualization
    - Model export for simulation

Key Concepts:
    - Curvature (κ): 1/radius for corners, 0 for straights
    - Arc length: Distance along the track centerline
    - Track mesh: Series of points with position, heading, and curvature

Author: Original development for TIPE project
Date: October 2023
Python Version: 3.7+

Dependencies:
    - numpy: Numerical computations and matrix operations
    - matplotlib: Track visualization
    - scipy: Interpolation functions

Usage:
    python CircuitModelingSimV2.py
    
    Modify the input file path to process different circuits.
    Output is saved as 'Model_<TrackName>_TrackSimV2.txt'
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

def ouvertureFileTrack(fichier):
    """
    Read and parse track/circuit data file.
    
    File format:
    - Line 1: Country<tab>TrackName
    - Subsequent lines: CornerType<tab>Length<tab>Radius
    
    Corner types:
    - 'Straight': Length in meters, Radius = 0
    - 'Left': Arc length in meters, Radius in meters (positive)
    - 'Right': Arc length in meters, Radius in meters (positive)
    
    Parameters
    ----------
    fichier : str
        Path to the track data file
        
    Returns
    -------
    tuple
        head[0] : str - Country name
        head[1] : str - Track name
        Circuit : list - List of [type, length, radius] for each segment
        
    Example
    -------
    >>> pays, nom, circuit = ouvertureFileTrack("monza.txt")
    >>> print(nom)
    'Monza'
    >>> print(circuit[0])
    ['Straight', '465.949', '0.000']
    """
    file = open(fichier, "r")
    Circuit = []
    head = (file.readline()).split("\t")
    
    for ligne in file:
        # Remove newline and split by tab
        Circuit.append(ligne[:-1].split("\t"))
    
    return head[0], head[1][:-1], Circuit


def somme(n):
    """
    Calculate sum of integers from 0 to n (triangular number).
    
    Used for calculating cumulative angle corrections.
    
    Parameters
    ----------
    n : int
        Upper limit of summation
        
    Returns
    -------
    int
        Sum: 0 + 1 + 2 + ... + n = n(n+1)/2
        
    Example
    -------
    >>> somme(5)
    15
    """
    res = 0
    for i in range(n + 1):
        res += i
    return res


def lenghtCircuit(Circuit):
    """
    Calculate total track length.
    
    Parameters
    ----------
    Circuit : list
        List of track segments [type, length, radius]
        
    Returns
    -------
    float
        Total track length in meters (rounded to 2 decimal places)
        
    Example
    -------
    >>> lenghtCircuit([['Straight', '100', '0'], ['Left', '50', '100']])
    150.0
    """
    dist = 0
    for elt in Circuit:
        dist += float(elt[1])
    return round(dist, 2)


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


def meshCircuit(Circuit, dx):
    """
    Discretize the track into uniform segments.
    
    This function converts the corner-by-corner track description into a fine
    mesh of points, each with:
    - Segment identifier (e.g., "Turn 1", "Straight 2")
    - Segment length (dx or remainder)
    - Cumulative distance from start
    - Curvature (1/radius, signed for direction)
    
    The curvature sign convention:
    - Positive: Left turn
    - Negative: Right turn
    - Zero: Straight
    
    Parameters
    ----------
    Circuit : list
        List of track segments [type, length, radius]
    dx : float
        Target segment length for discretization [m]
        
    Returns
    -------
    list
        Meshed track data: [[name, length, cumulative_dist, curvature], ...]
        
    Notes
    -----
    The function creates multiple mesh points per corner segment, with each
    point representing approximately dx meters of track. The last point of
    each segment handles any remainder length.
    """
    # Counters for naming segments
    cptDIR = {"Straight": 1, "Turn": 1}
    
    # Map corner types to standardized names
    transfoRLT = {"Straight": "Straight", "Left": "Turn", "Right": "Turn"}
    
    cptDIST = 0  # Cumulative distance
    tabRES = []  # Result array
    epsi = 0     # Sign for curvature (+1 for left, -1 for right)
    
    for elt in Circuit:
        # Determine curvature sign
        if elt[0] == "Right":
            epsi = -1
        else:
            epsi = 1
        
        # Create first point of segment (with zero length, marks segment start)
        tabRES.append([
            transfoRLT[elt[0]] + " " + str(cptDIR[transfoRLT[elt[0]]]),
            0,
            round(cptDIST, 2),
            epsi * float(elt[2])  # Signed radius
        ])
        
        # Create intermediate points at dx intervals
        for div in range(1, int(float(elt[1]) // dx) + 1):
            cptDIST += dx
            tabRES.append([
                transfoRLT[elt[0]] + " " + str(cptDIR[transfoRLT[elt[0]]]),
                dx,
                round(cptDIST, 2),
                epsi * float(elt[2])  # Signed radius
            ])
        
        # Handle remaining length (less than dx)
        cptDIST += float(elt[1]) % dx
        tabRES.append([
            transfoRLT[elt[0]] + " " + str(cptDIR[transfoRLT[elt[0]]]),
            round(float(elt[1]) % dx, 2),
            round(cptDIST, 2),
            epsi * float(elt[2])  # Signed radius
        ])
        
        # Increment segment counter
        cptDIR[transfoRLT[elt[0]]] += 1
    
    return tabRES


def finaldata(mesh):
    """
    Convert track mesh data to simulation format.
    
    Transforms the mesh data by:
    1. Removing zero-length segments (marker points)
    2. Converting radius to curvature (κ = 1/R) for corners
    3. Keeping radius as 0 for straights
    
    Parameters
    ----------
    mesh : list
        Raw mesh data from meshCircuit()
        
    Returns
    -------
    list
        Processed mesh: [[name, length, cumulative_dist, curvature], ...]
        where curvature = 1/radius for turns, 0 for straights
        
    Notes
    -----
    The curvature representation is more convenient for dynamics calculations
    since centripetal acceleration a_y = v² × κ
    """
    res = []
    for elt in mesh:
        # Skip zero-length marker points
        if elt[1] != 0:
            if "Straight" in elt[0]:
                res.append(elt)  # Keep as-is (curvature = 0)
            else:
                # Convert radius to curvature: κ = 1/R
                res.append(elt[:-1] + [1 / elt[3]])
    return res


def TraceCircuit(X, Y):
    """
    Visualize the track layout.
    
    Creates a 2D scatter plot of the track with:
    - Blue points showing the track centerline
    - Red triangle marker at start/finish
    
    Parameters
    ----------
    X : list or array
        X coordinates of track points [m]
    Y : list or array
        Y coordinates of track points [m]
    """
    plt.figure(figsize=(12, 8))
    plt.scatter(X, Y, s=4, c='blue', label='Track centerline')
    
    # Mark start/finish line with a red triangle
    plt.plot((X[-1] + X[0]) / 2, (Y[-1] + Y[0]) / 2, 
             marker="<", color="red", markersize=12, label='Start/Finish')
    
    plt.xlabel('X [m]', fontsize=14)
    plt.ylabel('Y [m]', fontsize=14)
    plt.title('Track Layout')
    plt.axis('equal')  # Equal aspect ratio for accurate representation
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


def alphacorr(alpha):
    """
    Normalize heading angle to [-π, π] range.
    
    Circuit closure requires that the total heading change equals exactly 2π
    (or -2π for clockwise tracks). This function normalizes the cumulative
    heading angle to identify the closure error.
    
    Parameters
    ----------
    alpha : float
        Cumulative heading angle [rad]
        
    Returns
    -------
    float
        Normalized angle in [-π, π]
        
    Notes
    -----
    For a properly closed circuit running counter-clockwise, the total heading
    change should be 2π. Any deviation is the closure error that needs to be
    distributed across all segments.
    """
    if abs(alpha % (2 * pi * (alpha / abs(alpha)))) < pi:
        return alpha % (2 * pi * (alpha / abs(alpha)))
    else:
        return (abs(alpha) / alpha) * (2 * pi) - alpha % (2 * pi * (alpha / abs(alpha)))


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    
    # =========================================================================
    # SECTION 1: READ TRACK DATA
    # =========================================================================
    print("=" * 60)
    print("LapTimeSim - Track Model Generator")
    print("=" * 60)
    
    # Read track data file
    # Modify this path to process different circuits
    pays, nom_circuit, circuit = ouvertureFileTrack("../data/tracks/monza.txt")
    
    print(f"\nProcessing circuit: {nom_circuit} ({pays})")
    
    # =========================================================================
    # SECTION 2: INITIALIZE PARAMETERS
    # =========================================================================
    
    # Export dictionary for storing all computed data
    export = {}
    
    # Calculate and store total track length
    L = lenghtCircuit(circuit)
    export["Lcircuit"] = L
    print(f"Track length: {L} m ({L/1000:.2f} km)")
    
    # Discretization frequency (segment length in meters)
    f = 5  # 5 meter segments
    
    # =========================================================================
    # SECTION 3: TRACK DISCRETIZATION
    # =========================================================================
    # Convert corner-by-corner data to fine mesh
    
    mesh = meshCircuit(circuit, f)
    print(f"✓ Track discretization complete ({len(mesh)} points)")
    
    # =========================================================================
    # SECTION 4: TRACK GEOMETRY RECONSTRUCTION
    # =========================================================================
    # Convert curvature data to X,Y coordinates
    #
    # The algorithm integrates the curvature to build the track geometry:
    # 1. At each point, update heading angle based on curvature and arc length
    # 2. Move in the current heading direction by the segment length
    # 3. Apply corrections to ensure track closure
    
    # Initialize coordinate and heading arrays
    X = [0]      # X coordinates [m]
    Y = [0]      # Y coordinates [m]
    alpha = [0]  # Heading angles [rad]
    
    # -------------------------------------------------------------------------
    # First Pass: Raw Geometry (uncorrected)
    # -------------------------------------------------------------------------
    # Integrate curvature to get heading, then heading to get position
    
    for i in range(1, len(mesh)):
        if "Turn" in mesh[i][0]:
            if mesh[i][1] != 0:  # Non-zero segment length
                # Update heading: Δα = arc_length / radius
                # For turns, mesh[i][3] is the signed radius
                alpha.append(alpha[-1] + mesh[i][1] / mesh[i][3])
                
                # Calculate position change using rotation matrix
                # [dx]   [cos(α) -sin(α)] [ds]
                # [dy] = [sin(α)  cos(α)] [0 ]
                temp = np.matmul(
                    np.array([[cos(alpha[-1]), -sin(alpha[-1])],
                              [sin(alpha[-1]),  cos(alpha[-1])]]),
                    np.array([mesh[i][1], 0])
                )
                
                # Update position (negative because we defined direction convention)
                X.append(-temp[0] + X[-1])
                Y.append(-temp[1] + Y[-1])
                
        elif "Straight" in mesh[i][0]:
            if mesh[i][1] != 0:
                # Straights: heading doesn't change, just move forward
                temp = np.matmul(
                    np.array([[cos(alpha[-1]), -sin(alpha[-1])],
                              [sin(alpha[-1]),  cos(alpha[-1])]]),
                    np.array([mesh[i][1], 0])
                )
                X.append(-temp[0] + X[-1])
                Y.append(-temp[1] + Y[-1])
    
    # -------------------------------------------------------------------------
    # Apply Angular Correction
    # -------------------------------------------------------------------------
    # The cumulative heading should be exactly ±2π for a closed circuit.
    # Any deviation is distributed linearly across all points.
    
    # Calculate angular correction per point
    alpha_corr = alphacorr(alpha[-1]) / (len(alpha) - 1)
    
    # Apply linear correction to all heading angles
    for i in range(len(alpha)):
        alpha[i] = alpha[i] - i * alpha_corr
    
    # Recompute positions with corrected headings
    X = [0]
    Y = [0]
    cpt_alpha = 0  # Index for corrected alpha array
    
    for i in range(1, len(mesh)):
        if "Turn" in mesh[i][0]:
            if mesh[i][1] != 0:
                temp = np.matmul(
                    np.array([[cos(alpha[cpt_alpha]), -sin(alpha[cpt_alpha])],
                              [sin(alpha[cpt_alpha]),  cos(alpha[cpt_alpha])]]),
                    np.array([mesh[i][1], 0])
                )
                X.append(-temp[0] + X[-1])
                Y.append(-temp[1] + Y[-1])
                cpt_alpha += 1
                
        elif "Straight" in mesh[i][0]:
            if mesh[i][1] != 0:
                temp = np.matmul(
                    np.array([[cos(alpha[cpt_alpha]), -sin(alpha[cpt_alpha])],
                              [sin(alpha[cpt_alpha]),  cos(alpha[cpt_alpha])]]),
                    np.array([mesh[i][1], 0])
                )
                X.append(-temp[0] + X[-1])
                Y.append(-temp[1] + Y[-1])
    
    # -------------------------------------------------------------------------
    # Apply Position Correction
    # -------------------------------------------------------------------------
    # After angular correction, there may still be a small gap between
    # start and end points. Distribute this error linearly.
    
    X_cor = (X[-1] - X[0]) / len(X)  # X correction per point
    Y_cor = (Y[-1] - Y[0]) / len(Y)  # Y correction per point
    
    # Final position calculation with both corrections
    X = [0]
    Y = [0]
    cpt_alpha = 0
    
    for i in range(1, len(mesh)):
        if "Turn" in mesh[i][0]:
            if mesh[i][1] != 0:
                temp = np.matmul(
                    np.array([[cos(alpha[cpt_alpha]), -sin(alpha[cpt_alpha])],
                              [sin(alpha[cpt_alpha]),  cos(alpha[cpt_alpha])]]),
                    np.array([mesh[i][1], 0])
                )
                X.append(-temp[0] + X[-1] - X_cor)
                Y.append(-temp[1] + Y[-1] - Y_cor)
                cpt_alpha += 1
                
        elif "Straight" in mesh[i][0]:
            if mesh[i][1] != 0:
                temp = np.matmul(
                    np.array([[cos(alpha[cpt_alpha]), -sin(alpha[cpt_alpha])],
                              [sin(alpha[cpt_alpha]),  cos(alpha[cpt_alpha])]]),
                    np.array([mesh[i][1], 0])
                )
                X.append(-temp[0] + X[-1] - X_cor)
                Y.append(-temp[1] + Y[-1] - Y_cor)
    
    # Store coordinates in export dictionary
    export["X"] = X
    export["Y"] = Y
    
    # Process mesh data to simulation format (convert radius to curvature)
    tr = finaldata(mesh)
    export['tr'] = tr
    
    print(f"✓ Geometry reconstruction complete")
    print(f"  Closure error (X): {abs(X[-1] - X[0]):.4f} m")
    print(f"  Closure error (Y): {abs(Y[-1] - Y[0]):.4f} m")
    
    # =========================================================================
    # SECTION 5: VISUALIZATION
    # =========================================================================
    TraceCircuit(X, Y)
    print("✓ Visualization complete")
    
    # =========================================================================
    # SECTION 6: SAVE MODEL
    # =========================================================================
    print("\n" + "=" * 60)
    print("Track Modeling Complete!")
    print("=" * 60)
    
    # Save model file
    output_filename = f"Model_{nom_circuit}_TrackSimV2.txt"
    save = open(output_filename, 'w')
    save.write(nom_circuit + "\t" + pays + '\n')
    save.write(str(export))
    save.close()
    
    print(f"\n✓ Model exported as: {output_filename}")
    print(f"\nTrack Summary:")
    print(f"  Name: {nom_circuit}")
    print(f"  Country: {pays}")
    print(f"  Length: {L} m ({L/1000:.2f} km)")
    print(f"  Segments: {len(circuit)}")
    print(f"  Mesh points: {len(tr)}")
