# Physics Model Documentation

This document provides a detailed explanation of the physics models used in LapTimeSim.

## Table of Contents

1. [Overview](#overview)
2. [Coordinate System](#coordinate-system)
3. [Force Model](#force-model)
4. [Tire Model](#tire-model)
5. [Aerodynamic Model](#aerodynamic-model)
6. [Powertrain Model](#powertrain-model)
7. [Simulation Algorithm](#simulation-algorithm)
8. [Equations Reference](#equations-reference)

---

## Overview

LapTimeSim uses a **quasi-static** approach to lap time simulation. This means we assume the vehicle is always in equilibrium at each point along the track - there are no transient dynamics like oscillations or weight transfer settling time.

The key assumptions are:
- Vehicle follows the racing line (track centerline)
- Driver uses 100% of available grip at all times
- Gear shifts are instantaneous (or with minimal time penalty)
- No driver error or reaction time

---

## Coordinate System

We use a vehicle-fixed coordinate system:

```
        Z (up)
        |
        |
        +----> X (forward/longitudinal)
       /
      /
     Y (left/lateral)
```

**Sign conventions:**
- Positive X: Forward acceleration
- Negative X: Braking (deceleration)
- Positive Y: Left turn
- Negative Y: Right turn
- Positive Z: Upward (lift)
- Negative Z: Downward (downforce/weight)

---

## Force Model

### Vertical Forces (Z-axis)

The total normal force on the tires determines available grip:

```
F_z = F_weight + F_downforce
```

Where:
- `F_weight = -M × g` (always pointing down, hence negative)
- `F_downforce = 0.5 × ρ × C_L × A × v²` (negative C_L = downforce)

The normal force on driven wheels specifically:

```
F_z_driven = (f_drive × F_weight + f_aero × F_downforce) / n_wheels
```

Where:
- `f_drive` = fraction of weight on driven wheels
- `f_aero` = fraction of aero load on driven wheels
- `n_wheels` = number of driven wheels (2 or 4)

### Longitudinal Forces (X-axis)

```
ΣF_x = F_engine - F_drag - F_rolling - F_brake
```

Or in terms of acceleration:

```
M × a_x = F_traction - F_drag - F_rolling
```

### Lateral Forces (Y-axis)

For cornering at velocity v with curvature κ = 1/R:

```
F_y = M × v² × κ  (centripetal force required)
```

This must equal the available lateral tire force.

---

## Tire Model

### Load-Sensitive Friction

Real tires exhibit decreasing friction coefficient with increasing load. We model this as:

```
μ = μ_0 + sensitivity × (N_rated - N_actual)
```

Where:
- `μ_0` = base friction coefficient
- `sensitivity` = load sensitivity [1/N]
- `N_rated` = tire's rated load
- `N_actual` = actual normal force on tire

### Friction Ellipse (Combined Grip)

When cornering and accelerating/braking simultaneously, the available grip is shared:

```
(a_x / a_x_max)² + (a_y / a_y_max)² ≤ 1
```

This is called the **friction ellipse** or **traction circle**.

In implementation:
```python
if ay != 0:
    ay_max = μ_y × N / M
    ellipse_multiplier = sqrt(1 - (ay / ay_max)²)
else:
    ellipse_multiplier = 1

ax_available = ax_max × ellipse_multiplier
```

### Maximum Accelerations

**Lateral (cornering):**
```
a_y_max = (1/M) × (μ_y + dμ_y × (N_rated - N/4)) × N
```

**Longitudinal (acceleration):**
```
a_x_max_acc = (1/M) × (μ_x + dμ_x × (N_rated - N_driven)) × N_driven × n_driven
```

**Longitudinal (braking):**
```
a_x_max_dec = -(1/M) × (μ_x + dμ_x × (N_rated - N/4)) × N
```

Note: Braking uses all 4 wheels, acceleration uses only driven wheels.

---

## Aerodynamic Model

### Drag Force

```
F_drag = 0.5 × ρ × C_D × A × v²
```

Where:
- `ρ` = air density [kg/m³] (typically 1.225)
- `C_D` = drag coefficient
- `A` = frontal area [m²]
- `v` = velocity [m/s]

### Downforce

```
F_downforce = 0.5 × ρ × C_L × A × v²
```

Where `C_L` is negative for downforce (using lift convention).

### Rolling Resistance

```
F_rolling = C_r × |F_z_total|
```

Where `C_r` is the rolling resistance coefficient (typically 0.01-0.02).

---

## Powertrain Model

### Engine Torque to Wheel Force

```
F_wheel = T_engine × i_gear × i_final × i_primary × η_total / r_wheel
```

Where:
- `T_engine` = engine torque [N·m]
- `i_gear` = current gear ratio
- `i_final` = final drive ratio
- `i_primary` = primary reduction ratio
- `η_total` = total drivetrain efficiency
- `r_wheel` = wheel radius [m]

### Gear Selection

At each velocity, the optimal gear is selected to maximize traction force:

```python
for each gear:
    F_x[gear] = interpolate(v_gear, T_gear × i / r) at current_v
optimal_gear = argmax(F_x)
```

### Vehicle Speed vs Engine Speed

```
v = ω_engine × r_wheel / (i_gear × i_final × i_primary × (2π/60))
```

Or:
```
RPM = v × (i_gear × i_final × i_primary × 60) / (r_wheel × 2π)
```

---

## Simulation Algorithm

### Step 1: Maximum Cornering Speed

For each track point, solve for the maximum velocity where available lateral grip equals required centripetal force:

```
M × v² × κ = F_y_max(v)
```

This is a quadratic in v² because F_y_max depends on downforce which is proportional to v².

### Step 2: Find Apex Points

Apex points are local minima in the v_max profile. These correspond to the slowest points through corners.

```python
apex_indices = find_peaks(v_top_speed - v_max)
```

### Step 3: Forward Integration (Acceleration)

From each apex, integrate forward:

```python
for each apex:
    v = v_apex
    for each point forward:
        a_x = min(a_engine_limit, a_tire_limit × ellipse_factor) - a_drag
        v_next = sqrt(v² + 2 × a_x × dx)
        if v_next > v_max[next_point]:
            break  # Would overshoot
        v = v_next
```

### Step 4: Backward Integration (Braking)

From each apex, integrate backward (simulating braking approach):

```python
for each apex:
    v = v_apex
    for each point backward:
        a_x = -a_brake_limit × ellipse_factor - a_drag
        v_prev = sqrt(v² - 2 × a_x × dx)
        if v_prev > v_max[prev_point]:
            break  # Would overshoot
        v = v_prev
```

### Step 5: Merge Profiles

At each point, take the minimum velocity from all integration passes:

```python
V[i] = min(all v[i] from all apex passes)
```

### Step 6: Calculate Lap Time

```python
lap_time = sum(dx[i] / V[i] for all i)
```

---

## Equations Reference

### Quick Reference Card

| Quantity | Equation |
|----------|----------|
| Downforce | `F_df = 0.5 × ρ × C_L × A × v²` |
| Drag | `F_dr = 0.5 × ρ × C_D × A × v²` |
| Centripetal | `F_c = M × v² / R = M × v² × κ` |
| Friction | `μ = μ_0 + sens × (N_rated - N)` |
| Max lateral a | `a_y = μ × F_z / M` |
| Ellipse | `(a_x/a_x_max)² + (a_y/a_y_max)² = 1` |
| Kinematics | `v² = v_0² + 2 × a × dx` |
| Lap time | `t = Σ(dx_i / v_i)` |

### Typical Values

| Parameter | Typical Range | Units |
|-----------|---------------|-------|
| C_L (GT car) | -2 to -4 | - |
| C_D (GT car) | -0.8 to -1.2 | - |
| C_L (F1) | -4 to -6 | - |
| μ_tire (slicks) | 1.5 to 2.0 | - |
| ρ (sea level) | 1.225 | kg/m³ |
| C_r (racing) | 0.008 to 0.015 | - |

---

## References

1. Milliken, W. & Milliken, D. (1995). Race Car Vehicle Dynamics
2. Pacejka, H. (2012). Tire and Vehicle Dynamics
3. Smith, C. (1978). Tune to Win
