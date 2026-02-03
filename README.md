# 🏎️ LapTimeSim - Quasi-Static Lap Time Simulator

[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A physics-based lap time simulation tool that models vehicle dynamics to predict optimal lap times on racing circuits. Originally developed as part of a TIPE (Travaux d'Initiative Personnelle Encadrés) project.

<p align="center">
  <img src="docs/circuit_speed_map.png" alt="Circuit Speed Map" width="600"/>
</p>

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Physics Model](#physics-model)
- [Installation](#installation)
- [Usage](#usage)
- [File Formats](#file-formats)
- [Project Structure](#project-structure)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## 🎯 Overview

LapTimeSim is a quasi-static lap time simulator that calculates the theoretical optimal lap time for a given car-track combination. The simulation uses fundamental physics principles including:

- **Tire friction modeling** with load sensitivity
- **Aerodynamic forces** (downforce and drag)
- **Powertrain modeling** with gear ratios and efficiency
- **Combined grip** using the friction ellipse (traction circle) concept
- **Optimal velocity profiling** through corners and straights

### What is a Quasi-Static Simulation?

Unlike dynamic simulations that solve differential equations in real-time, a quasi-static approach assumes the vehicle is always in equilibrium at each point. This simplification allows for faster computation while still capturing the essential physics of vehicle performance.

## ✨ Features

- 🚗 **Comprehensive Vehicle Model**: Mass distribution, aerodynamics, tires, powertrain
- 🛤️ **Track Discretization**: Converts track centerline data into simulation-ready mesh
- 📊 **GGV Diagram Generation**: 3D traction envelope visualization
- 🔄 **Optimal Gear Selection**: Automatic gear optimization based on engine curves
- 📈 **Rich Visualization**: Speed traces, circuit maps with velocity coloring, gear/RPM plots
- 💾 **Modular Data Format**: Separate car and track model files for easy mixing and matching

## 🔬 Physics Model

### Vehicle Dynamics

The simulator models the following forces:

#### Longitudinal Forces (X-axis)
```
F_x = F_engine - F_aero_drag - F_rolling - F_braking
```

- **Engine Force**: Interpolated from torque curve × gear ratio × efficiency
- **Aerodynamic Drag**: `F_drag = 0.5 × ρ × Cd × A × v²`
- **Rolling Resistance**: `F_roll = Cr × |F_z|`

#### Vertical Forces (Z-axis)
```
F_z = -M × g + F_downforce
```

- **Weight**: `F_weight = M × g`
- **Downforce**: `F_down = 0.5 × ρ × Cl × A × v²`

#### Lateral Forces (Y-axis)
```
F_y = M × v² / R  (centripetal force required for cornering)
```

### Tire Model

The tire friction model includes load sensitivity:

```
μ = μ_0 + sensitivity × (N_rated - N_actual)
```

Where:
- `μ_0` is the base friction coefficient
- `sensitivity` is the load sensitivity factor
- `N_rated` is the rated load for the tire
- `N_actual` is the actual normal load on the tire

### Friction Ellipse (Combined Grip)

The combined longitudinal and lateral grip is modeled using the friction ellipse:

```
(a_x / a_x_max)² + (a_y / a_y_max)² ≤ 1
```

This ensures that when the car is cornering at maximum lateral grip, longitudinal acceleration capacity is reduced accordingly.

### Simulation Algorithm

1. **Calculate Maximum Cornering Speeds**: For each point, compute the maximum velocity limited by lateral grip
2. **Find Apex Points**: Identify local minima in the velocity profile (braking zones)
3. **Forward Integration**: From each apex, simulate acceleration forward
4. **Backward Integration**: From each apex, simulate braking backward
5. **Merge Profiles**: Take the minimum velocity at each point from all integration passes

## 🛠️ Installation

### Prerequisites

- Python 3.7 or higher
- pip package manager

### Dependencies

```bash
pip install numpy matplotlib scipy
```

### Clone the Repository

```bash
git clone https://github.com/yourusername/LapTimeSim.git
cd LapTimeSim
```

## 🚀 Usage

### Quick Start

1. **Model your car** (or use an existing model):
```bash
python src/CarModelSimV2.py
```

2. **Model your track** (or use an existing model):
```bash
python src/CircuitModelingSimV2.py
```

3. **Run the lap simulation**:
```bash
python src/LapSimV2.py
```

### Configuration

Edit the file paths in each script to point to your data files:

```python
# In LapSimV2.py
nom_circuit, pays, tr = ouvertureFileTrack("data/tracks/Model_Monza_TrackSimV2.txt")
nom_voiture, veh = ouvertureFileCar("data/cars/Model_Alpine_A110_CarSimV2.txt")
```

## 📁 File Formats

### Car Input File (`.txt`)

Tab-separated values with parameter names and values:

```
Alpine A110 (cup)
Total Mass	1100	kg
Front Mass Distribution	42	%
Wheelbase	2420	mm
Lift Coefficient CL	-2.8	-
Drag Coefficient CD	-0.8	-
...
MARKER_info/data
1000	150
2000	200
...
```

The section after `MARKER_info/data` contains the engine torque curve (RPM vs Torque in N·m).

### Track Input File (`.txt`)

Tab-separated with corner type, length, and radius:

```
Country	TrackName
Straight	465.949	0.000
Left	5.364	1856.545
Right	4.046	787.304
...
```

- **Straight**: Length in meters, radius = 0
- **Left/Right**: Arc length in meters, radius in meters

### Model Output Files

Both car and track models are exported as Python dictionaries serialized to text files, containing all computed parameters needed for simulation.

## 📂 Project Structure

```
LapTimeSim/
├── README.md
├── LICENSE
├── requirements.txt
├── src/
│   ├── CarModelSimV2.py      # Vehicle parameter processing
│   ├── CircuitModelingSimV2.py   # Track geometry processing
│   └── LapSimV2.py           # Main lap time simulation
├── data/
│   ├── cars/                 # Car definition files
│   │   ├── AlpineA110V2.txt
│   │   └── essaieCARV2.txt
│   └── tracks/               # Track definition files
│       ├── castelletV2.txt
│       └── monza.txt
├── docs/
│   └── physics_model.md      # Detailed physics documentation
├── examples/
│   └── example_usage.py      # Example usage scripts
└── output/
    └── Model_*_SimV2.txt     # Generated model files
```

## 📊 Examples

### Sample Output

Running the simulation on Paul Ricard (Le Castellet) with an Alpine A110 Cup:

```
Découpage du circuit : OK
Simulation terminé
Temps au tour: 1:52.340
```

### Visualization

The simulator generates several plots:

1. **Speed vs Time**: Vehicle speed throughout the lap
2. **Circuit Map**: Track layout colored by velocity
3. **Engine RPM & Gear**: Gear selection and engine speed
4. **Yaw Rate**: Vehicle rotation rate through corners
5. **GGV Diagram**: 3D traction envelope (Car Model only)

## 🔧 Customization

### Adding a New Car

1. Create a new car definition file following the format in `data/cars/`
2. Run `CarModelSimV2.py` with the new file path
3. Use the generated model file in `LapSimV2.py`

### Adding a New Track

1. Create track data with corner-by-corner radius and length
2. Format as per `data/tracks/` examples
3. Run `CircuitModelingSimV2.py` to generate the model

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

### Development Ideas

- [ ] Add elevation changes (3D track model)
- [ ] Implement thermal tire model
- [ ] Add fuel consumption modeling
- [ ] Create GUI interface
- [ ] Support for telemetry data import
- [ ] Validate against real lap times

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Developed as part of a TIPE project on racing simulation
- Inspired by professional lap time simulation tools
- Physics model based on vehicle dynamics literature

## 📧 Contact

For questions or feedback, please open an issue on GitHub.

---

*Made with ❤️ for motorsport enthusiasts and engineering students*
