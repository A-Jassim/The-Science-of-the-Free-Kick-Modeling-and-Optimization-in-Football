# ⚽ The Science of the Free Kick: Modeling and Optimization in Football

Author: **Akram HALIMI**  
Academic Year: **2023–2024**

## Overview

This project explores the physics, modeling, and simulation of football free kicks, with a focus on improving player performance through scientific analysis. Using computational tools and physics-based equations, the project simulates realistic ball trajectories under various impact conditions — including aerodynamic forces like drag and the Magnus effect.

The approach involves:
- Analytical modeling of the physical forces acting on a spinning football
- Numerical resolution of the governing differential equations (Runge-Kutta method)
- Monte Carlo simulations for optimization of shooting parameters
- Reproduction of real-world free kicks, including Roberto Carlos' iconic 35m goal

## Project Components

- `rapport_en.pdf` — Full technical report detailing the theoretical framework, equations, simulation approach, and analysis of results.
- `Presentation.pdf` — Visual presentation summarizing the methodology, physical principles, and key findings.
- `courbe version test en legende.py` — Python simulation script:
  - Models the ball's motion using Newtonian dynamics
  - Calculates optimal impact angles and spin
  - Visualizes 3D trajectories of free kicks
  - Highlights valid goal trajectories vs. failed attempts

## Features

- Physics-based projectile modeling with drag and Magnus effects
- Accurate 3D goal and wall visualization using `matplotlib`
- Stochastic simulation using Monte Carlo methods to explore variability
- Reproducibility of professional-level kicks from input parameters
- Impact map visualization on the ball surface

## Technologies

- Python 3
- NumPy
- Matplotlib (2D + 3D plots)

## Sample Output

- Identification of optimal impact coordinates on the ball
- Simulation of hundreds of randomized shots to statistically optimize parameters
- Visual confirmation of goal entry based on official goalpost dimensions

## How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/A-Jassim/The-Science-of-the-Free-Kick-Modeling-and-Optimization-in-Football.git
   cd The-Science-of-the-Free-Kick-Modeling-and-Optimization-in-Football
   ```

2. Run the simulation:
   ```bash
   python -u "script.py"
   ```

3. View plots and read printed analysis in the terminal.

## Acknowledgments

- Inspired by physics research in sports biomechanics
- Validation based on historical football free kicks (e.g., Roberto Carlos)
- Supported by Regressi and academic simulations

---