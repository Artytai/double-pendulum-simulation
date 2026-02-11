# Double Pendulum Simulation

A real-time physics simulation of a double pendulum built from first principles — no physics engines used.

This project models one of the most famous chaotic systems in classical mechanics and numerically integrates the equations of motion using a fourth-order Runge–Kutta (RK4) method.

---

## Overview

The double pendulum is a nonlinear dynamical system that exhibits deterministic chaos. While governed by simple physical laws, its motion becomes highly sensitive to initial conditions.

This simulation:

- Implements the full coupled equations of motion
- Integrates them using RK4 for numerical stability
- Renders motion using HTML5 Canvas
- Displays live telemetry and total mechanical energy
- Allows real-time parameter adjustments

---

## Physics Model

The system consists of:

- Two masses: m₁ and m₂
- Two rigid rods of lengths ℓ₁ and ℓ₂
- Gravitational acceleration g
- Optional damping (air resistance)

The equations of motion are derived using Lagrangian mechanics and result in two coupled second-order nonlinear differential equations.

They are converted into four first-order equations:

θ₁' = ω₁  
θ₂' = ω₂  
ω₁' = f(θ₁, θ₂, ω₁, ω₂)  
ω₂' = g(θ₁, θ₂, ω₁, ω₂)

These are numerically integrated using:

### Numerical Method
**Fourth-Order Runge–Kutta (RK4)**

RK4 was chosen over Euler integration due to:

- Higher stability
- Reduced energy drift
- Better long-term accuracy

---

## Features

- Real-time simulation rendering
- Adjustable parameters:
  - Gravity
  - Masses
  - Rod lengths
  - Time step (dt)
  - Drag / damping
- Live telemetry display:
  - Time
  - Angular positions (θ₁, θ₂)
  - Angular velocities (ω₁, ω₂)
  - Total mechanical energy
- Energy conservation tracking
- Clean dark-mode UI
- Responsive layout

---

## Controls

| Control | Description |
|----------|------------|
| Play | Starts the simulation |
| Pause | Stops integration |
| Reset | Resets to initial conditions |
| Sliders | Adjust physical parameters in real time |
| dt | Adjust integration time step |

---

## Demonstrated Concepts

- Deterministic chaos
- Sensitivity to initial conditions
- Nonlinear coupled differential equations
- Numerical integration methods
- Energy conservation in mechanical systems

---

---

## Why This Is Interesting

The double pendulum is a classic example of how:

> Simple equations can produce incredibly complex behavior.

Even extremely small changes in initial angle can result in dramatically different motion over time — a hallmark of chaotic systems.

---

## Possible Extensions

- 3D double pendulum
- Triple pendulum
- Phase space visualization
- Lyapunov exponent calculation
- Collision constraints
- Elastic rods instead of rigid rods
- GPU-based rendering

---

## Author

Arthur Tai  
Built as a computational physics and numerical simulation project.

