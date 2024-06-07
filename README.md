# Flax and Teal: Solving the Burgers' Equation Using RK4 Method

This repository contains the implementation of the Burgers' equation solver using the Runge-Kutta method of fourth order (RK4) in Python. The project includes the mathematical formulation, numerical method, Python code, and output results.

## Table of Contents

- [Introduction](#introduction)
- [Mathematical Formulation](#mathematical-formulation)
- [Numerical Method](#numerical-method)
- [Python Code](#python-code)
- [Results](#results)
- [Usage](#usage)
- [Requirements](#requirements)
- [License](#license)

## Introduction

The Burgers' equation is a fundamental partial differential equation from fluid mechanics, combining nonlinear convection and diffusion. It serves as a model for various physical phenomena. This project demonstrates how to numerically solve the Burgers' equation using the RK4 method.

## Mathematical Formulation

The Burgers' equation is given by:

\[ u_t + u u_x = \nu u_{xx} \]

where:
- \(u\) is the velocity field
- \(t\) is time
- \(x\) is the spatial coordinate
- \(\nu\) is the viscosity coefficient

## Numerical Method

To solve this equation numerically:
1. Discretize the spatial domain into a grid.
2. Use finite difference approximations for the spatial derivatives.
3. Apply the RK4 method to update the solution in time.

## Python Code

The following Python code is used to solve the Burgers' equation:

```python
import numpy as np

n = 10  # Number of grid points
y0 = np.cos(np.linspace(0, 2*np.pi, n))
t0 = 0  # Initial time
t_end = 2  # Final time
nu = 0.1  # Viscosity coefficient

# Define the function f(y, t)
def f(y, t, nu=0.1):
    n = len(y)
    dydt = np.zeros_like(y)

    # Compute dy/dt for each point using finite differences
    for i in range(n):
        if i == 0:
            dydt[i] = -y[i] * (y[i+1] - y[n-1]) / 2 + nu * (y[n-1] - 2 * y[i] + y[i+1]) / (n-1)**2
        elif i == n-1:
            dydt[i] = -y[i] * (y[0] - y[i-1]) / 2 + nu * (y[i-1] - 2 * y[i] + y[0]) / (n-1)**2
        else:
            dydt[i] = -y[i] * (y[i+1] - y[i-1]) / 2 + nu * (y[i-1] - 2 * y[i] + y[i+1]) / (n-1)**2

    return dydt

# Run RK4 method
def RK4(f, n, y0, t0, t_end):
    t_values = [t0]
    y_values = [y0]
    t = t0
    y = y0
    h = (t_end - t0) / n  # Step size
    while t < t_end:
        k1 = f(y, t)
        k2 = f(y + 0.5 * h * k1, t + 0.5 * h)
        k3 = f(y + 0.5 * h * k2, t + 0.5 * h)
        k4 = f(y + h * k3, t + h)
        y = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        t = t + h
        t_values.append(t)
        y_values.append(y)
    return np.array(t_values), np.array(y_values)

t_values, y_values = RK4(f, n, y0, t0, t_end)
