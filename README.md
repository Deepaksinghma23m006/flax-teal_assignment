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



# Results

The following output was obtained from running the above Python code:

At t = 0 value of y is = [ 1.          0.76604444  0.17364818 -0.5        -0.93969262 -0.93969262
 -0.5         0.17364818  0.76604444  1.        ]
for t= 0.2 value of y is = [ 1.01896957  0.83179839  0.19828195 -0.56055162 -0.9767954  -0.89507841
 -0.44871143  0.15385997  0.70540544  0.97282155] 

for t= 0.4 value of y is = [ 1.02813156  0.9019579   0.2294544  -0.63219901 -1.00393829 -0.84584025
 -0.40532041  0.13777068  0.65039073  0.9395927 ] 

for t= 0.6000000000000001 value of y is = [ 1.0263321   0.97492956  0.26952071 -0.71701407 -1.01865531 -0.79445568
 -0.36852488  0.12450979  0.60092067  0.90243711] 

for t= 0.8 value of y is = [ 1.01317279  1.04818259  0.32185739 -0.81752954 -1.01882251 -0.74311209
 -0.33721425  0.11344028  0.55667433  0.86335101] 

for t= 1.0 value of y is = [ 0.98913423  1.11804652  0.39137643 -0.93700851 -1.00281008 -0.69363051
 -0.31044632  0.10409067  0.51719651  0.82405106] 

for t= 1.2 value of y is = [ 0.95560192  1.17952005  0.485344   -1.07997456 -0.96959634 -0.64743527
 -0.28742668  0.09610806  0.48197561  0.78588322] 

for t= 1.4 value of y is = [ 0.91478965  1.22610241  0.61473165 -1.25325385 -0.91882513 -0.60556558
 -0.26748945  0.08922538  0.45049654  0.74978839] 

for t= 1.5999999999999999 value of y is = [ 0.86957853  1.24966661  0.79657071 -1.46805224 -0.85078942 -0.56872015
 -0.25007916  0.08323831  0.42227408  0.71631273] 

for t= 1.7999999999999998 value of y is = [ 0.8233109   1.24042811  1.05840659 -1.74426719 -0.76632747 -0.53732428
 -0.23473413  0.07798882  0.39687169  0.68564696] 

for t= 1.9999999999999998 value of y is = [ 0.77959231  1.18714043  1.44780212 -2.12016886 -0.66662805 -0.51160953
 -0.22107119  0.07335332  0.37391008  0.65767937]


# Usage
