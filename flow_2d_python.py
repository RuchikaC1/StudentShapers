# =========================================================================
# Final Corrected Script: Transient Two-Phase Stokes Flow using FiPy
# =========================================================================
# This script simulates a heavy fluid (liquid) displacing a light fluid (gas)
# around a cylindrical obstacle at low Reynolds number (Stokes Flow).
# It uses the Level-Set method to track the interface and a projection
# method to solve for pressure and velocity.
#
# FINAL MODIFICATION: This version uses a non-staggered (colocated) grid
# for velocity (velocity is a CellVariable) and correctly applies all
# boundary conditions to resolve previous errors.
# =========================================================================

import os
import numpy as np
import matplotlib.pyplot as plt
from fipy import (CellVariable, FaceVariable, Grid2D, TransientTerm,
                  ImplicitDiffusionTerm, ConvectionTerm, Viewer,
                  ImplicitSourceTerm, numerix)

print("Setting up FiPy simulation...")

# --- 1. Simulation Parameters ---

# Domain and Grid
nx = 80
ny = 40
dx = 0.5
dy = 0.5
mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)

# Fluid Properties (Liquid and Gas)
rho1, rho2 = 100., 1.0  # Density
mu1, mu2 = 10., 0.1     # Viscosity

# Flow and Time Parameters
inlet_velocity = 1.0
dt = 0.05
steps = 400

# Interface Smoothing Parameter
epsilon = 2.0 * dx

# Cylinder Parameters
cylinder_radius = 8.0 * dx
cylinder_center = (nx * dx / 3.0, ny * dy / 2.0)

# --- 2. Define Variables ---
# MODIFIED: Velocity is a CellVariable to allow for the use of TransientTerm.

# Level-Set function (phi < 0 is liquid, phi > 0 is gas)
phi = CellVariable(name=r"$\phi$", mesh=mesh, hasOld=True)

# Pressure
pressure = CellVariable(name="p", mesh=mesh, hasOld=True)

# Velocity (defined on the cell centers)
velocity = CellVariable(name="u", mesh=mesh, rank=1, hasOld=True)

# Fluid Properties (defined on cells, dependent on phi)
density = CellVariable(name=r"$\rho$", mesh=mesh)
viscosity = CellVariable(name=r"$\mu$", mesh=mesh)


# --- 3. Initial Conditions and Cylinder Mask ---

# Initialize phi to define the initial vertical front
X, Y = mesh.cellCenters
phi.setValue(X - nx * dx / 5.0)

# Create a mask for the cylinder obstacle
cell_centers = mesh.cellCenters
x_cell, y_cell = cell_centers
cylinder_mask = (x_cell - cylinder_center[0])**2 + (y_cell - cylinder_center[1])**2 < cylinder_radius**2

# --- 4. Define the Physics (PDEs) ---

# -- Smoothed Fluid Properties --
H = 0.5 * (1 + numerix.tanh(phi / epsilon))
density.setValue(rho2 + (rho1 - rho2) * (1 - H))
viscosity.setValue(mu2 + (mu1 - mu2) * (1 - H))

# -- Boundary Conditions --
# Manually construct the inlet velocity vector to avoid ValueError
left_faces = mesh.facesLeft
inlet_y = mesh.faceCenters[1][left_faces.value]
inlet_vx = inlet_velocity * (1.0 - (2.0 * inlet_y / (ny * dy) - 1.0)**2)
inlet_vy = np.zeros_like(inlet_vx)
inlet_value = np.array([inlet_vx, inlet_vy])

# Apply the constraints
velocity.constrain(inlet_value, where=mesh.facesLeft)
phi.constrain(-1.0, where=mesh.facesLeft)

# Outlet (right face)
pressure.constrain(0., where=mesh.facesRight)

# Walls (top and bottom)
# MODIFIED: Correctly apply no-slip condition for a colocated grid
velocity.constrain([0., 0.], where=mesh.facesTop | mesh.facesBottom)


# -- PDE Definitions --
# Level-Set Advection Equation
phi_eq = (TransientTerm(var=phi) + ConvectionTerm(coeff=velocity.faceValue, var=phi) == 0)

# Pressure-Poisson Equation
# MODIFIED: Use velocity.faceValue.divergence
pressure_eq = (ImplicitDiffusionTerm(var=pressure) ==
               (1. / dt) * velocity.faceValue.divergence)

# Reformulated Momentum Predictor Equation
velocity_star = CellVariable(name="u*", mesh=mesh, rank=1)

# MODIFIED: Correctly define the convection coefficient as a pure FaceVariable
convection_coeff = density.arithmeticFaceValue * velocity.faceValue

predictor_eq = (TransientTerm(coeff=density, var=velocity) ==
                -ConvectionTerm(coeff=convection_coeff, var=velocity)
                + ImplicitDiffusionTerm(coeff=viscosity, var=velocity))


# --- 5. Visualization Setup ---
print("Setting up visualization...")
fig, ax = plt.subplots(figsize=(10, 5))
cax = None # Handle for the colorbar

# --- 6. Main Solver Loop ---
print("Starting solver loop...")
for step in range(steps):
    print(f"Step {step+1}/{steps}")

    # Store old values
    phi.updateOld()
    pressure.updateOld()
    velocity.updateOld()

    # Update fluid properties based on the current phi
    H = 0.5 * (1 + numerix.tanh(phi / epsilon))
    density.setValue(rho2 + (rho1 - rho2) * (1 - H))
    viscosity.setValue(mu2 + (mu1 - mu2) * (1 - H))
    
    # --- Solve the System ---
    # 1. Predict tentative velocity (v*)
    # MODIFIED: Update the convection coefficient inside the loop
    convection_coeff.setValue(density.arithmeticFaceValue * velocity.faceValue)
    predictor_eq.solve(var=velocity, dt=dt)
    velocity_star.setValue(velocity)

    # 2. Solve for pressure
    # MODIFIED: Use velocity_star.faceValue.divergence
    pressure_eq.RHS.setValue((1. / dt) * velocity_star.faceValue.divergence)
    pressure_eq.solve(var=pressure, dt=dt)

    # 3. Correct the velocity
    # MODIFIED: Interpolate pressure gradient back to cell centers
    grad_p = pressure.grad
    velocity.setValue(velocity_star - (dt / density) * grad_p.arithmeticCellValue)

    # 4. Enforce no-slip on the cylinder obstacle
    velocity.setValue([0.,0.], where=cylinder_mask)

    # 5. Advect the interface
    phi_eq.solve(var=phi, dt=dt)

    # --- Update plot ---
    if __name__ == '__main__':
        if cax is not None:
            cax.remove()

        vel_mag_data = velocity.mag.value.reshape((ny, nx))
        phi_data = phi.value.reshape((ny, nx))

        ax.clear()
        im = ax.imshow(vel_mag_data, origin='lower',
                       extent=[0, nx*dx, 0, ny*dy],
                       cmap='viridis', vmin=0, vmax=inlet_velocity * 1.5)
        ax.contour(X, Y, phi_data, levels=[0], colors='k', linewidths=2)
        
        # Draw cylinder
        circle = plt.Circle(cylinder_center, cylinder_radius, color='white', ec='black')
        ax.add_patch(circle)
        
        cax = fig.colorbar(im, ax=ax)
        ax.set_title(f"Velocity Magnitude at t = {step*dt:.2f} s")
        ax.set_aspect('equal')
        plt.pause(0.01)

print("Simulation finished.")
plt.show(block=True)
