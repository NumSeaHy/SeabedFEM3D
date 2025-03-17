"""
This file contains the configuration parameters for run the simulation.
"""

# Importing the necessary libraries for compute the physical properties of the porous domain
include("SedimentModels.jl")
include("BiotStollFuncs.jl")


using .SedimentModels
using .BiotStollFuncs

# Domain parametrization
L = 0.2 # Length of the domain [m]
t_F = 0.1 # Height of the domain [m]
t_P = 0.2 # Height of the porous domain [m]
w = 0.2
d_PML = 0.05 # Thickness of the PML [m]

# Scattering objects parametrization
    
    # Spheres parametrization
    N_rigid_spheres = 0 # Number of rigid spheres in the porous domain [-]
    N_porous_spheres = 0 # Number of porous spheres in the porous domain [-]
    r = 5.0e-2 # Radius of the sphere [m]
    σ_r = 0.005 # Standard deviation of the  radius of the sphere [m]
    tol_sphere = 0.5 * r # Tolerance for the spheres to avoid collisions with the boundaries of the physical domain
    
    # Cockle parametrization
    N_cockles = 2 # Number of cockles in the porous domain [-]
    by_default_radius =  0.075/2 # Default radius of the cockle measured in the real geometry[m]
    σ_r_cockle = 0.01 # Standard deviation of the radius of the cockle [m]
    cockle_brep_path = "./cockle_geometries/Cockle.brep"
    cockle_closed_brep_path = "./cockle_geometries/ClosedCockle.brep"
    max_iterations = 50000 # Maximum number of iterations to find a suitable position for the cockle
    

# Transducer pressure 
P_0 = 5e5im

# Frequency parameters
f = 15e3
ω = 2 * π * f

# Fluid domains properties
ρ_F(ω) = 4000.
c_F(ω) = 2000.
η_F = 1.0e-3 # Dynamic viscosity of the fluid [Pa s]


# Porous domain properties using the Biot-Stoll model
sediment(ω) = predefined_sediment("MediumSilt"; ρF=ρ_F(ω), KF=ρ_F(ω)*c_F(ω)^2, η=η_F) # Check list_sediment(ω)s() for more options of predefined sediment(ω)s
ρ_P(ω) = sediment(ω).β * sediment(ω).ρF + (1 - sediment(ω).β) * sediment(ω).ρr # Mass density of the porous domain [kg/m^3]
c_P(ω) = compute_wave_properties(ω, sediment(ω))[1] + 1im*compute_wave_properties(ω, sediment(ω))[2]/ω*compute_wave_properties(ω, sediment(ω))[1]^2 # [1] returns the real part of the phase velocity and [2] the attenuation coefficient α

# Wavenumbers associated with the fluid and porous domains
k_F(ω) = ω/c_F(ω)
k_P(ω) = ω/c_P(ω)

# PML parameters
RPML = 1e-5 # Reflection coefficient of the PML
σ_0 = -3/4 * log(RPML)/d_PML
