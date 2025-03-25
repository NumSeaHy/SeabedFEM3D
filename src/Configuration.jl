"""
This file contains the configuration parameters for run the simulation.
"""

using Revise

# Importing the necessary libraries for compute the physical properties of the porous domain
includet("SedimentModels.jl")
includet("BiotStollFuncs.jl")

using .SedimentModels
using .BiotStollFuncs

# Domain parametrization
L = 0.2 # Length of the domain [m]
t_F = 0.1 # Height of the domain [m]
t_P = 0.2 # Height of the porous domain [m]
w = 0.2
d_PML = 0.05 # Thickness of the PML [m]

# Scattering objects parametrization
    # Maximum number of iterations for the scenario generation
    max_iterations = 100000 # Maximum number of iterations to find a suitable position for the cockle
    
    # Spheres parametrization
    N_rigid_spheres = 0 # Number of rigid spheres in the porous domain [-]
    N_porous_spheres = 0 # Number of porous spheres in the porous domain [-]
    r = 5.0e-2 # Radius of the sphere [m]
    σ_r = 0.005 # Standard deviation of the  radius of the sphere [m]
    tol_sphere = 0.5 * r # Tolerance for the spheres to avoid collisions with the boundaries of the physical domain
    
    # Cockle parametrization
    N_open_cockles = 1 # Number of open cockles in the porous domain [-]
    N_closed_cockles = 0 # Number of closed cockles in the porous domain [-]
    by_default_cockle_radius =  6.1e-2/2 # Default radius of the cockle measured in the real geometry[m]
    σ_r_cockle = 1e-2 # Standard deviation of the radius of the cockle [m]

    # Queen scallop parametrization
    N_open_queenscallops = 0 # Number of open queen scallops in the porous domain [-]
    N_closed_queenscallops = 0 # Number of closed queen scallops in the porous domain [-]
    by_default_queenscallop_radius =  7.5e-2/2 # Default radius of the queen scallop measured in the real geometry[m]
    σ_r_queenscallop = 1e-2 # Standard deviation of the radius of the queen scallop [m]

    # Scallop parametrization
    N_open_scallops = 0 # Number of open scallops in the porous domain [-]
    N_closed_scallops = 0 # Number of closed scallops in the porous domain [-]
    by_default_scallop_radius =  10.7e-2/2 # Default radius of the scallop measured in the real geometry[m]
    σ_r_scallop = 1e-2 # Standard deviation of the radius of the scallop [m]
    

# Transducer pressure 
P0 = 5e5im

# Frequency parameters
f = 5e3
ω = 2 * π * f

# Fluid domains properties
ρ_F(ω) = 1000. + 0.0im # Mass density of the fluid [kg/m^3]
c_F(ω) = 1432. + 0.0im # Speed of sound in the fluid [m/s]
η_F = 1.0e-3 # Dynamic viscosity of the fluid [Pa s]


# Porous domain properties using the Biot-Stoll model
sediment(ω) = predefined_sediment("MediumSilt"; ρF=ρ_F(ω), KF=ρ_F(ω)*c_F(ω)^2, η=η_F) # Check list_sediment(ω)s() for more options of predefined sediment(ω)s
ρ_P(ω) = sediment(ω).β * sediment(ω).ρF + (1 - sediment(ω).β) * sediment(ω).ρr # Mass density of the porous domain [kg/m^3]
c_P(ω) = compute_wave_properties(ω, sediment(ω))[1] + 1im*compute_wave_properties(ω, sediment(ω))[2]/ω*compute_wave_properties(ω, sediment(ω))[1]^2 # [1] returns the real part of the phase velocity and [2] the attenuation coefficient α

# Wavenumbers associated with the fluid and porous domains
k_F(ω) = ω/c_F(ω)
k_P(ω) = ω/c_P(ω)

# Define the bulk modulus associated with the fluid and porous domains
K_F(ω) = ρ_F(ω) * c_F(ω)^2
K_P(ω) = ρ_P(ω) * c_P(ω)^2

# PML parameters
RPML = 1e-5 # Reflection coefficient of the PML
σ_0 = -3/4 * log(RPML)/d_PML
