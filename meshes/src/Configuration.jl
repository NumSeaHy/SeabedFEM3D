# Warning! Once the code is functional, is recommendable to define all the parameters with a 
# const value (i.e const a = 3), because this will optimize the use of Julia JIT

# Domain parametrization
L = 1 # Length of the domain [m]
t_F = 0.4 # Height of the domain [m]
t_P = 0.8 # Height of the porous domain [m]
w = 1
d_PML = 0.1 # Thickness of the PML [m]

# Scattering objects parametrization
    
    # Spheres parametrization
    N_rigid_spheres = 0 # Number of rigid spheres in the porous domain [-]
    N_porous_spheres = 0 # Number of porous spheres in the porous domain [-]
    r = 5.0e-2 # Radius of the sphere [m]
    σ_r = 0.005 # Standard deviation of the  radius of the sphere [m]
    tol_sphere = 0.5 * r # Tolerance for the spheres to avoid collisions with the boundaries of the physical domain
    
    # Cockle parametrization
    N_cockles = 1 # Number of cockles in the porous domain [-]
    by_default_radius =  0.075/2 # Default radius of the cockle measured in the real geometry[m]
    σ_r_cockle = 0.05 # Standard deviation of the radius of the cockle [m]
    cockle_brep_path = "./clam_geometries/EntireClam.brep"
    max_iterations = 1000 # Maximum number of iterations to find a suitable position for the cockle
    

# Transducer pressure 
P_0 = 5e5im

# Frequency parameters
f = 5e3
ω = 2 * π * f

# Fluid domains properties
ρ_F(ω) = 4000.
c_F(ω) = 2000.

# Porous domain properties
ρ_P = 1000.
c_P(ω) = 3000.

# Wavenumbers associated with the fluid and porous domains
k_F(ω) = ω/c_F(ω)
k_P(ω) = ω/c_P(ω)

# PML parameters
RPML = 1e-5 # Reflection coefficient of the PML
σ_0 = -3/4 * log(RPML)/d_PML
