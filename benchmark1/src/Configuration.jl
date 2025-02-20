# Warning! Once the code is functional, is recommendable to define all the parameters with a 
# const value (i.e const a = 3), because this will optimize the use of Julia JIT

# Domain parametrization
L = 1 # Length of the domain [m]
H = 0.8 # Height of the domain [m]
w = 1 # Width of the domain [m]
d_PML = 0.1 # Thickness of the PML [m]

# Sphere parametrization
x₀ = 0.0 # x-coordinate of the center of the sphere [m]
y₀ = H/2 # y-coordinate of the center of the sphere [m]
z₀ = 0.0 # z-coordinate of the center of the sphere [m]
r = 0.2 # Radius of the sphere [m]

# Transducer pressure
P_0 = 5e5im

# Frequency parameters
f = 15e3
ω = 2 * π * f

# Fluid domains properties
ρ_F(ω) = 2000.
c_F(ω) = 8000.

# Wavenumbers associated with the fluid and porous domains
k_F(ω) = ω/c_F(ω)

# PML parameters
RPML = 1e-5 # Reflection coefficient of the PML
σ_0 = -3/4 * log(RPML)/d_PML
