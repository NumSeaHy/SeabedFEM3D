# Warning! Once the code is functional, is recommendable to define all the parameters with a 
# const value (i.e const a = 3), because this will optimize the use of Julia JIT

# Domain parametrization
L = 1 # Length of the domain [m]
H = 0.4 # Height of the domain [m]
w = 1
d_PML = 0.1 # Thickness of the PML [m]
r = 0.1 # Radius of the rock [m]
x₀ = 0 # x-coordinate of the center of the source [m]
y₀ = H / 2 # y-coordinate of the center of the source [m]
z₀ = 0 # z-coordinate of the center of the source [m]

# Frequency parameters
f = 15e3
ω = 2 * π * f

# Domains properties
ρ = 2000.
c = 8000.
# Transducer pressure
P_0 = 5e5im

RPML = 1e-5 # Reflection coefficient of the PML
σ_0 = -3/4 * log(RPML)/d_PML