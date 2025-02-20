# Warning! Once the code is functional, is recommendable to define all the parameters with a 
# const value (i.e const a = 3), because this will optimize the use of Julia JIT

# Domain parametrization
L = 1 # Length of the domain [m]
t_F = 0.4 # Height of the domain [m]
t_P = 0.8 # Height of the porous domain [m]
w = 1
d_PML = 0.1 # Thickness of the PML [m]

# Frequency parameters
f = 15e3
ω = 2 * π * f

# Fluid domains properties
ρ_F(ω) = 2000.
c_F(ω) = 8000.

# Porous domain properties
ρ_P = 1000.
c_P(ω) = 4000.

# Wavenumbers associated with the fluid and porous domains
k_F(ω) = ω/c_F(ω)
k_P(ω) = ω/c_P(ω)

# PML parameters
RPML = 1e-5 # Reflection coefficient of the PML
σ_0 = -3/4 * log(RPML)/d_PML
