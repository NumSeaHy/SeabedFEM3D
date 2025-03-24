using Gridap
using GridapGmsh
using Revise

include("Configuration.jl")
includet("AnalyticalIncident.jl")
using .AnalyticalIncident


model = GmshDiscreteModel("./data/mesh.msh")

Ω = Triangulation(model)

xp = get_physical_coordinate(Ω)

u(x) = analytical_incident(x, ω, P0)

u_field = u ∘ xp


writevtk(Ω,"./results/result_interpolated.vtu", cellfields=[ "Re(u)"=>real(u_field), "Im(u)"=>imag(u_field)])