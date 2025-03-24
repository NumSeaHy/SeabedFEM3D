module AnalyticalIncident

using Gridap
using Gridap.Fields

function analytical_incident(x, ω, P0)
    ux = 1 * x[1] * ω^2
    uy = 2 * x[2] * ω^2
    uz = 3 * x[3] * ω^2

    return VectorValue(u_x, u_y, u_z)

end


end