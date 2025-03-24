module AnalyticalIncident

export analytical_incident

using Gridap
using Gridap.Fields

function analytical_incident(x, ω, P0)
    u_x = 1 * P0 * x[1] * ω^2
    u_y = 2 * P0 * x[2] * ω^2
    u_z = 3 * P0 * x[3] * ω^2

    return VectorValue(u_x, u_y, u_z)

end


end