module AnalyticalIncident

export analytical_incident, spherical_incident_field

using Gridap
using Gridap.Fields
using SpecialFunctions

function distance(x, x₀, y₀, z₀)
    return sqrt((x[1]-x₀)^2 + (x[2]-y₀)^2 + (x[3]-z₀)^2)
end

function analytical_incident(x, ω, P0)
    u_x = P0 + 0.0im
    u_y = P0 + 0.0im
    u_z = P0 + 0.0im

    return VectorValue(u_x, u_y, u_z)
end

function spherical_incident_field(x, x₀, y₀, z₀, r, kF, ρF, P0, ω)
    
    C = P0 * kF * 1/(2 * ρF * ω^2 * hankelh1(1/2, kF*r) * sqrt(pi/(2*kF*r)))
    
    ux = C * sqrt(pi/(2*kF*distance(x, x₀, y₀, z₀))) * (hankelh1(-1/2, kF*distance(x, x₀, y₀, z₀)) -
                                                        hankelh1(1/2, kF*distance(x, x₀, y₀, z₀))/(kF*distance(x, x₀, y₀, z₀)) -
                                                        hankelh1(3/2, kF*distance(x, x₀, y₀, z₀))) * 
                                                        (x[1]-x₀)/distance(x, x₀, y₀, z₀)

    uy = C * sqrt(pi/(2*kF*distance(x, x₀, y₀, z₀))) * (hankelh1(-1/2, kF*distance(x, x₀, y₀, z₀)) -
                                                        hankelh1(1/2, kF*distance(x, x₀, y₀, z₀))/(kF*distance(x, x₀, y₀, z₀)) -
                                                        hankelh1(3/2, kF*distance(x, x₀, y₀, z₀))) * 
                                                        (x[2]-y₀)/distance(x, x₀, y₀, z₀)
    

    uz = C * sqrt(pi/(2*kF*distance(x, x₀, y₀, z₀))) * (hankelh1(-1/2, kF*distance(x, x₀, y₀, z₀)) -
                                                        hankelh1(1/2, kF*distance(x, x₀, y₀, z₀))/(kF*distance(x, x₀, y₀, z₀)) -
                                                        hankelh1(3/2, kF*distance(x, x₀, y₀, z₀))) * 
                                                        (x[3]-z₀)/distance(x, x₀, y₀, z₀)
                                                        

    return VectorValue(ux, uy, uz)
    
end


end