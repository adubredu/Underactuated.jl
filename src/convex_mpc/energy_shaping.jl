function compute_energy(l, robot::SimplePendulum)
    m = robot.mass; g = robot.gravity 
    θ, θ̇ = get_generalized_coordinates(robot)
    E = 0.5*m*l^2*θ̇  - m*g*l*cos(θ)
    return E 
end

function energy_shaping_controller(l, robot::SimplePendulum)
    K = 1.0
    q, q̇ = get_generalized_coordinates(robot)
    E = compute_energy(robot.length, robot)
    Eᵈ = compute_energy(l, robot)
    Ē = E - Eᵈ
    u = -K*q̇^2*Ē
    return u 
end


