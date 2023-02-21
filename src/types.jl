# robot types 
struct SimplePendulumType end
struct AcrobotType end

# robots
abstract type Robot end

mutable struct SimplePendulum <: Robot
    model::Any 
    data::Any 
    viewer::Any
    visualize::Bool 
    mass::Float64 
    length::Float64 
    gravity::Float64 
    function SimplePendulum(model, data, viewer, visualize)
        mass = 1.0
        length = 0.5
        gravity = 9.81
        new(model, data, viewer, visualize, mass, length, gravity)
    end
end

mutable struct Acrobot <: Robot
    model::Any 
    data::Any 
    viewer::Any
    visualize::Bool 
    xd::Vector{Float64} 
    Nx::Int 
    Nu::Int 
    τ_min::Float64 
    τ_max::Float64
    function Acrobot(model, data, viewer, visualize)
        xd=zeros(4) 
        Nx = 4; Nu = 1; τ_min = -5e3; τ_max = 5e3
        new(model, data, viewer, visualize, xd, Nx, Nu, τ_min, τ_max)
    end
end