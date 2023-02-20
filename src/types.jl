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
    m1::Float64 
    m2::Float64 
    l1::Float64 
    l2::Float64 
    g::Float64 
    function Acrobot(model, data, viewer, visualize)
        m1 = 1
        m2 = 1 
        l1 = 0.05
        l2 = 0.05
        g = 9.81
        new(model, data, viewer, visualize, m1, m2, l1, l2, g)
    end
end