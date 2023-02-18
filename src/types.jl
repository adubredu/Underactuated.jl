# robot types 
struct SimplePendulumType
end

# robots
abstract type Robot end

mutable struct SimplePendulum <: Robot
    model::Any 
    data::Any 
    viewer::Any
    visualize::Bool 
end
