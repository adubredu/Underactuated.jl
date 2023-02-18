module Underactuated

using MuJoCo

include("types.jl")
include("step.jl")

include("robots/simple_pendulum.jl")

# Robot types
export SimplePendulumType 

# Robots 
export SimplePendulum

export create_robot,
       step
end
