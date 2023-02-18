module Underactuated

using MuJoCo
using MuJoCo.PythonCall

include("types.jl")
include("step.jl")

include("robots/simple_pendulum.jl")
include("controllers/energy_shaping.jl")

# Robot types
export SimplePendulumType 

# Robots 
export SimplePendulum

export create_robot,
       step,
       get_generalized_coordinates,
       compute_energy,
       apply_torque!

# Controllers 
export energy_shaping_controller

end
