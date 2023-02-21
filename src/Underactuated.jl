module Underactuated

using MuJoCo
using MuJoCo.PythonCall
using LinearAlgebra
using ControlSystems 
using JuMP 
using OSQP

include("types.jl")
include("step.jl")

include("robots/simple_pendulum.jl")
include("robots/acrobot.jl")

include("controllers/energy_shaping.jl")
include("controllers/lqr.jl")
include("controllers/convex_mpc.jl")

# Robot types
export SimplePendulumType,
        AcrobotType

# Robots 
export SimplePendulum,
        Acrobot

export create_robot,
       step,
       get_generalized_coordinates,
       compute_energy,
       apply_torque!,
       get_linearized_dynamics

# Controllers 
export energy_shaping_controller,
        lqr_controller,
        initialize_solver, 
        build_QP,
        convex_mpc_controller

end
