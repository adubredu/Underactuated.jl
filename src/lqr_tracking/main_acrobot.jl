using Revise 
using StaticArrays
using Underactuated
using LinearAlgebra
using ControlSystems

robot_model = create_robot(AcrobotType())
As, Bs = linearize_dynamics_along_trajectory(xtraj, utraj, robot_model)
As, Bs = linearize_dynamics_along_trajectory(xtraj, utraj, robot_model)
