using Revise 
using StaticArrays
using JuMP 
using OSQP 
using Underactuated
using LinearAlgebra
import Ipopt 
using Plots 
using RobotZoo
const RZ = RobotZoo 
using RobotDynamics 

include("mpc_acrobot.jl")
include("../trajectory_optimization/dircol_acrobot.jl")

# initial and goal states 
x0 = [-pi/2; 0; 0; 0]
xgoal = [pi/2; 0; 0; 0]
T = 10.0 # seconds
h = 0.1 # timestep 1KHz
horizon=10
robot_dynamics = RobotZoo.Acrobot() 
# xtraj, utraj, Thist, Nsteps = compute_trajectory(x0, xgoal, robot_dynamics, T;h=h)

robot_model = create_robot(AcrobotType())
mpc_solver = initialize_solver(robot_model; horizon=horizon)
As, Bs = linearize_dynamics_along_trajectory(xtraj, utraj, robot_model)

for i=1:Nsteps 
    τ = solve_mpc(xtraj, utraj, robot_model, As, Bs, i, mpc_solver, Nsteps; horizon=horizon)
    apply_torque!(robot_model, τ)
    step(robot_model) 
    robot_model.viewer.render()
    @show i
end

robot_model.viewer.close()