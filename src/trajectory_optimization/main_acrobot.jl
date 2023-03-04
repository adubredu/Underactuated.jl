using Revise
using TrajOptPlots
using MeshCat
using StaticArrays
using JuMP
using LinearAlgebra
import Ipopt 
using Plots 
using RobotZoo
const RZ = RobotZoo 
using RobotDynamics
using Underactuated

include("dircol_acrobot.jl")
include("visualize_acrobot.jl")

# initial and goal states 
x0 = [-pi/2; 0; 0; 0]
xgoal = [pi/2; 0; 0; 0]
T = 5.0 # seconds
robot = RobotZoo.Acrobot()


xtraj, utraj, Thist = compute_trajectory(x0, xgoal, robot, T)
visualize_acrobot_trajectory(robot, xtraj, utraj, Thist)