using Revise
using Underactuated 
using Underactuated.LinearAlgebra

T = 5000
robot = create_robot(AcrobotType())

A, B = get_linearized_dynamics(robot)
Q, R = 1e3*I(robot.Nx), 0.01*I(robot.Nu)
horizon = 100
robot.data.qpos[0] = 1e-1#π
robot.data.qpos[1] = 1e-1#2π

solver = initialize_solver(robot; horizon=horizon) 
solver = build_QP(A, B, Q, R, solver, robot; horizon=horizon)

step(robot)
# convex_mpc_controller(solver, robot)
for i = 1:T 
    convex_mpc_controller(solver, robot)
    step(robot) 
    robot.viewer.render()
end

robot.viewer.close()