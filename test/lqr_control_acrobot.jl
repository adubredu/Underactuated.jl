using Revise
using Underactuated 
using Underactuated.LinearAlgebra

robot = create_robot(AcrobotType())

A, B = get_linearized_dynamics(robot)
Q, R = 1e3*I(robot.Nx), 0.01*I(robot.Nu)
T = 2000
robot.data.qpos[0] = 1e-1#π
robot.data.qpos[1] = 1e-1#2π

step(robot)

for i = 1:T 
    lqr_controller(A, B, Q, R, robot) 
    step(robot) 
    robot.viewer.render()
end

robot.viewer.close()