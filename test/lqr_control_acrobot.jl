using Revise
using Underactuated 
using Underactuated.LinearAlgebra

robot = create_robot(AcrobotType())

A, B = get_linearized_dynamics(robot)
Q, R = I, I
Horizon = 5000
robot.data.qpos[0] = 1e-3#π
robot.data.qpos[1] = 1e-3#2π

step(robot)

for i = 1:Horizon 
    lqr_controller(A, B, Q, R, robot) 
    step(robot) 
    robot.viewer.render()
end

robot.viewer.close()