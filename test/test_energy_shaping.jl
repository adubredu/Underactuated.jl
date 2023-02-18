using Revise
using Underactuated 

robot = create_robot(SimplePendulumType())

Horizon = 5000
ld = robot.length*cos(1.57)
step(robot)

for _ = 1:Horizon
    τ = energy_shaping_controller(ld, robot)
    apply_torque!(τ,  robot)
    step(robot)
    robot.viewer.render()
end

robot.viewer.close()