using Revise
using Underactuated 

robot = create_robot(SimplePendulumType())

Horizon = 5000

for _ = 1:Horizon
    step(robot)
    robot.viewer.render()
end

robot.viewer.close()