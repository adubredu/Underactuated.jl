using Revise
using Underactuated 

robot = create_robot(SimplePendulumType())

Horizon = 10000

for _ = 1:Horizon
    step(robot)
    @show get_generalized_coordinates(robot)
    robot.viewer.render()
end

robot.viewer.close()