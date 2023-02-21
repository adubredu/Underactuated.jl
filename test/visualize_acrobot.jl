using Revise
using Underactuated 

robot = create_robot(AcrobotType())
robot.data.qpos[0] = π
robot.data.qpos[1] = 2π
Horizon = 500

for _ = 1:Horizon
    step(robot)
    @show get_generalized_coordinates(robot)
    robot.viewer.render()
end

robot.viewer.close()