import Base: step 
function step(robot::Robot)
    mujoco.mj_step(robot.model, robot.data)
end