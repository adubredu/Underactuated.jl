function create_robot(type::SimplePendulumType; visualize=true)
    xml_path = joinpath(@__DIR__, "../models/simple_pendulum.xml")
    model = mujoco.MjModel.from_xml_path(xml_path)
    data = mujoco.MjData(model) 
    viewer = visualize ? mujoco_viewer.MujocoViewer(model, data) : nothing
    robot = SimplePendulum(model, data, viewer, visualize)

    return robot
end 

function get_generalized_coordinates(robot::SimplePendulum)
    data = robot.data 
    θ = pyconvert(Float64, data.joint("pin").qpos[0])
    θ̇  = pyconvert(Float64, data.joint("pin").qvel[0])
    return [θ, θ̇ ]
end

function apply_torque!(τ::Float64,  robot::SimplePendulum)   
    robot.data.actuator("pin_torque").ctrl[0] = τ
end