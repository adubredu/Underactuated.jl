function create_robot(type::SimplePendulumType; visualize=true)
    xml_path = joinpath(@__DIR__, "../models/simple_pendulum.xml")
    model = mujoco.MjModel.from_xml_path(xml_path)
    data = mujoco.MjData(model) 
    viewer = visualize ? mujoco_viewer.MujocoViewer(model, data) : nothing
    robot = SimplePendulum(model, data, viewer, visualize)

    return robot
end 