function create_robot(type::AcrobotType; visualize=true)
    xml_path = joinpath(@__DIR__, "../models/acrobot.xml")
    model = mujoco.MjModel.from_xml_path(xml_path)
    data = mujoco.MjData(model) 
    viewer = visualize ? mujoco_viewer.MujocoViewer(model, data) : nothing
    robot = Acrobot(model, data, viewer, visualize)

    return robot
end 

function get_generalized_coordinates(robot::Acrobot)
    data = robot.data 
    θ1 = pyconvert(Float64, data.joint("shoulder").qpos[0])
    θ2 = pyconvert(Float64, data.joint("elbow").qpos[0])
    θ̇1  = pyconvert(Float64, data.joint("shoulder").qvel[0])
    θ̇2  = pyconvert(Float64, data.joint("elbow").qvel[0])
    return [θ1, θ2, θ̇1, θ̇2 ]
end

function get_linearized_dynamics(robot::Acrobot)
    θ₀ = [0.0, 0.0] 
    robot.xd = [θ₀; zeros(2)] 
    robot.data.qpos[0] = θ₀[1]
    robot.data.qpos[1] = θ₀[2]
    robot.data.qacc = 0
    mujoco.mj_inverse(robot.model, robot.data)
    qfrc0 = robot.data.qfrc_inverse
    ctrl0 = pyconvert(Vector{Float64}, qfrc0)' * 
                pinv(pyconvert(Matrix{Float64}, robot.data.actuator_moment))
    ctrl0 = ctrl0[1]
    robot.data.ctrl[0] = ctrl0  
    nv = pyconvert(Int, robot.model.nv); nu = pyconvert(Int, robot.model.nu)
    A_np = numpy.zeros((2*nv, 2*nv))
    B_np = numpy.zeros((2*nv, nu))
    ϵ = 1e-6
    centered = true
    mujoco.mjd_transitionFD(robot.model, robot.data, ϵ, centered, A_np, B_np, nothing, nothing)
    A = pyconvert(Matrix{Float64}, A_np)
    B = pyconvert(Matrix{Float64}, B_np)
    return A, B     
end

function get_linearized_dynamics(x, u, robot::Acrobot) 
    robot.xd = x 
    robot.data.qpos[0] = x[1]
    robot.data.qpos[1] = x[2]
    robot.data.qvel[0] = x[3]
    robot.data.qvel[1] = x[4]
    robot.data.qacc = 0
    mujoco.mj_inverse(robot.model, robot.data)
    qfrc0 = robot.data.qfrc_inverse
    ctrl0 = pyconvert(Vector{Float64}, qfrc0)' * 
                pinv(pyconvert(Matrix{Float64}, robot.data.actuator_moment))
    ctrl0 = ctrl0[1]
    robot.data.ctrl[0] = ctrl0  
    nv = pyconvert(Int, robot.model.nv); nu = pyconvert(Int, robot.model.nu)
    A_np = numpy.zeros((2*nv, 2*nv))
    B_np = numpy.zeros((2*nv, nu))
    ϵ = 1e-6
    centered = true
    mujoco.mjd_transitionFD(robot.model, robot.data, ϵ, centered, A_np, B_np, nothing, nothing)
    A = pyconvert(Matrix{Float64}, A_np)
    B = pyconvert(Matrix{Float64}, B_np)
    return A, B     
end

function apply_torque_limits(robot::Acrobot, u::Float64)
    return clamp(u, robot.τ_min, robot.τ_max)
end

function apply_torque!(robot::Acrobot, τ::Float64) 
    robot.data.actuator("elbow").ctrl[0] = -τ
end