function lqr_controller(A, B, Q, R, robot::Acrobot) 
    K = lqr(Discrete, A, B, Matrix{Float64}(Q), Matrix{Float64}(R)) 
    xc = get_generalized_coordinates(robot)
    xd = robot.xd 
    Δx = xd - xc 
    u = -K*Δx
    τ = apply_torque_limits(robot, u[1])
    robot.data.ctrl[0] = τ
end