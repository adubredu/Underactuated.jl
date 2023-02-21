function lqr_controller(A, B, Q, R, robot::Acrobot) 
    K = lqr(A, B, Q, R) 
    xc = get_generalized_coordinates(robot)
    xd = robot.xd 
    Δx = xd - xc 
    u = -K*Δx
    τ = apply_torque_limits(robot, u[1])
    robot.data.ctrl[0] = τ
end