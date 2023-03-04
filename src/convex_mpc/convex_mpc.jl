function initialize_solver(robot::Acrobot; horizon=1) 
    model = JuMP.Model(OSQP.Optimizer)
    set_silent(model)
    @variable(model, x[1:robot.Nx, 1:horizon])
    @variable(model, u[1:robot.Nu, 1:horizon])
    return model 
end

function build_QP(A, B, Q, R, model, robot::Acrobot; horizon=1) 
    x = model.obj_dict[:x]
    u = model.obj_dict[:u]
    A_h = [A for _=1:horizon]
    B_h = [B for _=1:horizon]
    Q_h = [Q for _=1:horizon]
    R_h = [R for _=1:horizon] 
    P = are(Discrete, A, B, Matrix{Float64}(Q), Matrix{Float64}(R))
    # display(P)
    P_h = [I(robot.Nx) for _=1:horizon]
    # P_h = [P for _=1:horizon]
    xd = robot.xd
    τ_mins = [robot.τ_min for _=1:horizon]
    τ_maxs = [robot.τ_max for _=1:horizon]
    @objective(model, Min,  sum([0.5*((x[:, i]-xd)'*Q_h[i]*(x[:, i]-xd)) + 0.5*(u[:, i]'*R_h[i]*u[:, i])  + (x[:, i]-xd)'*P_h[i]*(x[:, i]-xd) for i=1:horizon])) #+ (x[:, i]-xd)'*P_h[i]*(x[:, i]-xd)
    @constraint(model, [i = 1:horizon-1],  x[:, i+1] .== A_h[i]*x[:, i] + B_h[i]*u[:, i])
    @constraint(model, [i = 1:horizon],  τ_mins[i] .≤ u[:,i] .≤ τ_maxs[i])
    return model
end

function convex_mpc_controller(model, robot::Acrobot) 
    x = model.obj_dict[:x]
    u = model.obj_dict[:u]
    xc = get_generalized_coordinates(robot)
    try delete(model, model.ext[:init_state])  catch nothing; end 
    model.ext[:init_state] = @constraint(model, x[:, 1] .== xc)
    optimize!(model)
    τs = value.(u)
    τ = τs[1]  
    robot.data.ctrl[0] = -τ
end