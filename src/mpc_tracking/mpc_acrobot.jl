function initialize_solver(robot::Acrobot; horizon=1) 
    model = JuMP.Model(OSQP.Optimizer)
    set_silent(model)
    @variable(model, x[1:robot.Nx, 1:horizon])
    @variable(model, u[1:robot.Nu, 1:horizon])
    return model 
end

function linearize_dynamics_along_trajectory(xtraj, utraj, robot::Acrobot)
    As = []; Bs = []
    _, Ntraj = size(xtraj)
    for i=1:Ntraj 
        x = xtraj[:, i]
        u = utraj[:, i]
        A, B = get_linearized_dynamics(x, u, robot)
        push!(As, A)
        push!(Bs, B)
    end
    return As, Bs
end

function solve_mpc(xtraj, utraj, robot_model, As, Bs, i::Int, solver, Nsteps; horizon=5)
    x0 = get_generalized_coordinates(robot_model)
    x = solver.obj_dict[:x]
    u = solver.obj_dict[:u]
    try delete(solver, solver.ext[:init_state])  catch nothing; end 
    solver.ext[:init_state] = @constraint(solver, x[:,1] .== x0)
    N = horizon+i > Nsteps ? Nsteps-i : horizon 
    Q = Diagonal(100.0*ones(4)); R = 100
    @objective(solver, Min, sum([0.5*((x[:,k]-xtraj[:,k+i-1])'*Q*(x[:,k]-xtraj[:, k+i-1])) + 0.5*((u[:,k]-utraj[:,k+i-1])'*R*(u[:,k]-utraj[:,k+i-1])) for k=1:N ]))
    try delete(solver, solver.ext[:dynamics])  catch nothing; end 
    solver.ext[:dynamics] = @constraint(solver, [k = 1:N-1],  x[:, k+1] .== As[k+i-1]*x[:, k] + Bs[k+i-1]*u[:, k])
    optimize!(solver)
    τs = value.(u)
    τ = -τs[1]
    return τ
end