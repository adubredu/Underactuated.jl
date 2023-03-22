using LinearAlgebra
using Plots
using ForwardDiff
using Ipopt
using MathOptInterface
const MOI = MathOptInterface;


Nx = 8     # number of state
Nu = 2     # number of controls
Tfinal = 4.4 # final time
h = 0.1 #10 hz
Nm = 5 #number of steps in each mode
Nt = Int(ceil(Tfinal/h)+1)   # number of time steps
Nmodes = Int(ceil(Nt/Nm))
thist = Array(range(0,h*(Nt-1), step=h));
n_nlp = Nx*Nt + Nu*(Nt-1) # number of decision variables
c_init_inds = 1:Nx
c_term_inds = (c_init_inds[end]+1):(c_init_inds[end]+Nx)
c_dyn_inds = (c_term_inds[end]+1):(c_term_inds[end]+Nx*(Nt-1))
c_stance_inds = (c_dyn_inds[end]+1):(c_dyn_inds[end]+Int(ceil(Nmodes/2)*Nm))
c_length_inds = (c_stance_inds[end]+1):(c_stance_inds[end]+Nt)
m_nlp = c_length_inds[end]

#Hopper Dynamics
g = 9.81
m1 = 5.0 #body mass
m2 = 1.0 #foot mass
ℓ_min = 0.5 #minimum length
ℓ_max = 1.5 #maximum length

function flight_dynamics(x,u)
    M = Diagonal([m1 m1 m2 m2])
    
    r1 = x[1:2]
    r2 = x[3:4]
    v = x[5:8]
    
    ℓ1 = (r1[1]-r2[1])/norm(r1-r2)
    ℓ2 = (r1[2]-r2[2])/norm(r1-r2)
    B = [ℓ1  ℓ2;
         ℓ2 -ℓ1;
        -ℓ1 -ℓ2;
        -ℓ2  ℓ1]
    
    v̇ = [0; -g; 0; -g] + M\(B*u)
    
    ẋ = [v; v̇]
end

function stance_dynamics(x,u)
    m1 = 5.0 #body mass
    m2 = 1.0 #foot mass
    M = Diagonal([m1 m1 m2 m2])
    g = 9.81
    
    r1 = x[1:2]
    r2 = x[3:4]
    v = x[5:8]
    
    ℓ1 = (r1[1]-r2[1])/norm(r1-r2)
    ℓ2 = (r1[2]-r2[2])/norm(r1-r2)
    B = [ℓ1  ℓ2;
         ℓ2 -ℓ1;
         0   0;
         0   0]
    
    v̇ = [0; -g; 0; 0] + M\(B*u)
    
    ẋ = [v; v̇]
end

function flight_dynamics_rk4(x,u)
    #RK4 integration with zero-order hold on u
    f1 = flight_dynamics(x, u)
    f2 = flight_dynamics(x + 0.5*h*f1, u)
    f3 = flight_dynamics(x + 0.5*h*f2, u)
    f4 = flight_dynamics(x + h*f3, u)
    return x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
end

function stance_dynamics_rk4(x,u)
    #RK4 integration with zero-order hold on u
    f1 = stance_dynamics(x, u)
    f2 = stance_dynamics(x + 0.5*h*f1, u)
    f3 = stance_dynamics(x + 0.5*h*f2, u)
    f4 = stance_dynamics(x + h*f3, u)
    return x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
end

function jump_map(x)
    #Assume the foot experiences inelastic collisions
    xn = [x[1:6]; 0.0; 0.0]
    return xn
end


# Cost weights
Q = Diagonal([1.0*ones(4); 1.0*ones(4)]);
R = 0.001;
Qn = Q


function stage_cost(x,u,k)
    return 0.5*((x-xref[:,k])'*Q*(x-xref[:,k])) + 0.5*(u-uref)'*R*(u-uref)
end
function terminal_cost(x)
    return 0.5*((x-xref[:,end])'*Qn*(x-xref[:,end]))
end


function cost(ztraj)
    z = reshape(ztraj[1:(end-Nx)],Nx+Nu,Nt-1)
    xtraj = [z[1:Nx,:] ztraj[end-(Nx-1):end]]
    utraj = z[(Nx+1):(Nx+Nu),:]
    J = 0.0
    for k = 1:(Nt-1)
        J += stage_cost(xtraj[:,k],utraj[:,k],k)
    end
    J += terminal_cost(xtraj[:,end])
    return J
end
     

#Reference Trajectory
uref = [m1*g; 0.0]
xref = zeros(Nx,Nt)
xref[1,:] .= LinRange(-1.0,1.0,Nt)
xref[2,:] .= 1.0 .+ 0.5.*sin.(2*pi/10.0*(0:(Nt-1)));
xref[3,:] .= LinRange(-1.0,1.0,Nt)
xref[5,2:end-1] .= (2.0/Tfinal)*ones(Nt-2)
xref[7,2:end-1] .= (2.0/Tfinal)*ones(Nt-2);

# stance_indices = []; flight_indices = []; jump_indices = [];
stance_const_indices = [];

function dynamics_constraint!(c,ztraj)
    d = reshape(view(c,c_dyn_inds),Nx,Nt-1)
    z = reshape(ztraj[1:(end-Nx)],Nx+Nu,Nt-1)
    xtraj = [z[1:Nx,:] ztraj[end-(Nx-1):end]]
    utraj = z[(Nx+1):(Nx+Nu),:]
    for k = 1:(Nmodes-1)
        if mod(k,2) == 1
            for j = 1:Nm
                s = (k-1)*Nm + j
                d[:,s] = stance_dynamics_rk4(xtraj[:,s],utraj[:,s]) - xtraj[:,s+1]
                # push!(stance_indices, s)
            end
        else
            for j = 1:(Nm-1)
                s = (k-1)*Nm + j
                d[:,s] = flight_dynamics_rk4(xtraj[:,s],utraj[:,s]) - xtraj[:,s+1]
                #push!(flight_indices, s)
            end
            s = k*Nm
            d[:,s] = jump_map(flight_dynamics_rk4(xtraj[:,s],utraj[:,s])) - xtraj[:,s+1]
            #push!(jump_indices, s)
        end
    end
    if mod(Nmodes,2) == 1
        for j = 1:(Nm-1)
            s = (Nmodes-1)*Nm + j
            d[:,s] = stance_dynamics_rk4(xtraj[:,s],utraj[:,s]) - xtraj[:,s+1]
            #push!(stance_indices, s)
        end
    else
        for j = 1:(Nm-1)
            s = (Nmodes-1)*Nm + j
            d[:,s] = flight_dynamics_rk4(xtraj[:,s],utraj[:,s]) - xtraj[:,s+1]
            #push!(flight_indices, s)
        end
    end
    
    return nothing
end
     


function stance_constraint!(c,ztraj)
    d = view(c,c_stance_inds)
    z = reshape(ztraj[1:(end-Nx)],Nx+Nu,Nt-1)
    xtraj = [z[1:Nx,:] ztraj[(end-(Nx-1)):end]]
    t = 1
    for k = 1:Nmodes
        if mod(k,2) == 1
            for j = 1:Nm
                s = (k-1)*Nm + j
                d[t] = xtraj[4,s]
                push!(stance_const_indices, s)
                t += 1
            end
        end
            
    end
    return nothing
end

function length_constraint!(c,ztraj)
    d = view(c,c_length_inds)
    z = reshape(ztraj[1:(end-Nx)],Nx+Nu,Nt-1)
    xtraj = [z[1:Nx,:] ztraj[(end-(Nx-1)):end]]
    for k = 1:Nt
        d[k] = norm(xtraj[1:2,k] - xtraj[3:4,k])
    end
end
     
function con!(c,ztraj)
    c[c_init_inds] .= ztraj[1:Nx] - xref[:,1] #initial state constraint
    c[c_term_inds] .= ztraj[(end-(Nx-1)):end] - xref[:,end] #terminal state constraint
    @views dynamics_constraint!(c,ztraj)
    @views stance_constraint!(c,ztraj)
    @views length_constraint!(c,ztraj)
end

struct ProblemMOI <: MOI.AbstractNLPEvaluator
    n_nlp::Int
    m_nlp::Int
    idx_ineq
    obj_grad::Bool
    con_jac::Bool
    sparsity_jac
    sparsity_hess
    primal_bounds
    constraint_bounds
    hessian_lagrangian::Bool
end

function ProblemMOI(n_nlp,m_nlp;
        idx_ineq=c_length_inds,
        obj_grad=true,
        con_jac=true,
        sparsity_jac=sparsity_jacobian(n_nlp,m_nlp),
        sparsity_hess=sparsity_hessian(n_nlp,m_nlp),
        primal_bounds=primal_bounds(n_nlp),
        constraint_bounds=constraint_bounds(m_nlp,idx_ineq=idx_ineq),
        hessian_lagrangian=false)

    ProblemMOI(n_nlp,m_nlp,
        idx_ineq,
        obj_grad,
        con_jac,
        sparsity_jac,
        sparsity_hess,
        primal_bounds,
        constraint_bounds,
        hessian_lagrangian)
end

function primal_bounds(n)
    x_l = -Inf*ones(n)
    x_u = Inf*ones(n)
    return x_l, x_u
end

function constraint_bounds(m; idx_ineq=(1:0))
    c_l = zeros(m)
    c_l[idx_ineq] .= ℓ_min

    c_u = zeros(m)
    c_u[idx_ineq] .= ℓ_max
    return c_l, c_u
end

function row_col!(row,col,r,c)
    for cc in c
        for rr in r
            push!(row,convert(Int,rr))
            push!(col,convert(Int,cc))
        end
    end
    return row, col
end

function sparsity_jacobian(n,m)

    row = []
    col = []

    r = 1:m
    c = 1:n

    row_col!(row,col,r,c)

    return collect(zip(row,col))
end

function sparsity_hessian(n,m)

    row = []
    col = []

    r = 1:m
    c = 1:n

    row_col!(row,col,r,c)

    return collect(zip(row,col))
end

function MOI.eval_objective(prob::MOI.AbstractNLPEvaluator, x)
    cost(x)
end

function MOI.eval_objective_gradient(prob::MOI.AbstractNLPEvaluator, grad_f, x)
    ForwardDiff.gradient!(grad_f,cost,x)
    return nothing
end

function MOI.eval_constraint(prob::MOI.AbstractNLPEvaluator,g,x)
    con!(g,x)
    return nothing
end

function MOI.eval_constraint_jacobian(prob::MOI.AbstractNLPEvaluator, jac, x)
    ForwardDiff.jacobian!(reshape(jac,prob.m_nlp,prob.n_nlp), con!, zeros(prob.m_nlp), x)
    return nothing
end

function MOI.features_available(prob::MOI.AbstractNLPEvaluator)
    return [:Grad, :Jac]
end

MOI.initialize(prob::MOI.AbstractNLPEvaluator, features) = nothing
MOI.jacobian_structure(prob::MOI.AbstractNLPEvaluator) = prob.sparsity_jac

function solve(x0,prob::MOI.AbstractNLPEvaluator;
        tol=1.0e-6,c_tol=1.0e-6,max_iter=10000)
    x_l, x_u = prob.primal_bounds
    c_l, c_u = prob.constraint_bounds

    nlp_bounds = MOI.NLPBoundsPair.(c_l,c_u)
    block_data = MOI.NLPBlockData(nlp_bounds,prob,true)

    solver = Ipopt.Optimizer()
    solver.options["max_iter"] = max_iter
    solver.options["tol"] = tol
    solver.options["constr_viol_tol"] = c_tol

    x = MOI.add_variables(solver,prob.n_nlp)

    for i = 1:prob.n_nlp
        xi = MOI.SingleVariable(x[i])
        MOI.add_constraint(solver, xi, MOI.LessThan(x_u[i]))
        MOI.add_constraint(solver, xi, MOI.GreaterThan(x_l[i]))
        MOI.set(solver, MOI.VariablePrimalStart(), x[i], x0[i])
    end

    # Solve the problem
    MOI.set(solver, MOI.NLPBlock(), block_data)
    MOI.set(solver, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(solver)

    # Get the solution
    res = MOI.get(solver, MOI.VariablePrimal(), x)

    return res
end


#Guess
xguess = xref + 0.1*randn(Nx,Nt)
uguess = kron(ones(Nt-1)', uref) + 0.1*randn(Nu,Nt-1)
z0 = [reshape([xguess[:,1:(Nt-1)]; uguess],(Nx+Nu)*(Nt-1),1); xguess[:,end]];


prob = ProblemMOI(n_nlp,m_nlp)
z_sol = solve(z0,prob) # solve


# z = reshape(z_sol[1:(end-Nx)],Nx+Nu,Nt-1)
# xtraj = [z[1:Nx,:] z_sol[end-(Nx-1):end]]
# utraj = z[(Nx+1):(Nx+Nu),:];

# plot([x[1] for x in eachcol(xtraj)], [y[2] for y in eachcol(xtraj)])
# plot!([x[3] for x in eachcol(xtraj)], [y[4] for y in eachcol(xtraj)])

