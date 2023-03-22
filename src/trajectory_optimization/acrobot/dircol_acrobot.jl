function memoize(foo::Function, n_outputs::Int)
    last_x, last_f = nothing, nothing
    last_dx, last_dfdx = nothing, nothing
    function foo_i(i, x::T...) where {T<:Real}
        if T == Float64
            if x != last_x
                last_x, last_f = x, foo(x...)
            end
            return last_f[i]::T
        else
            if x != last_dx
                last_dx, last_dfdx = x, foo(x...)
            end
            return last_dfdx[i]::T
        end
    end
    return [(x...) -> foo_i(i, x...) for i in 1:n_outputs]
end

function compute_trajectory(x0, xgoal, robot, T; h=0.1, qcost=1.0, rcost=0.1)
    Nx = state_dim(robot)
    Nu = control_dim(robot)
    Nsteps = Int(T/h) + 1
    Thist = collect(0.0 : h : h*(Nsteps-1))
    Q = Diagonal(qcost*ones(Nx))
    R = rcost

    #Hermite-Simpson integration with first-order hold on u
    function dircol_dynamics(z...)
        x1=collect(z[1:Nx]) 
        u1=collect(z[(Nx+1):(Nx+Nu)]) 
        x2=collect(z[(Nx+Nu)+1 : (Nx+Nu)+Nx]) 
        u2 = collect(z[(Nx+Nu)+Nx+1 : 2(Nx+Nu)]) 
    
        f1 = RZ.dynamics(robot, x1, u1)
        f2 = RZ.dynamics(robot, x2, u2)
        xm = 0.5*(x1 + x2) + (h/8.0)*(f1 - f2)
        um = 0.5*(u1 + u2)
        ẋm = (-3/(2.0*h))*(x1 - x2) - 0.25*(f1 + f2)
        fm = RZ.dynamics(robot, xm, um)
        Δf = fm - ẋm
        return Δf
    end

    memoized_dynamics = memoize(dircol_dynamics, Nx)
    model = Model(Ipopt.Optimizer)
    register(model, :f1, 2Nx+2Nu, memoized_dynamics[1]; autodiff=true)
    register(model, :f2, 2Nx+2Nu, memoized_dynamics[2]; autodiff=true)
    register(model, :f3, 2Nx+2Nu, memoized_dynamics[3]; autodiff=true)
    register(model, :f4, 2Nx+2Nu, memoized_dynamics[4]; autodiff=true)
    @variable(model, x[1:Nx, 1:Nsteps])
    @variable(model, u[1:Nu, 1:Nsteps])
    @objective(model, Min, sum(0.5*((x[:,i]-xgoal)'*Q*(x[:,i]-xgoal)) + 
                            0.5*u[:,i]'*R*u[:,i] for i=1:Nsteps))
    @constraint(model, x[:, 1] .== x0)
    @constraint(model, x[:, Nsteps] .== xgoal)
    for i=1:Nsteps-1
        @NLconstraint(model, f1(x[1,i], x[2,i], x[3,i], x[4,i],
                                u[1,i],
                                x[1,i+1], x[2,i+1], x[3,i+1], x[4,i+1],
                                u[1, i+1]) 
                                == 0.0)
        @NLconstraint(model, f2(x[1,i], x[2,i], x[3,i], x[4,i],
                                u[1,i],
                                x[1,i+1], x[2,i+1], x[3,i+1], x[4,i+1],
                                u[1, i+1]) 
                                == 0.0)
        @NLconstraint(model, f3(x[1,i], x[2,i], x[3,i], x[4,i],
                                u[1,i],
                                x[1,i+1], x[2,i+1], x[3,i+1], x[4,i+1],
                                u[1, i+1]) 
                                == 0.0)
        @NLconstraint(model, f4(x[1,i], x[2,i], x[3,i], x[4,i],
                                u[1,i],
                                x[1,i+1], x[2,i+1], x[3,i+1], x[4,i+1],
                                u[1, i+1]) 
                                == 0.0)
    end
    optimize!(model)
    xtraj = value.(x)
    utraj = value.(u)
    return xtraj, utraj, Thist, Nsteps
end