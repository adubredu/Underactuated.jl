using Revise 
using StaticArrays
using JuMP 
using OSQP 
using Underactuated
using LinearAlgebra
import Ipopt 
using Plots 

# x = [pb_x, pb_y, pf_x, pf_y, vb_x, vb_y, vf_x, vf_y]
# u = [F, τ]
Nx = 8
Nu = 2 
Tfinal = 4.4
h = 0.1
Nm = 5
Nt = Int(ceil(Tfinal/h)+1)
Nmodes = Int(ceil(Nt/Nm))
thist = Array(range(0,h*(Nt-1), step=h))

# Hopper Dynamics 
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

function flight_dynamics_rk4(z...)
    x=collect(z[1:Nx]) 
    u=collect(z[(Nx+1):(Nx+Nu)])
    #RK4 integration with zero-order hold on u
    f1 = flight_dynamics(x, u)
    f2 = flight_dynamics(x + 0.5*h*f1, u)
    f3 = flight_dynamics(x + 0.5*h*f2, u)
    f4 = flight_dynamics(x + h*f3, u)
    return x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
end

function stance_dynamics_rk4(z...)
    x=collect(z[1:Nx]) 
    u=collect(z[(Nx+1):(Nx+Nu)])
    #RK4 integration with zero-order hold on u
    f1 = stance_dynamics(x, u)
    f2 = stance_dynamics(x + 0.5*h*f1, u)
    f3 = stance_dynamics(x + 0.5*h*f2, u)
    f4 = stance_dynamics(x + h*f3, u)
    return x + (h/6.0)*(f1 + 2*f2 + 2*f3 + f4)
end

function jump_map(z...)
    #Assume the foot experiences inelastic collisions
    x = flight_dynamics_rk4(z...)
    xn = [x[1:6]; 0.0; 0.0]
    return xn
end

function dist(x...)
    z = collect(x)
    return norm(z)
end

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


Q = Diagonal([1.0*ones(4); 1.0*ones(4)]);
R = 0.001


#Reference Trajectory
uref = [m1*g; 0.0]
xref = zeros(Nx,Nt)
xref[1,:] .= LinRange(-1.0,1.0,Nt)
xref[2,:] .= 1.0 .+ 0.5.*sin.(2*pi/10.0*(0:(Nt-1)));
xref[3,:] .= LinRange(-1.0,1.0,Nt)
xref[5,2:end-1] .= (2.0/Tfinal)*ones(Nt-2)
xref[7,2:end-1] .= (2.0/Tfinal)*ones(Nt-2);


model = Model(Ipopt.Optimizer)
@variable(model, x[1:Nx, 1:Nt])
@variable(model, u[1:Nu, 1:Nt])

xguess = xref + 0.1*randn(Nx,Nt)
uguess = kron(ones(Nt)', uref) + 0.1*randn(Nu,Nt)
set_start_value.(x, xguess)
set_start_value.(u, uguess)

@objective(model, Min, sum(0.5*((x[:,i]-xref[:,i])'*Q*(x[:,i]-xref[:,i])) + 
                            0.5*(u[:,i]-uref)'*R*(u[:,i]-uref) for i=1:Nt))


flight_memoized_dynamics = memoize(flight_dynamics_rk4, Nx)
stance_memoized_dynamics = memoize(stance_dynamics_rk4, Nx) 
jumpmap_memoized_dynamics = memoize(jump_map, Nx)

for i=1:Nx
    register(model, Symbol(:f, i), Nx+Nu, flight_memoized_dynamics[i]; autodiff=true) 
    register(model, Symbol(:s, i), Nx+Nu, stance_memoized_dynamics[i]; autodiff=true)
    register(model, Symbol(:j, i), Nx+Nu, jumpmap_memoized_dynamics[i]; autodiff=true)
end 
register(model, :dist, 2, dist; autodiff=true)

counter=0
mode = :stance

@constraint(model, [i=1:Nx], x[i, 1] == xref[i, 1])
@constraint(model, [i=1:Nx], x[i, Nt] == xref[i, Nt])
stance_indices = []; flight_indices = []; jump_indices = [];
for i=1:Nt-1 
    global counter += 1 
    global mode
    # @show counter, mode
    if mode == :stance   
        @NLconstraint(model, x[1, i+1] == s1(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
        @NLconstraint(model, x[2, i+1] == s2(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
        @NLconstraint(model, x[3, i+1] == s3(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
        @NLconstraint(model, x[4, i+1] == s4(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
        @NLconstraint(model, x[5, i+1] == s5(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
        @NLconstraint(model, x[6, i+1] == s6(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
        @NLconstraint(model, x[7, i+1] == s7(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
        @NLconstraint(model, x[8, i+1] == s8(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i])) 
        @constraint(model, x[4, i] == 0.0)
        push!(stance_indices, i)
    elseif mode == :flight 
        if counter == Nm
            @NLconstraint(model, x[1, i+1] == j1(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[2, i+1] == j2(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[3, i+1] == j3(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[4, i+1] == j4(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[5, i+1] == j5(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[6, i+1] == j6(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[7, i+1] == j7(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[8, i+1] == j8(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i])) 
            push!(jump_indices, i)
        else
            @NLconstraint(model, x[1, i+1] == f1(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[2, i+1] == f2(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[3, i+1] == f3(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[4, i+1] == f4(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[5, i+1] == f5(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[6, i+1] == f6(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[7, i+1] == f7(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i]))
            @NLconstraint(model, x[8, i+1] == f8(x[1, i], x[2, i], x[3, i], x[4, i], x[5, i], x[6, i], x[7, i], x[8, i], u[1, i], u[2, i])) 
            push!(flight_indices, i)
        end
    end

    if counter == Nm
        mode = mode == :stance ? :flight : :stance
        counter = 0
    end 
end


# for k=1:(Nmodes-1)
#     if mod(k,2) == 1
#         for j=1:Nm 
#             s = (k-1)*Nm + j  
#             @NLconstraint(model, x[1, s+1] == s1(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[2, s+1] == s2(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[3, s+1] == s3(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[4, s+1] == 0.0)#s4(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[5, s+1] == s5(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[6, s+1] == s6(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[7, s+1] == s7(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[8, s+1] == s8(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             # @NLconstraint(model, x[4, s] == 0.0)
#         end
#     else
#         for j=(Nm-1)
#             s = (k-1)*Nm + j 
#             @NLconstraint(model, x[1, s+1] == f1(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[2, s+1] == f2(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[3, s+1] == f3(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[4, s+1] == f4(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[5, s+1] == f5(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[6, s+1] == f6(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[7, s+1] == f7(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#             @NLconstraint(model, x[8, s+1] == f8(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         end
#         s = k*Nm 
#         @NLconstraint(model, x[1, s+1] == j1(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[2, s+1] == j2(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[3, s+1] == j3(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[4, s+1] == j4(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[5, s+1] == j5(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[6, s+1] == j6(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[7, s+1] == j7(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[8, s+1] == j8(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#     end
# end
# if mod(Nmodes, 2) == 1
#     for j=1:(Nm-1)
#         s = (Nmodes-1)*Nm + j
#         @NLconstraint(model, x[1, s+1] == s1(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[2, s+1] == s2(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[3, s+1] == s3(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[4, s+1] == 0.0) #s4(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[5, s+1] == s5(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[6, s+1] == s6(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[7, s+1] == s7(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[8, s+1] == s8(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#     end 
# else
#     for j=1:(Nm-1)
#         s = (Nmodes-1)*Nm + j 
#         @NLconstraint(model, x[1, s+1] == f1(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[2, s+1] == f2(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[3, s+1] == f3(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[4, s+1] == f4(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[5, s+1] == f5(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[6, s+1] == f6(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[7, s+1] == f7(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#         @NLconstraint(model, x[8, s+1] == f8(x[1, s], x[2, s], x[3, s], x[4, s], x[5, s], x[6, s], x[7, s], x[8, s], u[1, s], u[2, s]))
#     end
# end 

for i=1:Nt 
    @NLconstraint(model, ℓ_min ≤ dist(x[1, i] - x[3, i], x[2, i]-x[4, i]) ≤ ℓ_max)
end





optimize!(model)

xtraj = value.(x)
utraj = value.(u)
plot([x[1] for x in eachcol(xtraj)], [y[2] for y in eachcol(xtraj)])
plot!([x[3] for x in eachcol(xtraj)], [y[4] for y in eachcol(xtraj)])