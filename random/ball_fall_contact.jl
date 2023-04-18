using Revise
using DifferentialEquations
using Plots 
gr()
#=
State variable, X = [x, y, ẋ, ẏ]
Dynamics: Ẋ = [ẋ, ẏ, 0, -9.81]
Initial state: X₀ = [0.5, 5.0, 0.2, 0.0]
=#

function dynamics(Ẋ, X, g, t)
    Ẋ[1] = X[3]
    Ẋ[2] = X[4]
    Ẋ[3] = 0.0
    Ẋ[4] = -g 
end

function condition(X, t, integrator)
    X[2]-0.86 
end

function affect!(integrator)
    integrator.u[4] = -0.8*integrator.u[4]
end

cb = ContinuousCallback(condition, affect!)

X₀ = [0.5, 10.0, 2.0, 0.0]
tspan = (0.0, 6.0)
g = 9.806
prob = ODEProblem(dynamics, X₀, tspan, g)
sol = solve(prob, Tsit5(), dt=0.001, callback=cb) # Tsitouras 5/4 Runge-Kutta method

xs = [sol(t)[1] for t=0.0:0.025:tspan[2]]
ys = [sol(t)[2] for t=0.0:0.025:tspan[2]] 
num_frames = length(xs)
anim = Animation()
for i = 1:num_frames 
    plot([0.0, 10.0], [0.0, 0.0], linewidth=4, linecolor=:black, legend=false, 
        grid=false, axis=false)
    scatter!([xs[i]], [ys[i]], xlim=(0, 10), ylim=(-1.0, 12.0),markersize=30) 
    frame(anim) 
end
mp4(anim, "fall_contact.mp4", fps=30)