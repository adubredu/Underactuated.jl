using Revise
using Plots
using DifferentialEquations

gr()  # use the GR backend

# Constants and initial conditions
g = 9.81
m1, m2 = 1.0, 1.0
L1, L2 = 1.0, 1.0
θ1₀, θ2₀ = π/4, π/4
ω1₀, ω2₀ = 0.0, 0.0

# Define the equations of motion for the double pendulum
function double_pendulum!(dY, Y, _, t)
    θ1, θ2, ω1, ω2 = Y

    Δθ = θ2 - θ1
    c = cos(Δθ)
    s = sin(Δθ)

    dY[1] = ω1
    dY[2] = ω2
    dY[3] = (m2 * g * sin(θ2) * c - m2 * s * (L1 * ω1^2 * c + L2 * ω2^2) - (m1 + m2) * g * sin(θ1)) / L1 / (m1 + m2 * s^2)
    dY[4] = ((m1 + m2) * (L1 * ω1^2 * s - g * sin(θ2) + g * sin(θ1) * c) + m2 * L2 * ω2^2 * s * c) / L2 / (m1 + m2 * s^2)
end

# Define the initial state and time span
Y₀ = [θ1₀, θ2₀, ω1₀, ω2₀]
tspan = (0.0, 10.0)

# Solve the equations of motion using a numerical integration method
prob = ODEProblem(double_pendulum!, Y₀, tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# Convert the solution to Cartesian coordinates
x1 = L1 .* sin.(sol[1, :])
y1 = -L1 .* cos.(sol[1, :])
x2 = x1 + L2 .* sin.(sol[2, :])
y2 = y1 - L2 .* cos.(sol[2, :])

# Initialize the animation
animation = Animation()

# Loop over time steps, updating the plot for each frame
for i in 1:length(sol.t)
    p = plot([0, x1[i], x2[i]], [0, y1[i], y2[i]], xlim=(-(L1+L2), L1+L2), ylim=(-(L1+L2), L1+L2), lw=2, legend=false, markershape=:circle, markersize=10, aspect_ratio=:equal)
    frame(animation, p)
end

# Save the animation as an MP4 video file
mp4(animation, "double_pendulum.mp4", fps=25)
