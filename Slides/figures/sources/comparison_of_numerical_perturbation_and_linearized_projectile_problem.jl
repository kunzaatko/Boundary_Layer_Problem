using DifferentialEquations

tspan = (0.0, 2.0) # time span for the projectile problem solution (how long to evaluate the trajectory of the solution)
ε = 0.1

"""
    linearized_projectile(τ)

Evaluate the projectile at time `τ` using a solution of the DE, where the small term is not considered at all. (linearized DE)
"""
linearized_projectile(τ) = τ - (τ^2) / 2

function simple_projectile!(du, u, p, t)
    y = u[1]
    dy = u[2]
    du[1] = dy
    du[2] = -1 / ((1 + ε * y)^2)
end

projectile_problem = ODEProblem(simple_projectile!, [0.0, 1.0], tspan)
sol = solve(projectile_problem, Tsit5(), reltol=1e-10, abstol=1e-15)

asymptotic_projectile(τ) = τ * (1 - τ / 2) + ε * (τ^3) * (1 - τ / 4) / 3

f = Figure()
ax = Axis(f[1, 1], xlabel=L"\tau", ylabel=L"y(\tau)")

lines!(ax, sol.t, sol[1, :], label="Numerical solution")
lines!(ax, sol.t, t -> linearized_projectile(t), label="Linearized solution")
lines!(ax, sol.t, t -> asymptotic_projectile(t), label="Two-term asymptotic approximation")

axislegend(ax)

return f
