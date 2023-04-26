using DifferentialEquations

ε = 0.1
tspan = (0.0, 1.0) # time span (in our case the span of variable x)

@doc raw"""
    boundary_layer_problem!(du, u, p, x; ε=ε)

Specification of the boundary layer problem 
```math
    \varepsilon \cdot y'' = y \cdot y' - y
```
for the solution by `DifferentialEquations`.
"""
function boundary_layer_problem!(du, u, p, x; ε=ε)
    y = u[1]
    dy = u[2]
    du[1] = u[2]
    du[2] = 1 / ε * y * (dy - 1)
end

"Define the boundary layer problem with initial condition `y_prime_0`"
problem_ODE(y_prime_0) = ODEProblem(boundary_layer_problem!, [1.0, y_prime_0], tspan)

# solving the initial value ODE problem with predetermined initial conditions from perturbation method approximations
B1 = solve(problem_ODE(0.99990238), Tsit5(), reltol=1e-8, abstol=1e-10)
M = solve(problem_ODE(0.9051), Tsit5(), reltol=1e-8, abstol=1e-10)
B0 = solve(problem_ODE(-10.6942), Tsit5(), reltol=1e-8, abstol=1e-10)

f = Figure(title=L"Phase plane of the solutions $B0$, $M$ and $B1$")

axes = [Axis(f[1, i], xlabel=L"y", ylabel=L"z") for i in 1:3]

for (i, (ax, lab, sol)) in enumerate(zip(axes, [L"B0", L"M", L"B1"], [B0, M, B1]))
    ax.title[] = lab
    vlines!(ax, [-1, 1]; linestyle=:dash, color=Makie.wong_colors()[4], linewidth=1)    # initial conditions
    lines!(ax, sol[1, :], sol[2, :], color=Makie.wong_colors()[i])                      # phase portrait solutions
    scatter!(ax, sol[1, [1, end]], sol[2, [1, end]], color=Makie.wong_colors()[i])      # start and end points
end

return f
