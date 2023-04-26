using DifferentialEquations, LaTeXStrings

const g = 9.81                              # gravitational acceleration [m/s²]
L = 1.0                                     # length [m]

@doc raw"""
    simple_pendulum!(du, u, p, t)

Specification of the simple pendulum problem
```math
\varphi'' = -g/L \sin(\varphi)
```

# Parameters
- `du`: The derivative of the state vector.
- `u`: The current state vector.
- `p`: A tuple containing the values of g and L, respectively.
- `t`: The current time.
"""
function simple_pendulum!(du, u, p, t)
    φ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g / L) * sin(φ)
end

u0 = [pi / 4, 0.0]                          # initial conditions (φ, d/dt φ)

BASE_FIGURE_NAME = "comparison_of_numerical_solution_to_small_angles_pendulum_time_span_"

"""
    numerical_v_small_angles(tspan)

Create a figure comparing the numerical solution against the small angles solution.
"""
function numerical_v_small_angles(tspan)
    prob = ODEProblem(simple_pendulum!, u0, tspan)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-10)

    xticks = (tspan[1]:round(Int, tspan[end]), collect(map(i -> latexstring(i, L"\,s"), tspan[1]:round(Int, tspan[end]))))
    f = Figure()
    ax = Axis(f[1, 1],
        xlabel=L"t",
        ylabel=L"\varphi(t)",
        yticks=([pi / 4, pi / 8, 0, -pi / 4, -pi / 8], [L"\pi/4", L"\pi/8", L"0", L"-\pi/4", L"-\pi/8"]),
        xticks=xticks
    )

    lines!(ax, sol.t, sol[1, :], label="Numerical solution")
    lines!(ax, sol.t, u0[1] .* cos.(sqrt(g / L) .* sol.t), label=L"\varphi_0 \cdot \cos{\left(\sqrt{g} \cdot t\right)}")
    axislegend(ax)

    return (BASE_FIGURE_NAME .* join(string.(round.(tspan, sigdigits=2)), "_"), f)
end

return [numerical_v_small_angles(tspan) for tspan in [(0, pi / 2), (0, 2pi), (0, 6pi)]]
