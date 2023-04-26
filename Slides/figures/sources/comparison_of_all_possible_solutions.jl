using Elliptic, SpecialFunctions, DifferentialEquations, LaTeXStrings

const g = 9.81                              # gravitational acceleration [m/s²]
L = 1.0                                     # length [m]

@doc raw"""
    analytical_pendulum(t; φ_0 = pi / 4)

Compute the angle of a pendulum at time `t` using the analytical solution.

# Parameters:
 - `t`: time at which to compute the angle of the pendulum.
- `φ_0`: initial angle of the pendulum in radians

# Examples
Calculate the analytical solution of the pendulum's angle at time 2.5 seconds, with an initial angle of pi / 6:
```jldoctest
julia> analytical_pendulum(2.5, φ_0=pi / 6, N=6)
0.08271464342761728
```

Calculate the analytical solution of the pendulum's angle at time 1 second, with the default initial angle:
```jldoctest
julia> analytical_pendulum(1)
-0.7789538928851508
```
# Extended help
This function uses the analytical solution for the angle of a pendulum as derived by Duffing in 1918.
The solution is given by:

```math
\varphi(t) =  2 \cdot \arcsin(\sin(\varphi_0 / 2) \cdot \operatorname{cd}(t \cdot \omega,\, k^2))
```
where ``\omega = \sqrt(g / L)``, ``k = sin(\varphi_0 / 2)``, and ``\operatorname{cd}`` is the Jacobi elliptic function.
"""
function analytical_pendulum(t; φ_0=pi / 4)
    if φ_0 < -pi / 2 || φ_0 > pi / 2
        throw(DomainError("Initial angle must be within the range [-pi/2, pi/2]."))
    end

    ω = sqrt(g / L)
    k = sin(φ_0 / 2)
    2 * asin(k * Elliptic.Jacobi.cd(t * ω, k^2))
end

@doc raw"""
    perturbation_approximation_pendulum(t; φ_0=pi / 4, N=4)

Calculate the perturbation approximation of the pendulum's angle at time `t`.

# Parameters
    - `t`: The time at which to calculate the pendulum's angle.
    - `φ_0`: The initial angle of the pendulum in radians.
    - `N`: The number of terms to include in the perturbation approximation.

# Examples
Calculate the perturbation approximation of the pendulum's angle at time 2.5 seconds, with an initial angle of pi / 6 and 6 terms:
```jldoctest
julia> perturbation_approximation_pendulum(2.5, φ_0=pi / 6, N=6)
0.2202397634115413
```
Calculate the perturbation approximation of the pendulum's angle at time 1 second, with the default initial angle and 3 terms:
```jldoctest
julia> perturbation_approximation_pendulum(1)
-0.7797171023953037
```
"""
function perturbation_approximation_pendulum(t; φ_0=pi / 4, N=4)
    if φ_0 < 0 || φ_0 > pi
        throw(DomainError("φ_0 must be within the range of 0 to pi."))
    end

    ω_A = 4 * φ_0 * sqrt(2 * besselj(1, φ_0) / φ_0)
    coef(j, A) = ((-1)^j * besselj(2j + 1, A)) / (4j * (j + 1) * besselj(1, A))
    coefs = coef.(2 .* (1:N) .+ 1, φ_0)
    φ_0 * (1 - sum(coefs)) * cos(ω_A * t) + sum(coefs .* cos.((2 .* (1:N) .+ 1) .* ω_A .* t))
end

@doc raw"""
    simple_pendulum!(du, u, p, t)

Specification of the simple pendulum problem
```math
\varphi'' = -g/L \sin(\varphi)
```

# Parameters
- `du`: The derivative of the state vector.
- `u`: The current state vector.
- `t`: The current time.
"""
function simple_pendulum!(du, u, p, t)
    φ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g / L) * sin(φ)
end

u0 = [pi / 4, 0.0]                          # initial conditions (φ, d/dt φ)

BASE_FIGURE_NAME = "comparison_of_all_the_solutions_time_span_"

"""
    all_solutions(tspan)

Create a figure comparing all the solutions of the simple pendulum problem.
"""
function all_solutions(tspan)
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
    lines!(ax, sol.t, perturbation_approximation_pendulum.(sol.t), label="Perturbation Methods approximation")
    lines!(ax, sol.t, analytical_pendulum.(sol.t), label=L"2 \cdot \arcsin(\sin(\varphi_0 / 2) \cdot \text{cd} (t \cdot \sqrt{g / L},\, \sin^2(\varphi_0 / 2)))")
    axislegend(ax)

    return (BASE_FIGURE_NAME .* join(string.(round.(tspan, sigdigits=2)), "_"), f)
end

return [all_solutions(tspan) for tspan in [(0, pi / 2), (0, 2pi), (0, 6pi)]]
