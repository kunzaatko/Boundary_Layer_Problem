using LinearAlgebra

Δ = 1e-4 # step size
ε = 0.1 # epsilon of the DE

ylims = (-2, 2) # (ymin, ymax) to draw
zlims = (-15, 5) # (zmin, zmax) to draw

"Check if `y` and `z` are out predefined bounds."
function out_of_bounds(y, z; ylims=ylims, zlims=zlims)
    return y < ylims[1] || y > ylims[2] || z < zlims[1] || z > zlims[2]
end

"Velocity of the imaginary particle at points `y` and `z`"
velocity(y, z) = [z, 1 / ε * y * (z - 1)]

@doc raw"""
    trace_phase_portrait(y0, z0; Δ=Δ, min_traces=1000)

Generate a phase portrait by tracing the trajectory of the system of differential equations. Return a list of points that trace the trajectory of an imaginary particle in the phase portrait.

# Parameters
    - `y0`: The initial value of the first variable in the system.
    - `z0`: The initial value of the second variable in the system.
    - `Δ`: The step size used to trace the trajectory. Default value is 0.01.
    - `min_traces`: The minimum number of traces to generate. Default value is 1000. Note that this has to be big enough to "escape" the distance that makes the logic assume that the begining and the end of the path joined
"""
function trace_phase_portrait(y0, z0; Δ=Δ, min_traces=1000)
    trace = [[y0, z0]]
    for _ in 1:min_traces
        push!(trace, trace[end] .+ Δ .* velocity(trace[end]...))
    end
    while !(out_of_bounds(trace[end]...) || norm(trace[1] - trace[end]) < 100 * Δ)
        push!(trace, trace[end] .+ Δ .* velocity(trace[end]...))
    end
    if !(out_of_bounds(trace[end]...))
        push!(trace, trace[1])
    else
        while !(out_of_bounds(trace[1]...))
            pushfirst!(trace, trace[1] .- Δ .* velocity(trace[1]...))
        end
    end
    trace
end

# initial values to trace the track for
initial_points = [
    [0.2, 0.0],
    [0.5, 0.0],
    [0.75, 0.0],
    [1.0, 0.0],
    [1.3, 0.0],
    [1.5, 0.0],
    [0.0, 1.0],
    [0.0, 1.5],
    [0.0, 2.0],
    [0.0, 3.0],
    [0.0, 4.0],
    [-1.5, 2.5],
    [-1.0, 2.5],
]

"""
    figure_of_traces(initial_points)

Generate a figure of phase portrait trajectories of the system of differential equations.
"""
function figure_of_traces(initial_points)
    f = Figure()
    ax = Axis(f[1, 1], aspect=1.5, xlabel=L"y", ylabel=L"z")
    for point in initial_points
        trace = hcat(trace_phase_portrait(point...)...)

        lines!(ax, trace[1, :], trace[2, :],
            color=all(trace[:, 1] .== trace[:, end]) ? Makie.wong_colors()[1] : Makie.wong_colors()[2])

        # FIX: Why is this happening <24-04-23> 
        if !(point in [
            [-1.5, 2.5],
            [-1.0, 2.5],
        ])
            scatter!(ax, point[1], point[2], rotations=[atan(velocity(point...)[2], velocity(point...)[1])],
                color=all(trace[:, 1] .== trace[:, end]) ? Makie.wong_colors()[1] : Makie.wong_colors()[2], marker=:rtriangle, markersize=15)
        end
    end

    xlims!(ax, -2, 2)
    ylims!(ax, -15, 5)

    return f, ax
end

f_no_init, _ = figure_of_traces(initial_points)

f_init, ax = figure_of_traces(initial_points)
vlines!([-1, 1], linestyle=:dash, color=Makie.wong_colors()[4], linewidth=2, label="Boundary condition hyperplanes")
axislegend(ax)

return [("phase_portrait_imaginary_particle_trace_no_initial_vlines", f_no_init),
    ("phase_portrait_imaginary_particle_trace_initial_vlines", f_init)]



