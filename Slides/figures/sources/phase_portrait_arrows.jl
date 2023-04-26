arrow_number = 25 # number of arrows in a row and column
ys = LinRange(-2, 2, arrow_number) # ys to generate arrows on
ε = 0.1 # epsilon of the DE

lengthscale = 0.01 # scale the arrows by this number

# Generates an arrow for each point at the predefined positions of the derivatives of the phase field
f = Figure()
ax = Axis(f[1, 1], xlabel=L"y", ylabel=L"z")
for z in LinRange(-15, 5, arrow_number)
    zs = fill(z, arrow_number)
    zs′ = [1 / ε * y * (z - 1) for y in ys]
    ys′ = [z for z in zs]
    arrows!(ax,
        ys .- (lengthscale .* vec(ys′)) ./ 2,
        zs .- (lengthscale .* vec(zs′)) ./ 2,
        ys′,
        zs′,
        lengthscale=lengthscale, arrowsize=5)
end

xlims!(ax, -2, 2)
ylims!(ax, -15, 5)

return f
