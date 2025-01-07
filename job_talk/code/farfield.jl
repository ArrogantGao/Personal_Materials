using CairoMakie, LaTeXStrings
using FastSpecSoG

preset = 4
uspara = USeriesPara(preset)
M_mid = 5

f_mid(x) = sum(w * exp(- x^2 / s^2) for (s, w) in uspara.sw[1:M_mid])
f_long(x) = sum(w * exp(- x^2 / s^2) for (s, w) in uspara.sw[M_mid+1:end])

x = range(-0, 15, 100)

fig = Figure(size=(400, 300))
ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"f(x)")

lines!(ax, x, f_mid.(x), label="mid")
lines!(ax, x, f_long.(x), label="long")

# vlines!(ax, -5, color=:red, linestyle=:dash)
vlines!(ax, 10, color=:red, linestyle=:dash)

xlims!(ax, 0, 15)
ylims!(ax, 0, 0.32)

axislegend(ax, position=:rt)

fig
save("figs/farfield.svg", fig)