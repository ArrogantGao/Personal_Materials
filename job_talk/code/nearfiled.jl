using CairoMakie, LaTeXStrings
using FastSpecSoG
using Enzyme

preset = 4
uspara = USeriesPara(preset)

f = x -> abs(1 / x - U_series(x, uspara))

df = x -> abs(autodiff(Forward, f, Duplicated, Duplicated(x, 1.0))[1])

x = range(0.1, 10, 100)

fig = Figure(size=(800, 350), title = L"r_c = 10.0")
ax_1 = Axis(fig[1, 1], xlabel=L"r", ylabel=L"|U_{\text{near}}(r)|", yscale=log10)
ax_2 = Axis(fig[1, 2], xlabel=L"r", ylabel=L"|F_{\text{near}}(r)|", yscale=log10)

lines!(ax_1, x, f.(x), label="f")
lines!(ax_2, x, df.(x), label="df")

fig
save("figs/nearfield.svg", fig, px_per_unit=2)