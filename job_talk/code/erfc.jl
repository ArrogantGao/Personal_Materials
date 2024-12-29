using SpecialFunctions
using CairoMakie, LaTeXStrings

function plot_erfc()

    fig = Figure(size = (400, 300))
    ax = Axis(fig[1, 1])

    x = range(0, 10, length=5000)

    lines!(ax, x, erfc.(x) ./ x, label="erfc(x)/x")
    lines!(ax, x, erf.(x) ./ x, label="erf(x)/x")
    axislegend(ax)

    xlims!(0, 5)
    ylims!(-0.1, 2)

    return fig
end

fig = plot_erfc()
save("figs/erfc.png", fig)
