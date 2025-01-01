using CairoMakie
using ChebParticleMesh

L = 10.0, 10.0
periodicity = (true, false)
extra_pad = (0, 10) # no extra padding needed
N_real = (100, 100) # 16 grids in each direction in real space
w = (20, 20)
α = 0.5

gridinfo = GridInfo(N_real, w, periodicity, extra_pad, L)
gridbox = GridBox(gridinfo)

# using the Wkb window function, and generate its Chebyshev approximation
f_window = [x -> Wkb(x, (w[i] + 0.5) * gridinfo.h[i], 5.0 * w[i]) for i in 1:2]
cheb_coefs = tuple([ChebCoef(f_window[i], gridinfo.h[i], w[i], 10) for i in 1:2]...)

T = Float64

# using the Fourier transform of the window function and Green's function in Fourier space to generate the scaling factors
F_f_window = [x -> FWkb(x, (w[i] + 0.5) * gridinfo.h[i], 5.0 * w[i]) for i in 1:2]
func_scale = (kx, ky) -> iszero(kx^2 + ky^2) ? zero(T) : (F_f_window[1](kx) * F_f_window[2](ky))^(-2) * exp(-(kx^2 + ky^2) / (4*α^2)) / (kx^2 + ky^2)
scalefactor = ScalingFactor(func_scale, gridinfo)

n_atoms = 10
qs = [(-1.0)^i for i in 1:n_atoms]
poses = [L .* (rand(), rand()) for i in 1:n_atoms]


interpolate!(qs, poses, gridinfo, gridbox, cheb_coefs)

data_1 = copy(real.(gridbox.pad_grid))

fft!(gridbox)

data_2 = copy(real.(gridbox.pad_grid))

ChebParticleMesh.scale!(gridbox, scalefactor)

data_3 = copy(real.(gridbox.pad_grid))

ifft!(gridbox)

data_4 = copy(real.(gridbox.pad_grid))

E = gather(qs, poses, gridinfo, gridbox, cheb_coefs) / 8π

fig_1 = Figure(size = (300, 300))
ax_1 = Axis(fig_1[1, 1])
hidedecorations!(ax_1)
scatter!(ax_1, poses, color = [i % 2 == 0 ? :red : :blue for i in 1:n_atoms])
fig_1
save("figs/particlemesh_1.svg", fig_1)

s = 3 .* size(data_1)
style = :balance

fig_2 = Figure(size = s)
ax_2 = Axis(fig_2[1, 1])
hidedecorations!(ax_2)
heatmap!(ax_2, data_1, colormap = style)
fig_2
save("figs/particlemesh_2.svg", fig_2)

fig_3 = Figure(size = s)
ax_3 = Axis(fig_3[1, 1])
hidedecorations!(ax_3)
heatmap!(ax_3, data_2, colormap = style)
fig_3
save("figs/particlemesh_3.svg", fig_3)

fig_4 = Figure(size = s)
ax_4 = Axis(fig_4[1, 1])
hidedecorations!(ax_4)
heatmap!(ax_4, data_3, colormap = style)
fig_4
save("figs/particlemesh_4.svg", fig_4)

fig_5 = Figure(size = s)
ax_5 = Axis(fig_5[1, 1])
hidedecorations!(ax_5)
heatmap!(ax_5, data_4, colormap = style)
fig_5
save("figs/particlemesh_5.svg", fig_5)