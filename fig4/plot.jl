using DelimitedFiles
using CairoMakie

function plot()

    vortex = readdlm("vortex_energy.txt", header=true)[1]
    skx_vortex = readdlm("skx_vortex_energy.txt", header=true)[1]

    eV = 1.602176565e-19
    
    fig = Figure(; size=(500, 360), fontsize=18)
    ax = Axis(fig[1, 1]; xlabel="Q", ylabel="Energy (eV)")

    s2 = scatterlines!(ax, vortex[:, 1], vortex[:, 2]/eV; marker=:circle, markersize=10, color=:sienna1, strokecolor=:sienna1, label="Vortex")
    s1 = scatterlines!(ax, skx_vortex[:, 1], skx_vortex[:, 2]/eV; marker=:rect, markersize=10, color=:slateblue1, strokecolor=:slateblue1, label="n-skyrmion Vortex")
    
    axislegend(ax; position=(0.5, 0.8), labelsize=14)

    save("fig4.pdf", fig)
    return nothing
end

plot()