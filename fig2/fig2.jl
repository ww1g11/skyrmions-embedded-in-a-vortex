using CairoMakie
using DelimitedFiles
using Printf
using MicroMagnetic

function plot_skyrmion()

    data1, units = read_table("../scripts/1skx_vortex/R_1_sd.txt")
    data2, units = read_table("../scripts/1skx_vortex/R_3_sd.txt")
    data3, units = read_table("../scripts/vortex/R_1_sd.txt")
    data4, units = read_table("../scripts/vortex/R_3_sd.txt")


    fig = Figure(; size=(500, 360), fontsize=16)
    ax = Axis(fig[1, 1]; xlabel="Hz (mT)", ylabel="Q")

    scatterlines!(ax, data2["zeeman_Hz"]/mT, data2["Q"]; marker=:star5, label="1-skyrmion vortex, p=+1")
    scatterlines!(ax, data3["zeeman_Hz"]/mT, data3["Q"], marker=:circle, linestyle=:dot, label="DMI vortex, p=+1")
    scatterlines!(ax, data4["zeeman_Hz"]/mT, data4["Q"]; marker=:diamond, linestyle=:dot, label="DMI vortex, p=-1")
    scatterlines!(ax, data1["zeeman_Hz"]/mT, data1["Q"], marker=:rect, label="1-skyrmion vortex, p=-1")
    
    ylims!(ax, [-2, 2])
 
    axislegend()
    save("fig2.pdf", fig)  
    return fig
end

plot_skyrmion()



