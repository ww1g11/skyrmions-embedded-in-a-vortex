using CairoMakie
using DelimitedFiles
using MicroMagnetic
using Printf

mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=300, ny=300, nz=1);

function compute_Qs(Kus, Rs)
    m, n = length(Kus), length(Rs)
    Qs = zeros(m, n)
    shape = Cylinder(radius=100e-9)
    for i = 1:m, j=1:n
        vtkname = @sprintf("../scripts/phase_1skx_vortex/vtks/Ku_%g_R_%g.vts", Kus[i], Rs[j])
        !isfile(vtkname) && continue
        m0 = read_vtk(vtkname)
        Q = compute_skyrmion_number(m0, mesh, shape)
        println(i," ", j, " ", Q)
        Qs[i, j] = abs(Q + 1.5) < 0.1 ? 1 : 0
    end
    return Qs
end

function plot(Kus, Rs, Qs)

    fig = Figure(size=(500, 360), fontsize = 16)
    ax1 = Axis(fig[1, 1], xlabel=L"Ku (MJ/m$^3$)", ylabel="R (nm)")
    contourf!(ax1, Kus, Rs, Qs, levels=[-0.5, 0.5, 1.5], colormap=[(:blue,0.5), (:red, 0.5)])
    #heatmap!(ax1, Kus, Rs, Qs, colormap=[(:blue,0.5), (:red, 0.5)])
    tightlimits!(ax1)
    text!(3.4, 195, text = "1-skyrmion vortex")

    save("fig3.pdf", fig)
    fig
end

Kus = [ku for ku in range(2, stop=7, length=26)]
Rs = [r for r in 100:20:300]
Qs = compute_Qs(Kus, Rs)

plot(Kus, Rs, Qs)