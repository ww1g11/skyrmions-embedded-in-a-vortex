using MicroMagnetic
using CairoMakie
using Printf

mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=200, ny=200, nz=1);

geo = Cylinder(; radius=200e-9);

function skx1(x,y,z)
    r = sqrt((x)^2 + y^2)
    if r < 20e-9
        return (0, 0, -1)
    end
    return (0, 0, 1)
end

function skx2(x,y,z)
    r = sqrt((x-50e-9)^2 + y^2)
    if r < 20e-9
        return (0, 0, -1)
    end
    r = sqrt((x+50e-9)^2 + y^2)
    if r < 20e-9
        return (0, 0, -1)
    end
    return (0, 0, 1)
end

for (i,f) in enumerate([skx1, skx2])
    args = (
        task = "Relax",
        name = @sprintf("skx_%d", i),
        mesh = mesh,
        shape = geo,
        Ms = 1.1e6,
        A = 1.3e-11,
        Ku = 3e5, 
        D = 2e-3,
        demag = false,
        dmi_type = "interfacial",
        m0 = f,
        stopping_dmdt = 0.01
    );
    sim = sim_with(args);

    vtkname = @sprintf("skx_%d.vts", i)
    save_vtk(sim, vtkname)

    m0 = read_vtk(vtkname)
    m = reshape(m0, 3, 200, 200, 1)
    shape = Cylinder(; radius=100e-9);
    Q = compute_skyrmion_number(m0, mesh, shape)
    println(vtkname, " Q=", Q)

    m0 = read_vtk(vtkname)
    m = reshape(m0, 3, 200, 200, 1)
    fig = plot_m(m[:,75:126, 75:126, :]; figsize=(400, 400), arrows=(20, 20))
    save(@sprintf("vortex_%d.png", i), fig)
end
