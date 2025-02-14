using MicroMagnetic
using CairoMakie
using Printf

mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=200, ny=200, nz=1);

function two_skx1(x,y,z)
    r = sqrt((x-30e-9)^2 + y^2)
    if r < 10e-9
        return (0, 0, -1)
    end
    r = sqrt((x+30e-9)^2 + y^2)
    if r < 10e-9
        return (0, 0, -1)
    end
    return (0, 0, 1)
end

function two_skx2(x,y,z)
    r = sqrt((x-30e-9)^2 + y^2)
    if r < 10e-9
        return (0, 0, 1)
    end
    r = sqrt((x+30e-9)^2 + y^2)
    if r < 10e-9
        return (0, 0, 1)
    end
    return (0, 0, -1)
end

function compute_skyrmion_number_shape(sim)
    shape = Cylinder(; radius=100e-9);
    Q = compute_skyrmion_number(Array(sim.spin), sim.mesh, shape)
    return Q
end

for (i,f) in enumerate([two_skx1, two_skx2]), (j, D) in enumerate([2e-3, -2e-3])
    args = (
        task = "Relax",
        mesh = mesh,
        Ms = 1.1e6,
        A = 1.3e-11,
        Ku = 3.2e5, 
        D = D,
        demag = false,
        m0 = f,
        stopping_dmdt = 0.1
    );
    sim = sim_with(args);

    vtkname = @sprintf("skx_%d.vts", 2*(i-1)+j)
    save_vtk(sim, vtkname)

    Q = compute_skyrmion_number_shape(sim)
    println(vtkname, " Q=", Q)
end


for i in 1:4
    vtkname = @sprintf("skx_%d.vts", i)
    args = (
        task = "Relax",
        mesh = mesh,
        shape = Cylinder(radius=200e-9),
        Ms = 1.1e6,
        A = 1.3e-11,
        Ku = 4e5,
        D = 2e-3,
        dmi_type = "interfacial",
        demag = true,
        m0 = read_vtk(vtkname),
        stopping_dmdt = 0.01
    );
    sim = sim_with(args);

    vtkname = @sprintf("m1_%d.vts", i)
    save_vtk(sim, vtkname)

    Q = compute_skyrmion_number_shape(sim)
    println(vtkname, " Q=", Q)

    println("energy: ", sum(sim.energy))

    m0 = read_vtk(vtkname)
    m = reshape(m0, 3, 200, 200, 1)
    fig = plot_m(m[:, 75:125, 75:125, :]; figsize=(400, 400), arrows=(20, 20))
    save(@sprintf("m1_%d.png", i), fig)
end


"""
for i in 4:4
    vtkname = @sprintf("m1_%d.vts", i)
    args = (
        task = "Relax",
        driver = "LLG",
        alpha = 1.0,
        name = @sprintf("R_%d", i),
        mesh = mesh,
        Ms = 1.1e6,
        A = 1.3e-11,
        Ku = 3.2e5, 
        D = 2e-3,
        dmi_type = "interfacial",
        demag = true,
        m0 = read_vtk(vtkname),
        stopping_dmdt = 0.1,
        H_s = [(0, 0, i*20mT) for i=0:20],
        save_vtk = true,
        saver_item=SaverItem("Q", "<unitless>", compute_skyrmion_number_shape)
    );
    sim = sim_with(args);
end
"""