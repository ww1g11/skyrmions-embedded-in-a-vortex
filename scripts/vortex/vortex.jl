using MicroMagnetic
using CairoMakie
using Printf

mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=200, ny=200, nz=1);

geo = Cylinder(; radius=200e-9);

function up_cw(x, y, z)
    r = sqrt(x^2 + y^2)
    if r < 20e-9
        return (0, 0, 1)
    end
    return (y / r, -x / r, 0)
end

function up_acw(x, y, z)
    r = sqrt(x^2 + y^2)
    if r < 20e-9
        return (0, 0, 1)
    end
    return (-y / r, x / r, 0)
end

function down_cw(x, y, z)
    r = sqrt(x^2 + y^2)
    if r < 20e-9
        return (0, 0, -1)
    end
    return (y / r, -x / r, 0)
end

function down_acw(x, y, z)
    r = sqrt(x^2 + y^2)
    if r < 20e-9
        return (0, 0, -1)
    end
    return (-y / r, x / r, 0)
end

function compute_skyrmion_number_shape(sim)
    shape = Cylinder(; radius=100e-9);
    Q = compute_skyrmion_number(Array(sim.spin), sim.mesh, shape)
    return Q
end

for (i,f) in enumerate([up_cw, up_acw, down_cw, down_acw])
    args = (
        task = "Relax",
        mesh = mesh,
        Ms = 1.1e6,
        A = 1.3e-11,
        #Ku = 4e5, 
        shape = geo,
        demag = true,
        m0 = f,
        stopping_dmdt = 0.01
    );
    sim = sim_with(args);

    vtkname = @sprintf("vortex_%d.vts", i)
    save_vtk(sim, vtkname)

    Q = compute_skyrmion_number_shape(sim)
    println(vtkname, " Q=", Q)
    println("energy: ", sum(sim.energy))

    m0 = read_vtk(vtkname)
    m = reshape(m0, 3, 200, 200, 1)
    fig = plot_m(m[:,75:126, 75:126, :]; figsize=(400, 400), arrows=(20, 20))
    save(@sprintf("vortex_%d.png", i), fig)
end


for (i,f) in enumerate([up_cw, up_acw, down_cw, down_acw])
    vtkname = @sprintf("vortex_%d.vts", i)
    args = (
        task = "Relax",
        mesh = mesh,
        shape = geo,
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

    vtkname = @sprintf("dmi_vortex_%d.vts", i)
    save_vtk(sim, vtkname)

    Q = compute_skyrmion_number_shape(sim)
    println(vtkname, " Q=", Q)
    println("energy: ", sum(sim.energy))

    m = reshape(Array(sim.spin), 3, 200, 200, 1)
    fig = plot_m(m[:, 75:126, 75:126, :]; figsize=(400, 400), arrows=(20, 20))
    save(@sprintf("dmi_vortex_%d.png", i), fig)
end

"""
for i = 1:4
    vtkname = @sprintf("vortex_%d.vts", i)
    args = (
        task = "Relax",
        name = @sprintf("R_%d", i),
        mesh = mesh,
        shape = geo,
        Ms = 1.1e6,
        A = 1.3e-11,
        Ku = 4e5, 
        D = 2e-3,
        dmi_type = "interfacial",
        demag = true,
        m0 = read_vtk(vtkname),
        stopping_dmdt = 0.01,
        H_s = [(0, 0, i*20mT) for i=0:20],
        save_vtk = true,
        saver_item=SaverItem("Q", "<unitless>", compute_skyrmion_number_shape)
    );
    sim = sim_with(args);
end
"""