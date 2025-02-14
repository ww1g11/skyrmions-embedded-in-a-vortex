using MicroMagnetic
using CairoMakie
using Printf

mesh = FDMesh(; dx=2e-9, dy=2e-9, dz=2e-9, nx=300, ny=300, nz=1);

function two_skx(x,y,z)
    r = sqrt((x-50e-9)^2 + y^2)
    if r < 10e-9
        return (0, 0, -1)
    end
    r = sqrt((x+50e-9)^2 + y^2)
    if r < 10e-9
        return (0, 0, -1)
    end
    return (0, 0, 1)
end

function relax_skyrmion()
    args = (
        task = "Relax",
        mesh = mesh,
        Ms = 1.1e6,
        shape = Cylinder(radius=300e-9),
        A = 1.3e-11,
        Ku = 3.2e5, 
        D = 2e-3,
        demag = false,
        m0 = two_skx,
        stopping_dmdt = 0.1
    );
    sim = sim_with(args);
    save_vtk(sim, "skx")
end


function run_simulation(m0, Ku, R)
    vtkname = @sprintf("vtks/Ku_%g_R_%g.vts", Ku, R)
    if isfile(vtkname)
        return true
    end
    args = (task = "Relax", mesh = mesh, Ms = 1.1e6,
        shape = Cylinder(radius=R*1e-9),
        A = 1.3e-11, Ku = Ku*1e5, D = 2e-3, dmi_type = "interfacial",
        demag = true, m0 = m0, stopping_dmdt = 0.05)

    sim = sim_with(args);

    m = reshape(Array(sim.spin), 3, mesh.nx, mesh.ny, mesh.nz)
    fig = plot_m(m[:, 125:175, 125:175, :]; figsize=(400, 400), arrows=(20, 20))

    shape = Cylinder(radius=100e-9)
    Q = compute_skyrmion_number(Array(sim.spin), mesh, shape)
    
    println(vtkname, " Q=", Q)

    pngname = @sprintf("pngs/Ku_%g_R_%g.png", Ku, R)
    save(pngname, fig)

    if abs(Q + 1.5) < 0.1
        save_vtk(sim, vtkname)
        return true
    end
    return false
end

function compute_phase_digarm(Kus, Rs, I, J)
    m, n = length(Kus), length(Rs)
    visited = falses(m, n)
    directions = [(-1,0), (0,-1), (1,0)]  # R should always decrease
    
    queue = [(I, J)]
    visited[I,J] = true
    
    while !isempty(queue)
        i, j = popfirst!(queue)
        Ku, R =  Kus[i], Rs[j]
        vtkname = @sprintf("vtks/Ku_%g_R_%g.vts", Ku, R)
        m0 = read_vtk(vtkname)

        for (dx, dy) in directions
            ni, nj = i+dx, j+dy
            println("checking ", ni, "  ", nj)
            if checkbounds(Bool, Kus, ni) && checkbounds(Bool, Rs, nj) && !visited[ni,nj] && run_simulation(m0, Kus[ni], Rs[nj])
                visited[ni,nj] = true
                push!(queue, (ni,nj))
            end
        end
    end
    
end

mkpath("vtks")
mkpath("pngs")

Kus = [ku for ku in range(2, stop=7, length=26)]
Rs = [r for r in 100:20:300]
I = 11
J = length(Rs)

relax_skyrmion()
m0 = read_vtk("skx.vts")
run_simulation(m0, Kus[I], Rs[J])

compute_phase_digarm(Kus, Rs, I, J)