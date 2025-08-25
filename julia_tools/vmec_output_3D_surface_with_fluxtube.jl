using VMEC, PlasmaEquilibriumToolkit, GLMakie, DelimitedFiles

woutfile = "/u/nigu/gene-stuff/geometry/vmec/wout_kjm_beta-0.nc"

s=0.5

#α = 0
#α = 2*pi/5
α = pi/5

pol_turns = 1 
nz0 = 1000*pol_turns

step_θ = 100
step_ζ = 100

vmec = readVmecWout(woutfile)
vmec_surf = VmecSurface(s, vmec)

θ_r = range(0,2π, step_θ)
ζ_r = range(0,2π, step_ζ)

v = VmecCoordinates(s, 0, 0)

ψ = PestFromVmec()(v, vmec_surf).ψ

grid = MagneticCoordinateGrid(VmecCoordinates, s, θ_r, ζ_r)

cartGrid = CartesianFromVmec()(grid, vmec_surf)

Xs = [q[1] for q in cartGrid];
Ys = [q[2] for q in cartGrid];
Zs = [q[3] for q in cartGrid];

surfmat = [Xs Ys Zs]

writedlm("/u/nigu/gene-stuff/geometry/vmec/3D_vmec_surface_vis/surf_kjm_test_mat.dat", surfmat)

# plotting flux surface

p = surface(Xs, Ys, Zs, color = :blue)

# adding flux tube to surface

ι = -vmec_surf.iota[1]
ζ_range = vcat(range(α-pol_turns*π/ι,step=2*π*pol_turns/nz0/ι,length=nz0))

pest_points = MagneticCoordinateCurve(PestCoordinates, ψ, α, ζ_range);

cartCoords = CartesianFromVmec()(VmecFromPest()(pest_points, vmec_surf), vmec_surf)

Xf = [q[1] for q in cartCoords];
Yf = [q[2] for q in cartCoords];
Zf = [q[3] for q in cartCoords];

ftmat = [Xf Yf Zf]

writedlm("/u/nigu/gene-stuff/geometry/vmec/3D_vmec_surface_vis/flux_tube_kjm_test_mat.dat", ftmat)	

ft = lines!(Xf, Yf, Zf, color=:red, linewidth=6, label="Flux tube")

# calling specified surface above, named p 
print(p)
p

