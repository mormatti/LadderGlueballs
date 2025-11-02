group = "SU3"
lambda = 1//10
length = 11
jc = Int((length + 1)/2)
band = "1"
Lbig = 101
k0 = π/2
j0 = Int(round((Lbig + 1)/(3)))
sigma = 3
Lw = size(creator)[1]
maxdim = 150
cutoff = 1e-10

println("Reading data...")
constants = (JLD2.load("glueballs/data/0_constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("glueballs/data/1_smallsize.jld2"))["single_stored_object"]
wannier = (JLD2.load("glueballs/data/2_wannier.jld2"))["single_stored_object"]
biggroundstates = (JLD2.load("glueballs/data/3_biggroundstates.jld2"))["single_stored_object"]
vacua = (JLD2.load("glueballs/data/3_biggroundstates.jld2"))["single_stored_object"]

println("Computing the speed and the simulation time...")
energ = dispersion_interpolant(smallsize[group][lambda][length]["band$band"])
speed = dispersion_interpolant(smallsize[group][lambda][length]["band$band"], derivative = 1)
dispr = dispersion_interpolant(smallsize[group][lambda][length]["band$band"], derivative = 2)
result = optimize(k -> -abs(speed(k)), -0.1, π + 0.1)  # minimize the negative
kspeedmax = Optim.minimizer(result)
vmax = abs(speed(kspeedmax))
println("kspeedmax = ", kspeedmax / π, " π")
println("speedmax = ", vmax, " sites per unit time")


println("Importing the large vacuum...")
sitesL = siteinds(3, Lbig)
vacuumlarge = vacua[group][Lbig][lambda]["mps"]
replace_siteinds!(vacuumlarge, sitesL)

println("Computing the first creator and applying it...")
f(j) = exp(im * -k0 * j) * exp(-(j - (j0 - Lw/2))^2 / (2 * sigma^2))
psi1 = summation_local(creator, Lbig; convolution = f)
truncate!(psi1, cutoff = cutoff, maxdim = maxdim)
replace_siteinds!(psi1, sitesL)
psi = apply(psi1, vacuumlarge)
truncate!(psi, cutoff = cutoff, maxdim = maxdim)
normalize!(psi)

println("Computing the second creator and applying it...")
f(j) = exp(im * k0 * j) * exp(-(j - (j0 + Lbig - 2 * j0 - Lw/2))^2 / (2 * sigma^2))
psi2 = summation_local(creator, Lbig; convolution = f)
truncate!(psi2, cutoff = cutoff, maxdim = maxdim)
replace_siteinds!(psi2, sitesL)
psi = apply(psi2, psi)
truncate!(psi, cutoff = cutoff, maxdim = maxdim)
normalize!(psi)

println("Constructing the Hamiltonian MPO...")
H0 = lambda * constants[group]["Esq3"] - (1 - lambda) * (constants[group]["Up3"] + constants[group]["Up3"]')
H0mpo = mpo_from_matrix(Matrix(H0), d, cutoff, maxdim)
Hmpo = summation_local(H0mpo, Lbig, pbc = false, cutoff = cutoff)
replace_siteinds!(Hmpo, sitesL)

println("Deallocate data...")
constants = 0
smallsize = 0
wannier = 0
biggroundstates = 0
vacua = 0

# Evolve keeping into account the velocity
println("Computing time evolution...")
Δt = Lbig / (2 * vmax)
data = Dict()
gr()
tdvp_time_evolution!(data, Hmpo, psi, Δt / 100, Δt, H0mpo, maxdim = maxdim, cutoff = cutoff, groundstate = vacuumlarge)

# println("Computing local expectation values...")
# densitymat = local_expvals(data[:states], H0mpo)

# Measure the "density" and construct the matrix



println("Done!")

# data = Dict()

# We compute the time evolution and we observe local operators
# states, info = tdvp_time_evolution(H, psi, dt, Dt, maxdim = χmax, cutoff = ϵ)

# tdvp_time_evolution!(data, H, psi, dt, Dt, H0, maxdim = χmax, cutoff = ϵ)

# println("Computing local expectation values...")
# energymatrix = local_expvals(states, H0, d)
# Plots.heatmap(energymatrix)

constants = 0
smallsize = 0
wannier = 0
biggroundstates = 0
vacua = 0