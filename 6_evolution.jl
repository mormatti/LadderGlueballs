include("usings.jl")

let

println("Reading data...")
constants = (JLD2.load("data/constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("data/smallsize.jld2"))["single_stored_object"]
wannier = (JLD2.load("data/wannier.jld2"))["single_stored_object"]
vacua = (JLD2.load("data/biggroundstates.jld2"))["single_stored_object"]

# Parameters to choose
G = "SU3"
lambda = 1//10 # put 1//10 then
l = 11
B = "1"
L = 101 # put 100 then

savename = "data/evolution_pp_G$(G)_la$(denominator(lambda))_B$(B)_L$(L).jld2"

jc = Int((l + 1)/2)
j0_1 = Int(round((L + 1)/(3)))
j0_2 = L - j0_1
sigma = L / 30

maxdim = 100 # put 100 then
cutoff = 1e-10 # put the right value then

println("Importing the creator and truncating it...")
creator = (JLD2.load("data/creators.jld2"))["single_stored_object"][G][lambda][l][B][2]["creator"]
truncate!(creator, maxdim = maxdim, cutoff = cutoff)
Lw = size(creator)[1]

println("Computing the speed and the simulation time...")
energ = dispersion_interpolant(smallsize[G][lambda][l]["band$B"])
speed = dispersion_interpolant(smallsize[G][lambda][l]["band$B"], derivative = 1)
dispr = dispersion_interpolant(smallsize[G][lambda][l]["band$B"], derivative = 2)
result = optimize(k -> -abs(speed(k)), -0.1, π + 0.1)  # minimize the negative
k0 = Optim.minimizer(result)
vmax = abs(speed(k0))
println("kspeedmax = ", k0 / π, " π")
println("speedmax = ", vmax, " sites per unit time")

println("Importing the large vacuum...")
sitesL = siteinds(3, L)
vacuumlarge = vacua[G][L][lambda]["mps"]
replace_siteinds!(vacuumlarge, sitesL)

gauss_profile(j, km, jm) = exp(im * km * j) * exp(-(j - jm + (Lw - 1)/2)^2 / (2 * sigma^2))

println("Computing the first creator and applying it...")
psi1 = summation_local(creator, L, convolution = j -> gauss_profile(j, -k0, j0_1), cutoff = cutoff, maxdim = maxdim)
replace_siteinds!(psi1, sitesL)
psi = normalize(apply(psi1, vacuumlarge, cutoff = cutoff, maxdim = maxdim))

println("Computing the second creator and applying it...")
psi2 = summation_local(creator, L, convolution = j -> gauss_profile(j, k0, j0_2), cutoff = cutoff, maxdim = maxdim)
replace_siteinds!(psi2, sitesL)
psi = normalize(apply(psi2, psi, cutoff = cutoff, maxdim = maxdim))

println("Constructing the Hamiltonian MPO...")
H0 = lambda * constants[G]["Esq3"] - (1 - lambda) * (constants[G]["Up3"] + constants[G]["Up3"]')
H0mpo = mpo_from_matrix(Matrix(H0), d, cutoff, maxdim)
Hmpo = summation_local(H0mpo, L, pbc = false, cutoff = cutoff)
replace_siteinds!(Hmpo, sitesL)

println("Deallocate data...")
constants = 0
smallsize = 0
wannier = 0
vacua = 0
creator = 0

# Evolve keeping into account the velocity
println("Computing time evolution...")
Δt = L / (2 * vmax) * 1.5
println("Total time Δt = ", Δt)
data = Dict()
gr()
tdvp_time_evolution!(data, Hmpo, psi, Δt / 150, Δt, H0mpo, maxdim = maxdim, cutoff = cutoff, groundstate = vacuumlarge)

# println("Computing local expectation values...")
# densitymat = local_expvals(data[:states], H0mpo)

# Measure the "density" and construct the matrix
println("Done!")

JLD2.save_object(savename, data)

# data = Dict()

# We compute the time evolution and we observe local operators
# states, info = tdvp_time_evolution(H, psi, dt, Dt, maxdim = χmax, cutoff = ϵ)

# tdvp_time_evolution!(data, H, psi, dt, Dt, H0, maxdim = χmax, cutoff = ϵ)

# println("Computing local expectation values...")
# energymatrix = local_expvals(states, H0, d)
# Plots.heatmap(energymatrix)

end