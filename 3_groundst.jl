include("usings.jl")

# Local scope
let

# Parameters to set    
reset = false
Lbigs = [21, 31, 51, 101, 151, 201]
groups = ["ZZ3", "SU3"]
lambdas = [1//10, 3//10, 5//10, 7//10, 9//10]

dmrgargs = (nsweeps = 1000, maxdim = [50,100,200,400,500], cutoff = [1e-10], observer = CustomObserver(1e-10))

# Reading data
constants = (JLD2.load("glueballs/data/0_constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("glueballs/data/1_smallsize.jld2"))["single_stored_object"]
wannier = (JLD2.load("glueballs/data/2_wannier.jld2"))["single_stored_object"]
biggroundstates = (JLD2.load("glueballs/data/3_biggroundstates.jld2"))["single_stored_object"]

biggroundstates = Dict()

for group in groups
biggroundstates[group] = Dict()

# We loop over Ls
for Lbig in Lbigs
biggroundstates[group][Lbig] = Dict()

# We loop over lambdas
for λ in lambdas
biggroundstates[group][Lbig][λ] = Dict()
dict = biggroundstates[group][Lbig][λ]

# We generate the local Hamiltonian in matrix and MPO form
H0 = λ * constants[group]["Esq3"] - (1 - λ) * (constants[group]["Up3"] + constants[group]["Up3"]')
H0mpo = mpo_from_matrix(Matrix(H0), d; cutoff = cutoff, maxdim = maxdim)

# We consider a large system
sites = siteinds(d, Lbig)

# We construct the large Hamiltonian
Hmpo = summation_local(H0mpo, Lbig, pbc = false, cutoff = cutoff, maxdim = maxdim)
replace_siteinds!(Hmpo, sites)

# We perform DMRG
psir = random_mps(sites)
dict["energy"], dict["mps"] = dmrg(Hmpo, psir; dmrgargs...)

end # End λ loop

end # End L loop

end # End group loop

JLD2.save_object("glueballs/data/3_biggroundstates.jld2", biggroundstates)

end # End local scope