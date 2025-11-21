include("usings.jl")

# Local scope
let

# Parameters to set
reset = false
Lbigs = [21, 31, 51, 101, 201]
groups = ["ZZ3", "SU3"]
lambdas = [1//10, 5//10, 9//10]

# [25,50,100,200,400,500]

dmrgargs = (nsweeps = 1000, maxdim = [25,50,100,200,400,500], cutoff = [cutoff]) # observer = CustomObserver(1e-15))

# Reading data
constants = (JLD2.load("data/constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("data/smallsize.jld2"))["single_stored_object"]
wannier = (JLD2.load("data/wannier.jld2"))["single_stored_object"]
biggroundstates = (JLD2.load("data/biggroundstates.jld2"))["single_stored_object"]

biggroundstates = Dict()

for group in groups
biggroundstates[group] = reset ? Dict() : get!(biggroundstates, group, Dict())
print("G = ", group)

# We loop over Ls
for Lbig in Lbigs
biggroundstates[group][Lbig] = reset ? Dict() : get!(biggroundstates[group], Lbig, Dict())
print("L = ", Lbig)

# We loop over lambdas
for λ in lambdas
biggroundstates[group][Lbig][λ] = reset ? Dict() : get!(biggroundstates[group][Lbig], λ, Dict())
dict = biggroundstates[group][Lbig][λ]
print("λ = ", λ)

if λ == 9//10
    cutoff = 1e-15
    maxdim = 50
end
if λ == 5//10
    cutoff = 1e-10
    maxdim = 70
end
if λ == 1//10
    cutoff = 1e-10
    maxdim = 100
end

# We generate the local Hamiltonian in matrix and MPO form
H0 = λ * constants[group]["Esq3"] - (1 - λ) * (constants[group]["Up3"] + constants[group]["Up3"]')
H0mpo = mpo_from_matrix(Matrix(H0), d, cutoff, maxdim)

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

JLD2.save_object("data/biggroundstates.jld2", biggroundstates)

end # End local scope