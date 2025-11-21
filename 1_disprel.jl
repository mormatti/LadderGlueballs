include("usings.jl")

# Local scope
let

# Parameters to set
reset = false
groups = ["ZZ3", "SU3"]
# λ_list = sort(unique(vcat([1//10, 3//10, 5//10, 7//10, 9//10], [a//30 for a in 1:29])))
λ_list = [1//10, 3//10, 5//10, 7//10, 9//10]
lengths = [5, 7, 9, 11]
nlevels = 50 # should be 30 or more

# Reading data
constants = (JLD2.load("data/constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("data/smallsize.jld2"))["single_stored_object"]

# Group loop
for group in groups
println("Processing group ", group)

smallsize[group] = reset ? Dict() : get!(smallsize, group, Dict())

# Coupling loop
for λ ∈ λ_list
println("Processing λ = ", λ)

smallsize[group][λ] = reset ? Dict() : get!(smallsize[group], λ, Dict())

# Small system size loop
for l ∈ lengths
println("Processing system size L = ", l)

smallsize[group][λ][l] = reset ? Dict() : get!(smallsize[group][λ], l, Dict())
dct = smallsize[group][λ][l]

# Local operators
dct["h"] = λ * constants[group]["Esq3"] - (1 - λ) * (constants[group]["Up3"] + constants[group]["Up3"]')

# Operators
dct["C"] = kron_power(SparseMatrixCSC(constants[group]["C1"]), l)
dct["T"] = operator_translation(SparseMatrixCSC, d, l)
dct["R"] = operator_reflection(SparseMatrixCSC, d, l) # The 180 deg rot, for this model is like that
dct["H"] = summation_local(dct["h"], d, l; pbc = true)
dct["P"] = dct["R"] * dct["C"]

dct["disprel"] = Dict()

# Dispersion relation
drel = dispersion_relation(dct["H"], dct["T"], l, nlevels = nlevels)

# State quantum numbers
dct["disprel"]["E"] = [energy(el) for el ∈ drel]
dct["disprel"]["k"] = [momentum(el) for el ∈ drel]
dct["disprel"]["C"] = [real(wavefunction(el)' * dct["C"] * wavefunction(el)) for el ∈ drel]

# State parity. Zero if not defined (for non-zero and non-pi momentum)
parities = []
for el in drel
    if el.koverpi == 0//1 || el.koverpi == 1//1 || el.koverpi == -1//1
        parity = real(wavefunction(el)' * dct["P"] * wavefunction(el))
    else
        parity = 0
    end
    push!(parities, parity)
end
dct["disprel"]["P"] = parities

# Get band properties
gs = pop_groundstate!(drel)
dct["groundstate"] = gs
dct["gsenergy"] = energy(gs)

# The first band of the first sector (the 0++ band)
band1 = get_firstband(drel)
dct["band1"] = band1

# We get the second band (the 0-- band) deleting all the states with with C > 0 and getting the first band
tempdrel = filter!(x -> real(wavefunction(x)' * dct["C"] * wavefunction(x)) ≤ 0, drel)
band2 = get_firstband(tempdrel)
dct["band2"] = band2

# We compute the average of the energies of the first band
energiesband1 = [energy(el) for el in band1 for _ in 1:((el.koverpi in (0//1, 1//1, -1//1)) ? 1 : 2)]
dct["centroid1"] = sum(energiesband1) / length(energiesband1)
energiesband2 = [energy(el) for el in band2 for _ in 1:((el.koverpi in (0//1, 1//1, -1//1)) ? 1 : 2)]
dct["centroid2"] = sum(energiesband2) / length(energiesband2)

end # system size loop

end # coupling loop

end # group loop

JLD2.save_object("data/smallsize.jld2", smallsize)

end # End local scope