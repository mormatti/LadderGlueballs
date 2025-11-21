include("usings.jl")

let

# Parameters to set
reset = false
groups = ["ZZ3", "SU3"]
lambdas = [1//10, 5//10, 9//10]
lengths = [11]
bands = ["1", "2"]
supporthalfs = [0, 1, 2]
cutoff = 1e-15
maxdim = 100

# Reading data
# constants = (JLD2.load("data/constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("data/smallsize.jld2"))["single_stored_object"]
wannier = (JLD2.load("data/wannier.jld2"))["single_stored_object"]
# biggroundstates = (JLD2.load("data/biggroundstates.jld2"))["single_stored_object"]
creators = (JLD2.load("data/creators.jld2"))["single_stored_object"]
# creators = Dict()

for G in groups
creators[G] = reset ? Dict() : get!(creators, G, Dict())
println("Group G = ", G)
println()

for λ ∈ lambdas
creators[G][λ] = reset ? Dict() : get!(creators[G], λ, Dict())
println("Coupling λ = ", λ)
println()

#=
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
=#

for l in lengths
creators[G][λ][l] = reset ? Dict() : get!(creators[G][λ], l, Dict())
println("Chain length = ", l)
println()

for B in bands
creators[G][λ][l][B] = reset ? Dict() : get!(creators[G][λ][l], B, Dict())
println("Band number = ", B)
println()

for ℓ in supporthalfs
creators[G][λ][l][B][ℓ] = reset ? Dict() : get!(creators[G][λ][l][B], ℓ, Dict())
thedict = creators[G][λ][l][B][ℓ]
println("Halfsupport = ", ℓ)

jc = Int((l + 1)/2)

# Pick the small groundstate and wannier and construct the matrix and mpo
sitesl = siteinds(3, l)
ωvector = wavefunction(smallsize[G][λ][l]["groundstate"])
ωmps = mps_from_vector(ωvector, 3, cutoff = cutoff)
replace_siteinds!(ωmps, sitesl)
wvector = wannier[G][λ][l]["wannier$B"]
wmps = mps_from_vector(wvector, 3, cutoff = cutoff)
replace_siteinds!(wmps, sitesl)

# We define the A and the rho operators in matrix and mpo form
Ampo = partial_trace(truncate(outer(wmps', ωmps), cutoff = cutoff, maxdim = maxdim), jc - ℓ, jc + ℓ)
ρmpo  = partial_trace(truncate(outer(ωmps', ωmps), cutoff = cutoff, maxdim = maxdim), jc - ℓ, jc + ℓ)
Amatrix = matrix_from_mpo(Ampo)
ρmatrix = matrix_from_mpo(ρmpo)

# Method 1: not imposing unitarity
creatormatrix = Amatrix * pinv(ρmatrix, 1e-10)
thedict["creatormatrix"] = creatormatrix
creatormpo = mpo_from_matrix(creatormatrix, d, cutoff, maxdim)
truncate!(creatormpo, cutoff = cutoff, maxdim = maxdim)
thedict["creatormpo"] = creatormpo
# tests
latspace = Int((l - (2ℓ+1))/2)
extendedcreatormpo = insert_local(latspace, creatormpo, latspace)
replace_siteinds!(extendedcreatormpo, sitesl)
createdwannier = apply(extendedcreatormpo, ωmps)
normcreated = norm(createdwannier)
println("Norm of the created wannier = ", normcreated)
creatormatrix /= normcreated
creatormpo /= normcreated
extendedcreatormpo /= normcreated
fidelity = abs(inner(wmps', extendedcreatormpo, ωmps))
println("Fidelity = ", fidelity)
thedict["fidelity"] = fidelity

# Method 2: imposing unitarity
Svd = svd(matrix_from_mpo(Ampo))
creatormatrix = Svd.U * Svd.Vt
thedict["unitarycreatormatrix"] = creatormatrix
creatormpo = mpo_from_matrix(creatormatrix, 3, cutoff, maxdim)
truncate!(creatormpo, cutoff = cutoff, maxdim = maxdim)
thedict["unitarycreatormpo"] = creatormpo
# tests
svdsum = sum(Svd.S)
println("Sum of svd = ", svdsum)
thedict["svdsum"] = svdsum
println("1 - sum of svd = ", 1 - svdsum)
latspace = Int((l - (2ℓ+1))/2)
extendedcreatormpo = insert_local(latspace, creatormpo, latspace)
replace_siteinds!(extendedcreatormpo, sitesl)
createdwannier = apply(extendedcreatormpo, ωmps)
normcreated = norm(createdwannier)
println("Norm of the created wannier = ", normcreated)
creatormatrix /= normcreated
creatormpo /= normcreated
extendedcreatormpo /= normcreated
fidelity = abs(inner(wmps', extendedcreatormpo, ωmps))
println("Fidelity = ", fidelity)
thedict["unitaryfidelity"] = fidelity

# Creating the projectors of vacuum and of the particle
Pmpo = partial_trace(truncate(outer(ωmps', ωmps), cutoff = cutoff, maxdim = maxdim), jc - ℓ, jc + ℓ)
Pmatrix = matrix_from_mpo(Pmpo)
truncate!(Pmpo, cutoff = cutoff)
thedict["vacuumproj"] = Pmpo
Nmpo = partial_trace(truncate(outer(wmps', wmps), cutoff = cutoff, maxdim = maxdim), jc - ℓ, jc + ℓ)
Nmatrix = matrix_from_mpo(Nmpo)
truncate!(Nmpo, cutoff = cutoff)
thedict["numberop"] = Nmpo

end # end supporthalf loop

end # end band loop

end # end length loop

end # end lambda loop

end # end group loop

JLD2.save_object("data/creators.jld2", creators)

end # end local scope
