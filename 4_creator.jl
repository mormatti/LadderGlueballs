include("usings.jl")

# Parameters to set
reset = false
group = "SU3"
lambda = 1//10
length = 11
jc = Int((length + 1)/2)
supporthalf = 2
band = "1"

# Reading data
# constants = (JLD2.load("glueballs/data/0_constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("glueballs/data/1_smallsize.jld2"))["single_stored_object"]
wannier = (JLD2.load("glueballs/data/2_wannier.jld2"))["single_stored_object"]
# biggroundstates = (JLD2.load("glueballs/data/3_biggroundstates.jld2"))["single_stored_object"]

# Pick the small groundstate and wannier and construct the operator (in some way)
sitesl = siteinds(3, length)

ωvector = wavefunction(smallsize[group][lambda][length]["groundstate"])
ωmps = mps_from_vector(ωvector, 3, cutoff = cutoff)
replace_siteinds!(ωmps, sitesl)

wvector = wannier[group][lambda][length]["wannier$band"]
wmps = mps_from_vector(wvector, 3, cutoff = cutoff)
replace_siteinds!(wmps, sitesl)

ρωw = partial_trace(truncate(outer(ωmps', wmps), cutoff = cutoff), jc - supporthalf, jc + supporthalf)
Svd = svd(matrix_from_mpo(ρωw))
precision = 1 - sum(Svd.S)
println("Precision = 1 - Nuclear norm = ", precision)
println("Singular values greater than precision:")
tokeep = filter(x -> (x >= precision), Svd.S)
println(tokeep)

creatormatrix = sum((Svd.V)[:,n] * (Svd.U)[:,n]' for n ∈ eachindex(tokeep))
creator = mpo_from_matrix(creatormatrix, 3, 1e-10, 100)
truncate!(creator, cutoff = cutoff)

Nmpo = partial_trace(truncate(outer(wmps', wmps), cutoff = cutoff), jc - supporthalf, jc + supporthalf)
truncate!(Nmpo, cutoff = cutoff)