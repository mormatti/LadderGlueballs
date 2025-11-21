include("usings.jl")

let

reset = false
Lbig = 101
groups = ["ZZ3", "SU3"]
lambdas = [9//10]
lengths = [11]
bands = ["1"]
supporthalfs = [2]
cutoff = 1e-15
maxdim = 100

# Reading data
constants = (JLD2.load("data/constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("data/smallsize.jld2"))["single_stored_object"]
wannier = (JLD2.load("data/wannier.jld2"))["single_stored_object"]
vacua = (JLD2.load("data/biggroundstates.jld2"))["single_stored_object"]
creators = (JLD2.load("data/creators.jld2"))["single_stored_object"]
overlaps = (JLD2.load("data/overlaps.jld2"))["single_stored_object"]

for G in groups
overlaps[G] = reset ? Dict() : get!(overlaps, G, Dict())
println("Group G = ", G)
println()

for λ ∈ lambdas
overlaps[G][λ] = reset ? Dict() : get!(overlaps[G], λ, Dict())
println("Coupling λ = ", λ)
println()

sitesL = siteinds(3, Lbig)
vacuumlarge = vacua[G][Lbig][λ]["mps"]
replace_siteinds!(vacuumlarge, sitesL)
H0 = λ * constants[G]["Esq3"] - (1 - λ) * (constants[G]["Up3"] + constants[G]["Up3"]')
H0mpo = mpo_from_matrix(Matrix(H0), d, cutoff, maxdim)
Hmpo = summation_local(H0mpo, Lbig, pbc = false, cutoff = cutoff)
replace_siteinds!(Hmpo, sitesL)

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
overlaps[G][λ][l] = reset ? Dict() : get!(overlaps[G][λ], l, Dict())
println("Chain length = ", l)
println()
jc = Int((l + 1)/2)
Jc = Int((Lbig + 1)/2)

for B in bands
overlaps[G][λ][l][B] = reset ? Dict() : get!(overlaps[G][λ][l], B, Dict())
println("Band number = ", B)
println()

for ℓ in supporthalfs
overlaps[G][λ][l][B][ℓ] = reset ? Dict() : get!(overlaps[G][λ][l][B], ℓ, Dict())
thedict = overlaps[G][λ][l][B][ℓ]
println("Halfsupport = ", ℓ)

println("")
println("Group: ", G, " λ = ", λ, " band: ", B)
creator = (JLD2.load("data/creators.jld2"))["single_stored_object"][G][λ][l][B][ℓ]["creator"]
creatorU = creator
# sgn = (B == "1" ? 1 : -1)
# creator = mpo_from_matrix(Matrix(constants[G]["Up3"] + sgn * constants[G]["Up3"]') * 1/sqrt(2), 3, cutoff, maxdim)

eps = (λ - 1)/λ
if G == "ZZ3" && B == "1"
    Ebar = 27 / π^2 + eps + 4 * π^2/27 * eps^2
    vmaxpert = 8 * π^2 / 81 * λ*eps^2
    thedict["tperturbative"] = [-8 * π^2 / 81 * eps^2, 0, 0, 0]
elseif G == "ZZ3" && B == "2"
    Ebar = 27 / π^2
    vmaxpert = 4 * π^2 / 81* λ*eps^2
    thedict["tperturbative"] = [-4 * π^2 / 81 * eps^2, 0, 0, 0]
elseif G == "SU3" && B == "1"
    Ebar = 16/3 + eps + 3/4 * eps^2
    vmaxpert = -1/2 * λ*eps^2
    thedict["tperturbative"] = [eps^2 / 2, 0, 0, 0]
elseif G == "SU3" && B == "2"
    Ebar = 16/3
    vmaxpert = -1/12 * λ*eps^2
    thedict["tperturbative"] = [eps^2 / 12, 0, 0, 0]
end

energ = dispersion_interpolant(smallsize[G][λ][l]["band$B"])
speed = dispersion_interpolant(smallsize[G][λ][l]["band$B"], derivative = 1)
dispr = dispersion_interpolant(smallsize[G][λ][l]["band$B"], derivative = 2)
result = optimize(k -> -abs(speed(k)), -0.1, π + 0.1)  # minimize the negative
kspeedmax = Optim.minimizer(result)
vmaxreal = speed(kspeedmax)
thedict["disprelspeedmax"] = vmaxreal

function creatoratsite(ctr, j::Int)
    creatorlength = size(ctr)[1]
    halfcs = Int(floor(creatorlength / 2))
    minj = halfcs + 1
    maxj = Lbig - halfcs
    if j < minj || j > maxj
        error("Position $j is out of bounds for creator of length $creatorlength in lattice
    of length $Lbig")
    end
    return insert_local(j - halfcs - 1, ctr, Lbig - creatorlength - (j - halfcs - 1))
end

function wannieratsite(ctr, j::Int)
    creatormpo = creatoratsite(ctr, j)
    replace_siteinds!(creatormpo, sitesL)
    wannierbig = apply(creatormpo, vacuumlarge)
    normalize!(wannierbig)
    truncate!(wannierbig, cutoff = cutoff, maxdim = maxdim)
    replace_siteinds!(wannierbig, sitesL)
    return wannierbig
end

thedict["gram"] = []
thedict["gramU"] = []
thedict["tcreator"] = []
thedict["tcreatorU"] = []

wannierzero = wannieratsite(creator, Jc)
wannierUzero = wannieratsite(creatorU, Jc)
for jp in 1:5
    wannierjp = wannieratsite(creator, Jc + jp)
    wannierUjp = wannieratsite(creatorU, Jc + jp)
    push!(thedict["gram"], inner(wannierjp, wannierzero))
    push!(thedict["gramU"], inner(wannierUjp, wannierUzero))
    push!(thedict["tcreator"], inner(wannierjp', Hmpo, wannierzero))
    push!(thedict["tcreatorU"], inner(wannierUjp', Hmpo, wannierUzero))
    # println("Epsilon = ", eps)
    # println("Ebar = ", Ebar, ", t0 = ", ovlpt0)
    # println("⟨j|j⟩ = ", gram0, ", ⟨j+1|j⟩ = ", gram1)
    # println("vmaxpert = ", vmaxpert, ", t1 = ", -2 * ovlpt1, ", vmaxreal = ", vmaxreal)
    # println("relative error = ", ((-2 * ovlpt1) - vmaxreal) / vmaxreal)
end

end # end supporthalfs

end # end bands

end # end lengths

end # end lambdas

end # end groups

JLD2.save_object("data/overlaps.jld2", overlaps)

end # end local scope