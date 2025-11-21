include("usings.jl")

let

reset = false
Lbig = 51
groups = ["SU3"]
lambdas = [9//10]
lengths = [11]
bands = ["1", "2"]
supporthalfs = [2]

# Reading data
constants = (JLD2.load("data/constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("data/smallsize.jld2"))["single_stored_object"]
wannier = (JLD2.load("data/wannier.jld2"))["single_stored_object"]
vacua = (JLD2.load("data/biggroundstates.jld2"))["single_stored_object"]
creators = (JLD2.load("data/creators.jld2"))["single_stored_object"]
cones = (JLD2.load("data/conesL$Lbig.jld2"))["single_stored_object"]

for G in groups
cones[G] = reset ? Dict() : get!(cones, G, Dict())
println("Group: ", G)

for λ ∈ lambdas
cones[G][λ] = reset ? Dict() : get!(cones[G], λ, Dict())
println("   Lambda: ", λ)

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

for l in lengths
cones[G][λ][l] = reset ? Dict() : get!(cones[G][λ], l, Dict())

jc = Int((l + 1)/2)
Jc = Int((Lbig + 1)/2)

for B in bands
cones[G][λ][l][B] = reset ? Dict() : get!(cones[G][λ][l], B, Dict())
println("       Band: ", B)

for ℓ in supporthalfs
cones[G][λ][l][B][ℓ] = reset ? Dict() : get!(cones[G][λ][l][B], ℓ, Dict())
println("           Half-support: ", ℓ)

thisdict = cones[G][λ][l][B][ℓ]

creator = (JLD2.load("data/creators.jld2"))["single_stored_object"][G][λ][l][B][ℓ]["creator"]
creatorlength = size(creator)[1]
halfcs = Int(floor(creatorlength / 2))
minx = halfcs + 1
maxx = Lbig - halfcs

println("Computing point of maximum speed...")
energ = dispersion_interpolant(smallsize[G][λ][l]["band$B"])
speed = dispersion_interpolant(smallsize[G][λ][l]["band$B"], derivative = 1)
dispr = dispersion_interpolant(smallsize[G][λ][l]["band$B"], derivative = 2)
result = optimize(k -> -abs(speed(k)), -0.1, π + 0.1)  # minimize the negative
kspeedmax = Optim.minimizer(result)
thisdict["kspeedmax"] = kspeedmax
vmax = abs(speed(kspeedmax))
thisdict["speedmax"] = vmax
println("kspeedmax = ", kspeedmax / π, " π")
println("speedmax = ", vmax, " sites per unit time")

function creatoratsite(x::Int)
    if x < minx || x > maxx
        error("Position $x is out of bounds for creator of length $creatorlength in lattice
    of length $Lbig")
    end
    return insert_local(x - halfcs - 1, creator, Lbig - creatorlength - (x - halfcs - 1))
end

println("Applying the operator to large vacuum...")
sitesL = siteinds(3, Lbig)
creatormpo = creatoratsite(Jc)
replace_siteinds!(creatormpo, sitesL)
vacuumlarge = vacua[G][Lbig][λ]["mps"]
replace_siteinds!(vacuumlarge, sitesL)
wannierbig = normalize(apply(creatormpo, vacuumlarge))
truncate!(wannierbig, cutoff = cutoff, maxdim = maxdim)
replace_siteinds!(wannierbig, sitesL)

println("Constructing the Hamiltonian MPO...")
H0 = λ * constants[G]["Esq3"] - (1 - λ) * (constants[G]["Up3"] + constants[G]["Up3"]')
H0mpo = mpo_from_matrix(Matrix(H0), d, cutoff, maxdim)
Hmpo = summation_local(H0mpo, Lbig, pbc = false, cutoff = cutoff)
replace_siteinds!(Hmpo, sitesL)

# Evolve keeping into account the velocity
println("Computing time evolution...")
Δt = Lbig / (2 * vmax)
println("Δt = ", Δt)
thisdict["evolutiontime"] = Δt
steps = 100
dt = Δt / steps
println("dt = ", dt)
thisdict["timestep"] = dt
data = Dict()
println("Performing the simulation...")
tdvp_time_evolution!(data, Hmpo, wannierbig, dt, Δt, H0mpo, maxdim = maxdim, cutoff = cutoff, groundstate = vacuumlarge)
# println("Computing local expectation values...")
# densitymat = local_expvals(data[:states], H0mpo)
thisdict["evolutiondata"] = data

states = data[:states]
times = data[:times]

corrfunc = zeros(ComplexF64, size(times)[1], Lbig)
for n in eachindex(states)
    for x in minx:maxx
        creatoratx = creatoratsite(x)
        replace_siteinds!(creatoratx, sitesL)
        tooverlap = apply(creatoratx, vacuumlarge)
        corrfunc[n, x] = inner(states[n], tooverlap)
    end
end

thisdict["correlationfunction"] = corrfunc

end # end supporthalfs

end # end bands

end # end lengths

end # end lambdas

end # end groups

JLD2.save_object("data/conesL$Lbig.jld2", cones)

end # end local scope