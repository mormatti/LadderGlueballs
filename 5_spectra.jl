group = "SU3"
lambda = 1//10
length = 11
jc = Int((length + 1)/2)
band = "1"
Lbig = 51
Jc = Int((Lbig + 1)/2)
creatorlength = size(creator)[1]
halfcs = Int(floor(creatorlength / 2))
minx = halfcs + 1
maxx = Lbig - halfcs
maxdim = 20

constants = (JLD2.load("glueballs/data/0_constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("glueballs/data/1_smallsize.jld2"))["single_stored_object"]
wannier = (JLD2.load("glueballs/data/2_wannier.jld2"))["single_stored_object"]
biggroundstates = (JLD2.load("glueballs/data/3_biggroundstates.jld2"))["single_stored_object"]
vacua = (JLD2.load("glueballs/data/3_biggroundstates.jld2"))["single_stored_object"]

println("Computing point of maximum speed...")
energ = dispersion_interpolant(smallsize[group][lambda][length]["band$band"])
speed = dispersion_interpolant(smallsize[group][lambda][length]["band$band"], derivative = 1)
dispr = dispersion_interpolant(smallsize[group][lambda][length]["band$band"], derivative = 2)
result = optimize(k -> -abs(speed(k)), -0.1, π + 0.1)  # minimize the negative
kspeedmax = Optim.minimizer(result)
vmax = abs(speed(kspeedmax))
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
vacuumlarge = vacua[group][Lbig][lambda]["mps"]
replace_siteinds!(vacuumlarge, sitesL)
wannierbig = normalize(apply(creatormpo, vacuumlarge))
replace_siteinds!(wannierbig, sitesL)

println("Constructing the Hamiltonian MPO...")
H0 = lambda * constants[group]["Esq3"] - (1 - lambda) * (constants[group]["Up3"] + constants[group]["Up3"]')
H0mpo = mpo_from_matrix(Matrix(H0), d, cutoff, maxdim)
Hmpo = summation_local(H0mpo, Lbig, pbc = false, cutoff = cutoff)
replace_siteinds!(Hmpo, sitesL)

# Evolve keeping into account the velocity
println("Computing time evolution...")
Δt = Lbig / (2 * vmax) * 0.7
steps = 50
dt = Δt / steps
data = Dict()
println("First simulation...")
tdvp_time_evolution!(data, Hmpo, wannierbig, dt, Δt, H0mpo, maxdim = maxdim, cutoff = cutoff, groundstate = vacuumlarge)
# println("Computing local expectation values...")
# densitymat = local_expvals(data[:states], H0mpo)

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

constants = 0
smallsize = 0
wannier = 0
biggroundstates = 0
vacua = 0

import Pkg
using FFTW

gr()
Plots.heatmap(abs.(fftshift(fft(corrfunc), 1)))

# Measure the "density" and construct the matrix

# In a second file plot the fourier transform and overimpose the dispersion relation