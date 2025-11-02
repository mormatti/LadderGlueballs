include("usings.jl")

# Local scope
let

# Parameters to set
reset = false
groups = ["ZZ3", "SU3"]
lambdas = [a//30 for a in 1:29]
lengths = [11]

# Reading data
constants = (JLD2.load("glueballs/data/0_constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("glueballs/data/1_smallsize.jld2"))["single_stored_object"]
wannier = Dict()
# wannier = (JLD2.load("glueballs/data/2_wannier.jld2"))["single_stored_object"]

for group in groups
wannier[group] = reset ? Dict() : get!(wannier, group, Dict())

for λ ∈ lambdas
wannier[group][λ] = reset ? Dict() : get!(wannier[group], λ, Dict())

for Li in lengths
wannier[group][λ][Li] = reset ? Dict() : get!(wannier[group][λ], Li, Dict())

dct1 = smallsize[group][λ][Li]
dct2 = wannier[group][λ][Li]

if Li % 2 == 0
    error("Even length not supported yet for wannier minimization")
end
dct1["jc"] = Int(floor(7/2 + 1))

# Defining the local hamiltonian centered in the reflection site
IdmLat = operator_identity(SparseMatrixCSC, 3^Int((Li - 3)/2))
dct1["Hc"] = kron(IdmLat, dct1["h"], IdmLat)

dct2["wannier1"], dct2["info1"] = wannier_symmetric(
                        dct1["band1"],
                        dct1["gsenergy"],
                        x -> dct1["Hc"] * x,
                        x -> dct1["P"] * x,
                        (x,y) -> x' * y
                        )

dct2["wannier2"], dct2["info2"] = wannier_symmetric(
                        dct1["band2"],
                        dct1["gsenergy"],
                        x -> dct1["Hc"] * x,
                        x -> dct1["P"] * x,
                        (x,y) -> x' * y
                        )

end # end length loop

end # end lambda loop

end # end group loop

JLD2.save_object("glueballs/data/1_smallsize.jld2", smallsize)
JLD2.save_object("glueballs/data/2_wannier.jld2", wannier)

end # end local scope