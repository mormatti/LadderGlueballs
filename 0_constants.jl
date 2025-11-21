include("usings.jl")

let

# Defiing the local operators and constants
constants = Dict()
constants["ZZ3"] = Dict()
constants["SU3"] = Dict()
ZZ3 = constants["ZZ3"]
SU3 = constants["SU3"]

# Number of states in the local Hilbert space
ZZ3["localdim"] = 3
SU3["localdim"] = 3

# The identity operators
ZZ3["idm"] = operator_identity(SparseMatrixCSC, 3)
SU3["idm"] = operator_identity(SparseMatrixCSC, 3)
ZZ3["I1"] = ZZ3["idm"]
SU3["I1"] = SU3["idm"]
ZZ3["I2"] = kron(ZZ3["idm"], ZZ3["idm"])
SU3["I2"] = kron(SU3["idm"], SU3["idm"])
ZZ3["I3"] = kron(ZZ3["idm"], ZZ3["idm"], ZZ3["idm"])
SU3["I3"] = kron(SU3["idm"], SU3["idm"], SU3["idm"])

# The reflection operators
ZZ3["R2"] = operator_reflection(SparseMatrixCSC, 3, 2)
SU3["R2"] = operator_reflection(SparseMatrixCSC, 3, 2)
ZZ3["R3"] = operator_reflection(SparseMatrixCSC, 3, 3)
SU3["R3"] = operator_reflection(SparseMatrixCSC, 3, 3)

# The raising operator tau
ZZ3["tau"] = spzeros(3, 3)
ZZ3["tau"][1, 2] = 1
ZZ3["tau"][2, 3] = 1
ZZ3["tau"][3, 1] = 1
SU3["tau"] = ZZ3["tau"]

# The charge conjugation operator C
cc = spzeros(3, 3)
cc[1, 1] = 1
cc[2, 3] = 1
cc[3, 2] = 1
ZZ3["C1"] = cc
ZZ3["C2"] = kron(cc, cc)
ZZ3["C3"] = kron(cc, cc, cc)
SU3["C1"] = cc
SU3["C2"] = kron(cc, cc)
SU3["C3"] = kron(cc, cc, cc)

# A shortcut for ketbra
kb(n, m) = ketbra(n, m, 3)

# diagmat (old version) = SparseMatrixCSC(Diagonal([0.0, 1.0, 1.0, 2.0, 2.0, 2.5, 2.0, 2.5, 2.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.5, 2.5, 3.0, 2.5, 1.0, 2.0, 2.0, 2.5, 2.5, 3.0, 2.0, 2.5, 2.0]))
diagmat = SparseMatrixCSC(Diagonal([0, 0.5, 0.5, 3, 2.5, 3, 3, 3, 2.5, 0.5, 1, 1, 2.5, 2, 2.5, 3, 3, 2.5, 0.5, 1, 1, 3, 2.5, 3, 2.5, 2.5, 2]))
ZZ3["Esq3"] = 27/(4*π^2) * diagmat
SU3["Esq3"] = 4//3 * diagmat

# Plaquette term for Z_3
ZZ3["Up"] = ZZ3["tau"]
ZZ3["Up3"] = kron(ZZ3["I1"], ZZ3["tau"], ZZ3["I1"])

# Plaquette term for SU(3)
Dg1 = Diagonal([1, 1/√3, 1/3])
Dg2 = Diagonal([1, 1/3, 1/√3])
Dg3 = Diagonal([1, 1/√3, 1/√3])
Up = kron(Dg1, kb(3,1), Dg1) + kron(Dg2, kb(1,2), Dg2) + kron(Dg3, kb(2,3), Dg3)
SU3["Up"] = Up
SU3["Up3"] = Up

JLD2.save_object("data/constants.jld2", constants)

end