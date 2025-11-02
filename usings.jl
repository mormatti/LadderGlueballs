using Revise
using Scattensor
using ITensors, ITensorMPS
using LinearAlgebra
using SparseArrays
using Optim
using KrylovKit
using Logging
using JLD2
using Measures

using PGFPlotsX
PGFPlotsX.latexengine!(PGFPlotsX.PDFLATEX)  # or PGFPlotsX.LUALATEX
PGFPlotsX.latexengine()                     # shows the current engine
PGFPlotsX.CUSTOM_PREAMBLE = [
    raw"\usepackage{amsfonts}",  # provides \mathbb via amsfonts
    raw"\usepackage{amssymb}",   # also ok for \mathbb
    raw"\pgfplotsset{every axis legend/.style={",
    raw"    cells={anchor=west},",
    raw"    legend columns=-1,",
    raw"    at={(0.5,-0.2)},",
    raw"    anchor=north,",
    raw"    rounded corners=3pt,",
    raw"    draw=gray!60,",
    raw"    very thin,",
    raw"    fill=white,",
    raw"    font=\small",
    raw"}}"
]

using Plots, LaTeXStrings
pgfplotsx()

# Global parameters
d = 3 # The local dimension
cutoff = 1e-10 # Cutoff for everything
maxdim = 200 # Maximum svd truncation dimension in MPS

function ketbra(m::Integer, n::Integer, q::Integer)
    M = zeros(Int64, q, q)
    M[m, n] = 1
    return SparseMatrixCSC(M)
end

function g_from_lambda(λ::Real)
    if λ < 0 || λ >= 1
        error("Invalid λ value")
    end
    return (λ/(1 - λ))^(1/4)
end

function color_from_lambda(λ::Real)
    pE = [0.2,0.3,0.7]
    pB = [0.4,0.7,0.8]
    pI = λ * pE + (1 - λ) * pB
    return RGB(pI...)
end

A = operator_identity(MPO, 3, 3)