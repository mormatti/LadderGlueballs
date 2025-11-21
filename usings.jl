# Basic packages
using LinearAlgebra
using SparseArrays
using Optim
using JLD2
using Measures
using ITensors, ITensorMPS
using LaTeXStrings

# For plotting
using Plots
using PGFPlotsX

# Scattensor and Revise for development
using Revise
using Scattensor

using Base.Threads
println("Julia nthreads()         = ", nthreads())
println("BLAS.get_num_threads()  = ", BLAS.get_num_threads())
# opzionale ma utile: sincronizza BLAS con i thread Julia
BLAS.set_num_threads(nthreads())
println("BLAS threads after set  = ", BLAS.get_num_threads())

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

function bandcolor(G, B)
    if G == "ZZ3" && B == "1"
        color = colorant"#a00000"
    elseif G == "ZZ3" && B == "2"
        color = colorant"#800074"
    elseif G == "SU3" && B == "1"
        color = colorant"#1a80bb"
    elseif G == "SU3" && B == "2"
        color = colorant"#298c8c"
    else
        color = colorant"black"
    end
end

function cmaptransparent(color)
    cmap = cgrad([
        RGBA(red(color), green(color), blue(color), 0.0),
        RGBA(red(color), green(color), blue(color), 1.0)
    ])
    return cmap
end