include("usings.jl")
using Pkg
using FFTW
gr()

l = 11
L = 51
ℓ = 2

smallsize = (JLD2.load("data/smallsize.jld2"))["single_stored_object"]
cones = (JLD2.load("data/conesL$L.jld2"))["single_stored_object"]
overlaps = (JLD2.load("data/overlaps.jld2"))["single_stored_object"]

function zero_pad_center(A::AbstractArray, newsize::Tuple{Int,Int})
    oldsize = size(A)
    Nx_old, Ny_old = oldsize
    Nx_new, Ny_new = newsize
    @assert Nx_new ≥ Nx_old "newsize[1] must be ≥ size(A,1)"
    @assert Ny_new ≥ Ny_old "newsize[2] must be ≥ size(A,2)"
    B = zeros(eltype(A), Nx_new, Ny_new)
    x0 = div(Nx_new - Nx_old, 2) + 1
    y0 = div(Ny_new - Ny_old, 2) + 1
    B[x0:x0+Nx_old-1, y0:y0+Ny_old-1] .= A
    return B
end

function disprelfromcoeffs(coeffs, precoeff)
    coeffs = real.(coeffs)
    lc = length(coeffs)
    f(x) = sum([precoeff * coeffs[n] * cos(n * x) for n in 1:1])
    return f
end

function plotdisprels!(plt, subplot, G, λ, B)
    ovdict = overlaps[G][λ][11][B]
    for ℓ in [0,1,2]
        if ℓ == 0
            linestyle = :dot
            linewidth = 1.2
        elseif ℓ == 1
            linestyle = :dashdot
            linewidth = 1.4
        elseif ℓ == 2
            linestyle = :solid
            linewidth = 1.6
        end
        Plots.plot!(plt, disprelfromcoeffs(ovdict[ℓ]["tcreator"], 2), subplot = subplot, linestyle = linestyle, linewidth = linewidth, color = :blue)
        Plots.plot!(plt, disprelfromcoeffs(ovdict[ℓ]["tcreatorU"], 2), subplot = subplot, linestyle = linestyle, linewidth = linewidth, color = :green)
    end
    Plots.plot!(plt, disprelfromcoeffs(ovdict[ℓ]["tperturbative"], 1), subplot = subplot, linestyle = :dot, linewidth = 1, color = :red)
end

function plotspectra!(plt, subplot, G, λ, B)
    Nlimt = 80
    thedict = cones[G][λ][l][B][2]
    Mat = thedict["correlationfunction"][1:Nlimt, :]
    Nbigsize = 1000
    Mat = zero_pad_center(Mat, (Nbigsize, Nbigsize))
    Mat = abs.(fftshift(fft(Mat)))
    Acent = map(eachcol(Mat)) do r
        argmax(r)
    end
    centroid1 = sum(Acent) / length(Acent)
    centroid1 = centroid1 - Nbigsize / 2

    centroid = smallsize[G][λ][l]["centroid$B"]
    Elist = [energy(el) - centroid for el ∈ smallsize[G][λ][l]["band$B"]]
    dlt = (maximum(Elist) - minimum(Elist)) / 2
    ctst = (maximum(Elist) + minimum(Elist)) / 2
    ymin = ctst - 3 * dlt
    ymax = ctst + 3 * dlt

    dt = thedict["timestep"]
    shift = -2 * π * centroid1 / (Nbigsize * dt)
    ny, nx = size(Mat)
    x = range(-π, π; length = nx)

    yrangemin = -π/dt + shift
    yrangemax = π/dt + shift
    yrangew = yrangemax - yrangemin

    y = range(yrangemin, yrangemax; length = ny)
    Plots.heatmap!(plt, x, y, Mat, xlims = [-π,π], ylims = [ymin, ymax], c = cmaptransparent(bandcolor(G, B)), legend = false, subplot = subplot)
    y = range(yrangemin + yrangew, yrangemax + yrangew; length = ny)
    Plots.heatmap!(plt, x, y, Mat, xlims = [-π,π], ylims = [ymin, ymax], c = cmaptransparent(bandcolor(G, B)), legend = false, subplot = subplot)
    y = range(yrangemin - yrangew, yrangemax - yrangew; length = ny)
    Plots.heatmap!(plt, x, y, Mat, xlims = [-π,π], ylims = [ymin, ymax], c = cmaptransparent(bandcolor(G, B)), legend = false, subplot = subplot)
end

function plotdisprelsmall!(plt, subplot, G, λ, B)
    energ = dispersion_interpolant(smallsize[G][λ][l]["band$B"])
    centroid = smallsize[G][λ][l]["centroid$B"]
    Elist = [energy(el) - centroid for el ∈ smallsize[G][λ][l]["band$B"]]
    dlt = (maximum(Elist) - minimum(Elist)) / 2
    ctst = (maximum(Elist) + minimum(Elist)) / 2
    # ymin = ctst - 3 * dlt
    # ymax = ctst + 3 * dlt
    energcentered(x) = energ(x) - centroid
    Klist = [momentum(el) for el ∈ smallsize[G][λ][l]["band$B"]]
    # Plots.plot!(plt, energcentered, linewidth = 1, ylims = [ymin, ymax], linestyle = :dash, color = :red, subplot = subplot)
    Plots.scatter!(plt, Klist, Elist, color = :dodgerblue4, subplot = subplot)
end

let

plt = Plots.plot(layout=(4,3),
    legend=false,
    legendfontsize = 4,
    plot_padding = -4mm,
    margin = -4mm,
    framestyle = :box,
    guidefont = font(10),
    tickfont = font(10),
    legendfont = font(10),
    titlefont = font(10),
    foreground_color_legend = nothing,
    background_color_legend = :transparent,
    xlims = [-π,π],
    xticks = ([-π/2, 0, π/2], [L"-\frac{\pi}{2}", L"0", L"\frac{\pi}{2}"])
)

indx = 0
for G in ["ZZ3", "SU3"]
    for B in ["1", "2"]
        for λ in [1//10, 5//10, 9//10]
            indx += 1
            plotspectra!(plt, indx, G, λ, B)
            plotdisprels!(plt, indx, G, λ, B)
            plotdisprelsmall!(plt, indx, G, λ, B)
        end
    end
end

Plots.plot!(plt, size=(570,600), grid=true)
Plots.savefig("plots/spectra.pdf")

end