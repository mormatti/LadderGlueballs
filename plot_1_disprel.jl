include("usings.jl")

let

# Parameters to set
l = 11
groups = ["ZZ3", "SU3"]
λ_list = [1//10, 3//10, 5//10, 7//10, 9//10]
ncols = length(λ_list)

# Reading data
constants = (JLD2.load("data/constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("data/smallsize.jld2"))["single_stored_object"]

# Group cycle
for group in groups

leb = L"\langle\mathcal{E}\rangle"

texgroup = (group == "ZZ3") ? L"\mathbb{Z}_3" : L"SU(3)"

# Plots setup
plt = Plots.plot(layout=(1,ncols),
    legend = false,
    plot_padding = 0mm,
    margin = -7.5mm,
    framestyle = :box,
    guidefont = font(10),
    tickfont = font(10),
    legendfont = font(10),
    titlefont = font(10),
    xlims  = [0-0.3, π+0.3],
    xticks = ([0, π/2, π], [L"0", L"\frac{\pi}{2}", L"\pi"])
    # title = L"G = " * texgroup * L"\\, \\, \ell = $(l)"
    )

# λ cycle
for ii in 1:ncols
    λ = λ_list[ii]

    titlelatex = L"\lambda = \frac{%$(numerator(λ))}{%$(denominator(λ))}"

    Egs = smallsize[group][λ][l]["gsenergy"]
    Ebar = smallsize[group][λ][l]["centroid1"]
    ΔE = Ebar - Egs

    Elist = smallsize[group][λ][l]["disprel"]["E"]
    Klist = smallsize[group][λ][l]["disprel"]["k"]
    Clist = smallsize[group][λ][l]["disprel"]["C"]
    Plist = smallsize[group][λ][l]["disprel"]["P"]

    # Properties of points
    list_color = []
    list_markerstrokecolor = []
    list_marker = []
    list_markersize = []
    list_markerstrokewidth = []

    nstates = length(Elist)

    for indx in 1:nstates

        # Choose the color
        if Clist[indx] > 0 && group == "SU3"
            push!(list_color, :dodgerblue1)
            push!(list_markerstrokecolor, :dodgerblue4)
        elseif Clist[indx] < 0 && group == "SU3"
            push!(list_color, :purple1)
            push!(list_markerstrokecolor, :purple4)
        elseif Clist[indx] > 0 && group == "ZZ3"
            push!(list_color, :firebrick1)
            push!(list_markerstrokecolor, :firebrick4)
        elseif Clist[indx] < 0 && group == "ZZ3"
            push!(list_color, :darkorange1)
            push!(list_markerstrokecolor, :darkorange4)
        end

        if Plist[indx] == 0
            push!(list_marker, :circle)
            push!(list_markersize, 1.6)
            push!(list_markerstrokewidth, 0.5)
        elseif Plist[indx] > 0
            push!(list_marker, :square)
            push!(list_markersize, 1.4)
            push!(list_markerstrokewidth, 0.5)
        elseif Plist[indx] < 0
            push!(list_marker, :diamond)
            push!(list_markersize, 1.9)
            push!(list_markerstrokewidth, 0.5)
        end
    end

    ytickstoput = ([Egs, Egs + 0.5 * ΔE, Egs + ΔE, Egs + 1.5 * ΔE, Egs + 2ΔE, Egs + 2.5 * ΔE], [L"0", L"\frac{1}{2}" * leb, leb, L"\frac{3}{2}" * leb, L"2" * leb, L"\frac{5}{2}" * leb ])
    leftmargin = 0mm
    if ii != 1
        ytickstoput = ([Egs, Egs + 0.5 * ΔE, Egs + ΔE, Egs + 1.5 * ΔE, Egs + 2ΔE, Egs + 2.5 * ΔE], ["", "", "", "", "", ""])
        leftmargin = -5mm
    end

    titlelatex = L"\lambda = \frac{%$(numerator(λ))}{%$(denominator(λ))}"

    if group == "ZZ3"
        color1 = :firebrick
        color2 = :darkorange3
    else
        color1 = :dodgerblue3
        color2 = :purple3
    end

    Plots.plot!(plt, dispersion_interpolant(smallsize[group][λ][l]["band1"]), subplot = ii, color = color1)
    Plots.plot!(plt, dispersion_interpolant(smallsize[group][λ][l]["band2"]), subplot = ii, color = color2)

    Plots.scatter!(
        plt, 
        Klist, 
        Elist, 
        subplot = ii,
        title = titlelatex,
        color = list_color,
        marker = list_marker,
        markersize = list_markersize,
        markerstrokecolor = list_markerstrokecolor,
        markerstrokewidth = list_markerstrokewidth,
        ylims  = [Egs - 0.05 * (3ΔE), Egs + 2.5 * ΔE + 0.05 * (3ΔE)],
        yticks = ytickstoput,
        xlabel = L"k"
        )
end

Plots.plot!(plt, size=(270,230), grid=true)
Plots.savefig("plots/disprel$(group).pdf")

end # end group cycle

end # let