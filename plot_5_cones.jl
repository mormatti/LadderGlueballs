include("usings.jl")
include("plot_simulation.jl")
gr()

L = 31
cones = (JLD2.load("data/conesL$L.jld2"))["single_stored_object"]
ℓ = 2

let

plt = Plots.plot(layout = (4, 3),
        size = (300,450),
        legend = false,
        legendfontsize = 4,
        framestyle = :box,
        guidefont = font(10),
        tickfont = font(10),
        legendfont = font(10),
        titlefont = font(10),
        foreground_color_legend = nothing,
        background_color_legend = :transparent
)

graphind = 0
for G in ["ZZ3", "SU3"]
    for B in ["1", "2"]
        for λ in [9//10, 5//10, 9//10]

            graphind += 1
            matrixenergies = reduce(hcat, cones[G][λ][11][B][ℓ]["evolutiondata"][:energies])'

            Nt = size(matrixenergies, 1)
            Nx = size(matrixenergies, 2)

            # We compute the perturbative velocity

            margin = -4mm
            left_margin = 0mm
            bottom_margin = 0mm
            if λ == 9//10
                ylabel = L"\mathrm{time} \, \, t"
                yticks = ([0, 1/2], [L"0", L"T/2"])
                left_margin = 2mm
            else
                ylabel = ""
                yticks = ([0, 1/2], ["", ""])
            end
            if G == "SU3" && B == "2"
                xlabel = L"\mathrm{position} \, \, x"
                xticks = ([0, L/2], [L"0", L"L/2"])
                bottom_margin = 2mm
            else
                xlabel = ""
                xticks = ([0, L/2], ["", ""])
            end


            plotsimulation!(plt, 
                    "energy", 
                    matrixenergies, 
                    L, 
                    1, 
                    subplot = graphind, 
                    ylabel = ylabel, 
                    xlabel = xlabel,
                    xticks = xticks,
                    yticks = yticks,
                    margin = margin,
                    left_margin = left_margin,
                    bottom_margin = bottom_margin
            )
        end
    end
end

Plots.plot!(plt)

end