include("usings.jl")

# Local scope
let

# Parameters to set
groups = ["ZZ3", "SU3"]
lambdas = [a//30 for a in 1:29]
lengths = [11]

# Defining the plots
legnd = Plots.plot(size = (100, 100), legend_column = 2, legend_foreground_color = :white, xlims = (0, 1), ylims = (0, 1), xaxis=false, yaxis=false, framestyle=:none, legend=:topleft, legend_font_halign = :left)

# Define common marker settings
commonmarkersettings = (marker = :circle, markersize = 2, markerstrokewidth = 1)

# Define markers for legend
marker0ppzz3 = (color = :firebrick3, label = L"G = \mathbb{Z}_3; \, \, 0^{++}", seriescolor = :firebrick3, markercolor = :firebrick1, markerstrokecolor = :firebrick3)
Plots.plot!(legnd, [-1], [-1]; marker0ppzz3..., commonmarkersettings...)
marker0mmzz3 = (color = :darkorange3, label = L"G = \mathbb{Z}_3; \, \, 0^{--}", seriescolor = :darkorange3, markercolor = :darkorange1, markerstrokecolor = :darkorange3)
Plots.plot!(legnd, [-1], [-1]; marker0mmzz3..., commonmarkersettings...)
marker0ppsu3 = (color = :dodgerblue3, label = L"G = SU(3); \, \, 0^{++}", seriescolor = :dodgerblue3, markercolor = :dodgerblue1, markerstrokecolor = :dodgerblue3)
Plots.plot!(legnd, [-1], [-1]; marker0ppsu3..., commonmarkersettings...)
marker0mmsu3 = (color = :purple3, label = L"G = SU(3); \, \, 0^{--}", seriescolor = :purple3, markercolor = :purple1, markerstrokecolor = :purple3)
Plots.plot!(legnd, [-1], [-1]; marker0mmsu3..., commonmarkersettings...)

Plots.savefig("glueballs/plots/spread_legend.pdf")

pltmin = Plots.plot(size = (300, 200), framestyle = :box, grid=true, legend = false, guidefont = font(10), tickfont = font(10), titlefont = font(10), xlabel = L"\lambda", ylabel = L"\sigma^2(\theta^*)")

# Reading data
constants = (JLD2.load("glueballs/data/0_constants.jld2"))["single_stored_object"]
smallsize = (JLD2.load("glueballs/data/1_smallsize.jld2"))["single_stored_object"]
wannier = (JLD2.load("glueballs/data/2_wannier.jld2"))["single_stored_object"]

mindict = Dict()
mindict["ZZ3"] = Dict()
mindict["SU3"] = Dict()
mindict["ZZ3"]["band1"] = Dict()
mindict["ZZ3"]["band2"] = Dict()
mindict["SU3"]["band1"] = Dict()
mindict["SU3"]["band2"] = Dict()

# Creating the dictionaries of minimums
for group ∈ groups
    for λ ∈ lambdas
        for Li ∈ lengths
            mindict[group]["band1"][λ] = wannier[group][λ][Li]["info1"]["minimum"]
            mindict[group]["band2"][λ] = wannier[group][λ][Li]["info2"]["minimum"]
        end
    end
end

Plots.plot!(pltmin, mindict["ZZ3"]["band1"]; marker0ppzz3..., commonmarkersettings...)
Plots.plot!(pltmin, mindict["ZZ3"]["band2"]; marker0mmzz3..., commonmarkersettings...)
Plots.plot!(pltmin, mindict["SU3"]["band1"]; marker0ppsu3..., commonmarkersettings...)
Plots.plot!(pltmin, mindict["SU3"]["band2"]; marker0mmsu3..., commonmarkersettings...)

Plots.savefig("glueballs/plots/spread.pdf")

end # end local scope