include("usings.jl")

# Local scope
let

# Parameters to set
reset = false
lambdas = [1//10, 5//10, 9//10]
length = 11

lambdalabels = [L"\lambda = \frac{1}{10}", L"\lambda = \frac{1}{2}", L"\lambda = \frac{9}{10}"]

function colortochoose(G, B, a)
    clr = bandcolor(G, B)
    rc = red(clr)
    gc = green(clr)
    bc = blue(clr)
    return RGB(rc / a, gc / a, bc / a)
end

colorsZZ3band1 = [colortochoose("ZZ3", "1", a) for a in [1,1.5,10]]
colorsZZ3band2 = [colortochoose("ZZ3", "2", a) for a in [1,1.5,10]]
colorsSU3band1 = [colortochoose("SU3", "1", a) for a in [1,1.5,10]]
colorsSU3band2 = [colortochoose("SU3", "2", a) for a in [1,1.5,10]]

# Reading data
# constants = (JLD2.load("glueballs/data/0_constants.jld2"))["single_stored_object"]
# smallsize = (JLD2.load("glueballs/data/1_smallsize.jld2"))["single_stored_object"]
wannier = (JLD2.load("data/wannier.jld2"))["single_stored_object"]

plt = Plots.plot(layout=(2,2),
    size=(280,200),
    legend=:topleft,
    legendfontsize = 4,
    plot_padding = 0mm,
    margin = -7.5mm,
    framestyle = :box,
    guidefont = font(10),
    tickfont = font(10),
    legendfont = font(10),
    titlefont = font(10),
    xlims  = [0.5, length + 0.5],
    foreground_color_legend = nothing,
    background_color_legend = :transparent,
    ylims  = [-0.1, 1.1],
    xticks = ([1, (length+1)/2, length], [L"1", L"j_c", L"\ell"])
    )

for λi ∈ reverse(eachindex(lambdas))

λ = lambdas[λi]

energiesZZ3band1 = wannier["ZZ3"][λ][length]["info1"]["density"]
energiesZZ3band2 = wannier["ZZ3"][λ][length]["info2"]["density"]
energiesSU3band1 = wannier["SU3"][λ][length]["info1"]["density"]
energiesSU3band2 = wannier["SU3"][λ][length]["info2"]["density"]

energiesZZ3band1 = Dict(jj => abs(energiesZZ3band1[abs(Int(round(jj - (length+1)/2)))]) for jj in 1:length)
energiesZZ3band2 = Dict(jj => abs(energiesZZ3band2[abs(Int(round(jj - (length+1)/2)))]) for jj in 1:length)
energiesSU3band1 = Dict(jj => abs(energiesSU3band1[abs(Int(round(jj - (length+1)/2)))]) for jj in 1:length)
energiesSU3band2 = Dict(jj => abs(energiesSU3band2[abs(Int(round(jj - (length+1)/2)))]) for jj in 1:length)

Plots.plot!(plt, energiesZZ3band1, legend = false, subplot = 1, label = "", yticks = ([0, 0.25, 0.5, 0.75, 1], [L"0.00", L"0.25", L"0.50", L"0.75", L"1.00"]), ylabel = L"\mathcal{E}_j/{\langle \mathcal{E} \rangle}, \, \, \mathbb{Z}_3", color = colorsZZ3band1[λi], title = L"\mathrm{First \, \, band} \, \, 0^{++}_1")
Plots.plot!(plt, energiesZZ3band2, legend = false, subplot = 2, label = "", yticks = ([0, 0.25, 0.5, 0.75, 1], ["", "", "", "", ""]), color = colorsZZ3band2[λi], title = L"\mathrm{Second \, \, band} \, \, 0^{--}_1")
Plots.plot!(plt, energiesSU3band1, legend = false, subplot = 3, label = "", yticks = ([0, 0.25, 0.5, 0.75, 1], [L"0.00", L"0.25", L"0.50", L"0.75", L"1.00"]), xlabel = L"j", ylabel = L"\mathcal{E}_j/{\langle \mathcal{E} \rangle}, \, \, SU(3)", color = colorsSU3band1[λi])
Plots.plot!(plt, energiesSU3band2, legend = false, subplot = 4, label = "", yticks = ([0, 0.25, 0.5, 0.75, 1], ["", "", "", "", ""]), xlabel = L"j", color = colorsSU3band2[λi])

energiesZZ3band1forlegend = Dict(k => v - 100 for (k, v) in energiesZZ3band1)
energiesZZ3band2forlegend = Dict(k => v - 100 for (k, v) in energiesZZ3band2)
energiesSU3band1forlegend = Dict(k => v - 100 for (k, v) in energiesSU3band1)
energiesSU3band2forlegend = Dict(k => v - 100 for (k, v) in energiesSU3band2)

Plots.scatter!(plt, legend = :topleft, marker = :square, markerstrokewidth = 0, energiesZZ3band1forlegend, subplot = 1, color = colorsZZ3band1[λi], label = lambdalabels[λi])
Plots.scatter!(plt, legend = :topleft, marker = :square, markerstrokewidth = 0, energiesZZ3band2forlegend, subplot = 2, color = colorsZZ3band2[λi], label = lambdalabels[λi])
Plots.scatter!(plt, legend = :topleft, marker = :square, markerstrokewidth = 0, energiesSU3band1forlegend, subplot = 3, color = colorsSU3band1[λi], label = lambdalabels[λi])
Plots.scatter!(plt, legend = :topleft, marker = :square, markerstrokewidth = 0, energiesSU3band2forlegend, subplot = 4, color = colorsSU3band2[λi], label = lambdalabels[λi])

end # end lambda loop

Plots.savefig("plots/wanniers.pdf")

end # end local scope