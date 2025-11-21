function plotsimulation!(toplot, observable, matrix, systemlength, simulationtime; kwargs...)

    support = systemlength + 1 - size(matrix, 2)

    color = colorant"black"
    if observable == "energy"
        color = colorant"red"
    elseif observable == "entanglement"
        color = colorant"black"
    elseif observable == "particleZZ3_1"
        color = bandcolor("ZZ3", "1")
    elseif observable == "particleZZ3_2"
        color = bandcolor("ZZ3", "2")
    elseif observable == "particleSU3_1"
        color = bandcolor("SU3", "1")
    elseif observable == "particleSU3_2"
        color = bandcolor("SU3", "2")
    elseif observable == "otherparticles"
        color = colorant"blue"
    else
        color = colorant"black"
    end
    
    cmap = cgrad([
        RGBA(red(color), green(color), blue(color), 0.0),   # :dodgerblue, alpha=0
        RGBA(red(color), green(color), blue(color), 1.0)    # :dodgerblue, alpha=1
    ])

    startpoint = (support - 1)/2
    xrange = range(startpoint, systemlength - startpoint, length = size(matrix, 2))
    yrange = range(0, simulationtime, length = size(matrix, 1))
    Plots.heatmap!(toplot,
            xrange,
            yrange,
            matrix,
            xlims = [0, systemlength], 
            ylims = [0, simulationtime],
            c = cmap;
            kwargs...
    )
end