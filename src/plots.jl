using Plots; pyplot(border = :box)

function reproduce_plots(p::Param, hh::Household)
    for k in eachindex(p.gZ)
        for (o,t) in Dict([:c       => ("Consumption", (240,13)),
                            :d      => ("Deposit", (250,13)),
                            :bdot   => ("Liquid savings",(120,13)),
                            :adot   => ("Illiquid savings",  (250,13)),
                            :μ      => ("Stationary dist.", (352,30))])
            if o == :μ
                surface(p.gB[1:25], p.gA[1:20], getfield(hh, o)[1:25,1:20,k]')
            else
                surface(p.gB, p.gA, getfield(hh, o)[:,:,k]')
            end
            plot!(xlab = "Illiquid wealth", ylab = "Liquid wealth", title = t[1],
                c = :matter, colorbar = false, camera = t[2])
            # plots!(camera = t[2])
            savefig(joinpath(dirname(@__FILE__),"/figs/$(o)_$k"))
        end
    end
    return nothing
end
