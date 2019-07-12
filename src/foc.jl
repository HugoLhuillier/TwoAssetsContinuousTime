function consumption!(p::Param, hh::Household, ind::Tuple{Int64,Int64,Int64})
    hh.cF[ind...]  = inv_∂u(p, max(p.ε, hh.VbF[ind...]))
    hh.cB[ind...]  = inv_∂u(p, max(p.ε, hh.VbB[ind...]))
    hh.scB[ind...] = p.gInc[ind[3]] + p.gBrB[ind[1]] - hh.cB[ind...]
    hh.scF[ind...] = p.gInc[ind[3]] + p.gBrB[ind[1]] - hh.cF[ind...]
    upwind!(hh.c, hh.cB[ind...], hh.cF[ind...], hh.scB[ind...], hh.scF[ind...], ind)
    # TODO: understand better the origin of this... could think of it as "stay put" = zero saving
    if hh.scB[ind...] >= -p.ε && hh.scF[ind...] <= p.ε
        hh.c[ind...] = p.gInc[ind[3]] + p.gBrB[ind[1]]
    end
    return nothing
end

function deposit!(p::Param, hh::Household, ind::Tuple{Int64,Int64,Int64})
    VaF, VaB, VbF, VbB = hh.VaF[ind...], hh.VaB[ind...], hh.VbF[ind...], hh.VbB[ind...]
    hh.dFF[ind...]   = (min(0.,VaF/VbF - 1 + p.χ0) + max(0.,VaF/VbF - 1 - p.χ0)) * p.gA[ind[2]]/p.χ1
    hh.dBF[ind...]   = (min(0.,VaF/VbB - 1 + p.χ0) + max(0.,VaF/VbB - 1 - p.χ0)) * p.gA[ind[2]]/p.χ1
    hh.dBB[ind...]   = (min(0.,VaB/VbB - 1 + p.χ0) + max(0.,VaB/VbB - 1 - p.χ0)) * p.gA[ind[2]]/p.χ1
    hh.dFB[ind...]   = (min(0.,VaB/VbF - 1 + p.χ0) + max(0.,VaB/VbF - 1 - p.χ0)) * p.gA[ind[2]]/p.χ1
    upwind!(hh.dB, hh.dBB[ind...], hh.dBF[ind...], ind)
    upwind!(hh.dF, hh.dFB[ind...], hh.dFF[ind...], ind)
    # boundary conditions
    if ind[2] == 1
        # NOTE: have this in their code, but not sure why + not necessary for algo to converge
        # seems to work without it
        hh.dB[ind...]    = (hh.dBF[ind...] > p.ε) * hh.dBF[ind...]
        if ind[1] == 1; hh.dB[ind...]  = max(hh.dBB[ind...], 0.0); end
        hh.dF[ind...]    = (hh.dFF[ind...] > p.ε) * hh.dFF[ind...]
    elseif ind[2] == p.nJ
        hh.dB[ind...]    = (hh.dBB[ind...] < - p.ε) * hh.dBB[ind...]
        hh.dF[ind...]    = (hh.dFB[ind...] < - p.ε) * hh.dFB[ind...]
    end
    hh.sdB[ind...]   = - hh.dB[ind...] - χ(p, hh.dB[ind...], p.gA[ind[2]])
    hh.sdF[ind...]   = - hh.dF[ind...] - χ(p, hh.dF[ind...], p.gA[ind[2]])
    # boundary conditions
    if ind[1] == p.nI
        # NOTE: have this in their code, but not sure why
        # seems to work without it
        hh.sdF[ind...] = min(hh.sdF[ind...],0.0)
    end
    upwind!(hh.d, hh.dB[ind...], hh.dF[ind...], hh.sdB[ind...], hh.sdF[ind...], ind)
    # boudary condition
    if ind[1] == p.nI
        # ensure that use the backward solution at the upper bound of the grid
        hh.d[ind...]   = hh.dB[ind...]
    end
    return nothing
end

function check(hh::Household)
    if any(hh.dBF[:,1,:] .< 0.); @warn "negative deposite rate (bf) at amin"; end
    if any(hh.dBB[:,end,:] .> 0.); @warn "positive deposite rate (bb) at amax"; end
    if any(hh.dFF[:,1,:] .< 0.); @warn "negative deposite rate (bf) at amin"; end
    if any(hh.dFB[:,end,:] .> 0.); @warn "positive deposite rate (bb) at amax"; end
end

function loms!(p::Param, hh::Household)
    @inbounds for bi in eachindex(p.gB), ai in eachindex(p.gA), ki in eachindex(p.gZ)
        hh.adot[bi,ai,ki] = hh.d[bi,ai,ki] + p.gIncA[ki] + p.gArA[ai]
        hh.bdot[bi,ai,ki] = p.gInc[ki] + p.gBrB[bi] - hh.d[bi,ai,ki] -
                            χ(p, hh.d[bi,ai,ki], p.gA[ai]) - hh.c[bi,ai,ki]
    end
    return nothing
end
