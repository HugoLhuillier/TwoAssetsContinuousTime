# FIXME: be very cautious about boundary conditions. so far, we've been very handwaivy about it.

function consumption!(p, hh, k::Int)
    # boundary conditions
    # for the upper bound, this should not be a problem. if the upper bound is large enough,
    # households will decumulate for large a/b and therefore only the backward finite
    # difference is needed.
    # for the lower bound = borrowing constraint, the backward finite difference is
    # defined as if the borrowing constraint were binding, e.g. the drift is set to zero
    hh.cF[:,:,k] = inv_∂u(p, hh.VbF[:,:,k])
    hh.cB[:,:,k] = inv_∂u(p, hh.VbB[:,:,k])
    hh.scB[:,:,k] = (1 - p.ξ) * p.w * p.gZ[k] .+ rb(p,p.gB) .* p.gB .- hh.cB[:,:,k]
    hh.scF[:,:,k] = (1 - p.ξ) * p.w * p.gZ[k] .+ rb(p,p.gB) .* p.gB .- hh.cF[:,:,k]
    hh.c[:,:,k] = upwind(p, hh.cB[:,:,k], hh.cF[:,:,k], hh.scB[:,:,k], hh.scF[:,:,k])
    # TODO: understand better the origin of this... could think of it as "stay put" = zero saving
    hh.c[:,:,k] = ifelse.(hh.scF[:,:,k] .<= 0.0, ifelse.(hh.scB[:,:,k] .>= 0.0,
                        (1 - p.ξ) * p.w * p.gZ[k] .+ rb(p,p.gB) .* p.gB .+ zeros(p.nJ)',
                        hh.c[:,:,k]), hh.c[:,:,k])
    return nothing
end

deposit(p, Va, Vb) = (min.(0., Va ./ Vb .- 1 .+ p.χ0) + max.(0., Va ./ Vb .- 1 .- p.χ0)) .* p.gA' ./ p.χ1

function deposit!(p, hh, k::Int)
    hh.dFF[:,:,k]   = deposit(p, hh.VaF[:,:,k], hh.VbF[:,:,k])
    hh.dBF[:,:,k]   = deposit(p, hh.VaF[:,:,k], hh.VbB[:,:,k])
    hh.dBB[:,:,k]   = deposit(p, hh.VaB[:,:,k], hh.VbB[:,:,k])
    hh.dFB[:,:,k]   = deposit(p, hh.VaB[:,:,k], hh.VbF[:,:,k])
    if any(hh.dBF[:,1,k] .< 0.); @warn "negative deposite rate (bf) at amin for k = $k"; end
    if any(hh.dBB[:,end,k] .> 0.); @warn "positive deposite rate (bb) at amax for k = $k"; end
    if any(hh.dFF[:,1,k] .< 0.); @warn "negative deposite rate (bf) at amin for k = $k"; end
    if any(hh.dFB[:,end,k] .> 0.); @warn "positive deposite rate (bb) at amax for k = $k"; end
    hh.dB[:,:,k]    = upwind(p, hh.dBB[:,:,k], hh.dBF[:,:,k])
    hh.dF[:,:,k]    = upwind(p, hh.dFB[:,:,k], hh.dFF[:,:,k])
    # NOTE: have this in their code, but not sure why
    # seems to work without it
    # hh.dB[:,1,k]    = (hh.dBF[:,1,k] .> p.ε) .* hh.dBF[:,1,k]
    # hh.dB[:,end,k]  = (hh.dBB[:,end,k] .< - p.ε) .* hh.dBB[:,end,k]
    # hh.dB[1,1,k]    = max(hh.dBB[1,1,k], 0.0)
    # hh.dF[:,1,k]    = (hh.dFF[:,1,k] .> p.ε) .* hh.dFF[:,1,k]
    # hh.dF[:,end,k]  = (hh.dFB[:,end,k] .< - p.ε) .* hh.dFB[:,end,k]

    hh.sdB[:,:,k]   = -hh.dB[:,:,k] .-  χ(p, hh.dB[:,:,k], reshape(repeat(p.gA, inner = p.nI), p.nI, p.nJ))
    hh.sdF[:,:,k]   = -hh.dF[:,:,k] .-  χ(p, hh.dF[:,:,k], reshape(repeat(p.gA, inner = p.nI), p.nI, p.nJ))
    # NOTE: have this in their code, but not sure why
    # seems to work without it
    # hh.sdF[end,:,k] = min.(hh.sdF[end,:,k],0.0)
    hh.d[:,:,k]     = upwind(p, hh.dB[:,:,k], hh.dF[:,:,k], hh.sdB[:,:,k], hh.sdF[:,:,k])
    #NOTE: here, if intermediary case, then d = 0
    # ensure that use the backward solution at the upper bound of the grid
    hh.d[end,:,k]   = hh.dB[end,:,k]
    return nothing
end

function loms!(p,hh)
    for (k,z) in enumerate(p.gZ)
        hh.adot[:,:,k] = hh.d[:,:,k] .+ p.ξ * p.w * z .+ (ra(p,p.gA) .* p.gA)'
        hh.bdot[:,:,k] = (1 - p.ξ) * p.w * z .+ rb(p,p.gB) .* p.gB .- hh.d[:,:,k] .-
                        χ(p, hh.d[:,:,k], reshape(repeat(p.gA, inner = p.nI), p.nI, p.nJ)) .-
                        hh.c[:,:,k]
    end
    return nothing
end
