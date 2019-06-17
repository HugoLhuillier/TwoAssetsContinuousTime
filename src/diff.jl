function upwind!(O, cB::T, cF::T, dB::T, dF::T, k::Integer) where T <: Array{Float64}
    O[:,:,k] = cB .* (dB .< 0.0) .+ cF .* (dF .> 0.0)
    return nothing
end
function upwind!(O, cB::T, cF::T, k::Integer) where T <: Array{Float64}
    upwind!(O, cB, cF, cB, cF, k)
end

function ∂V!(p::Param, hh::Household, k::Integer)
    # NOTE: boundary condition
    # impose b >= blow and b <= bmax. then VbB(blow) = u'(binding) and id. for VbF(bmax)
    # NOTE bis: no boundary conditions for Va...
    hh.VbF[1:(p.nI - 1),:,k] = diff(hh.V[:,:,k], dims = 1) ./ p.db
    hh.VbF[p.nI,:,k]         = ∂u(p, (1 - p.ξ) * p.w * p.gZ[k] + rb(p, p.gB[end]) * p.gB[end]) .* ones(p.nJ)
    hh.VbB[2:p.nI,:,k]       = hh.VbF[1:(p.nI - 1),:,k]
    hh.VbB[1,:,k]            = ∂u(p, (1 - p.ξ) * p.w * p.gZ[k] + rb(p, p.gB[1]) * p.gB[1]) .* ones(p.nJ)
    hh.VaF[:,1:(p.nJ - 1),k] = diff(hh.V[:,:,k], dims = 2) ./ p.da
    hh.VaB[:,2:p.nJ,k]       = hh.VaF[:,1:(p.nJ - 1),k]
    return nothing
end
