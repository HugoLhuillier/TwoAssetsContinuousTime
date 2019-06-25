function upwind!(O, cB::T, cF::T, dB::T, dF::T, k::Integer) where T <: Array{Float64}
    O[:,:,k] = cB .* (dB .< -1e-12) .+ cF .* (dF .> 1e-12)
    return nothing
end
function upwind!(O, cB::T, cF::T, k::Integer) where T <: Array{Float64}
    upwind!(O, cB, cF, cB, cF, k)
end
function upwind_mat!(hh::Household, dot::Array{Float64}, d::Float64, k::Integer)
    hh.X[:,:,k] = -min.(dot, 0) ./ d
    hh.Y[:,:,k] = (min.(dot, 0) .- max.(dot, 0)) ./ d
    hh.Z[:,:,k] = max.(dot, 0) ./ d
    return nothing
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
