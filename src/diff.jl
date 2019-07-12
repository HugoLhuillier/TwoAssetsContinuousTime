function upwind!(O, cB::T, cF::T, dB::T, dF::T, ind::Tuple{Int64,Int64,Int64}) where T <: Float64
    O[ind...] = cB * (dB < -1e-12) + cF * (dF > 1e-12)
    return nothing
end
function upwind!(O, cB::T, cF::T, ind::Tuple{Int64,Int64,Int64}) where T <: Float64
    upwind!(O, cB, cF, cB, cF, ind)
end
# function upwind_mat!(hh::Household, dot::Array{Float64}, d::Float64, k::Integer)
#     hh.X[:,:,k] = -min.(dot, 0) ./ d
#     hh.Y[:,:,k] = (min.(dot, 0) .- max.(dot, 0)) ./ d
#     hh.Z[:,:,k] = max.(dot, 0) ./ d
#     return nothing
# end

function ∂V!(p::Param, hh::Household, k::Int)
    hh.VbF[1:(p.nI - 1),:,k] = diff(hh.V[:,:,k], dims = 1) ./ p.db
    hh.VbB[2:p.nI,:,k]       = hh.VbF[1:(p.nI - 1),:,k]
    hh.VaF[:,1:(p.nJ - 1),k] = diff(hh.V[:,:,k], dims = 2) ./ p.da
    hh.VaB[:,2:p.nJ,k]       = hh.VaF[:,1:(p.nJ - 1),k]
    return nothing
end

function ∂V_boundary!(p::Param, hh::Household, ind::Tuple{Int64,Int64,Int64})
    # impose b >= blow and b <= bmax. then VbB(blow) = u'(binding) and id. for VbF(bmax)
    if ind[1] == 1
        hh.VbB[ind...]  = ∂u(p, (1 - p.ξ) * p.w * p.gZ[ind[3]] + p.gBrB[ind[1]])
    elseif ind[1] == p.nI
        hh.VbF[ind...]  = ∂u(p, (1 - p.ξ) * p.w * p.gZ[ind[3]] + p.gBrB[ind[1]])
    end
    return nothing
end
