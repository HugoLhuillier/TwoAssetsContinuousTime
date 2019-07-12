function B_hjb!(p::Param, hh::Household, ind::Tuple{Int64,Int64,Int64})
    I                   = p.lindex[ind...]
    scF, sdF, scB, sdB  = hh.scF[ind...], hh.sdF[ind...], hh.scB[ind...], hh.sdB[ind...]
    if ind[1] == 1
        hh.B[I,I]   = min(scB, 0.0) / p.db - (max(scF, 0.0) + max(sdF, 0.0)) / p.db
        hh.B[I,I+1] = (max(scF, 0.0) + max(sdF, 0.0)) / p.db
    elseif ind[1] == p.nI
        hh.B[I,I-1] = - (min(0, scB) + sdB) / p.db
        hh.B[I,I]   = (min(scB, 0.0) - max(scF, 0.0) + sdB) / p.db
    else
        hh.B[I,I]   = (min(scB, 0.0) + min(sdB, 0.0)) / p.db - (max(scF, 0.0) + max(sdF, 0.0)) / p.db
        hh.B[I,I+1] = (max(scF, 0.0) + max(sdF, 0.0)) / p.db
        hh.B[I,I-1] = - (min(0, scB) + min(0, sdB)) / p.db
    end
    return nothing
end

# function B!(p::Param, hh::Household, ind::Tuple{Int64,Int64,Int64})
#     # make sure that we use backward difference at bmax and forward difference at bmin for deposit
#     if ind[1] == 1
#         hh.X[ind...]   = - min(0, hh.scB[ind...]) / p.db
#         hh.Y[ind...]   = min(hh.scB[ind...], 0.0) / p.db -
#             (max(hh.scF[ind...], 0.0) + max(hh.sdF[ind...], 0.0)) / p.db
#         hh.Z[ind...]   = (max(hh.scF[ind...], 0.0) + max(hh.sdF[ind...], 0.0)) / p.db
#     elseif ind[1] == p.nI
#         hh.X[ind...]   = - (min(0, hh.scB[ind...]) + hh.sdB[ind...]) / p.db
#         hh.Y[ind...]   = (min(hh.scB[ind...], 0.0) - max(hh.scF[ind...], 0.0) + hh.sdB[ind...]) / p.db
#         hh.Z[ind...]   = max(hh.scF[ind...], 0.0) / p.db
#     else
#         hh.X[ind...]   = - (min(0, hh.scB[ind...]) + min(0, hh.sdB[ind...])) / p.db
#         hh.Y[ind...]   = (min(hh.scB[ind...], 0.0) + min(hh.sdB[ind...], 0.0)) / p.db -
#             (max(hh.scF[ind...], 0.0) + max(hh.sdF[ind...], 0.0)) / p.db
#         hh.Z[ind...]   = (max(hh.scF[ind...], 0.0) + max(hh.sdF[ind...], 0.0)) / p.db
#     end
#     # put the diagonals together in a tridiagonal matrix
#     # TODO: thing of a more clever way to compute this matrix
#     # concB!(p,hh,k)
#     return nothing
# end

function B_kde!(p::Param, hh::Household, ind::Tuple{Int64,Int64,Int64})
    I           = p.lindex[ind...]
    bdot        = hh.bdot[ind...]
    hh.B[I,I]   = (min(bdot, 0) - max(bdot, 0)) / p.db
    if ind[1] == 1
        hh.B[I,I+1] = max(bdot, 0) / p.db
    elseif ind[1] == p.nI
        hh.B[I,I-1] = - min(bdot, 0) / p.db
    else
        hh.B[I,I-1] = - min(bdot, 0) / p.db
        hh.B[I,I+1] = max(bdot, 0) / p.db
    end
    return nothing
end

# function concB!(p::Param, hh::Household, k::Integer)
#     # do not need X[1,:,:] and Z[end,:,:]
#     hh.X[:,:,k]     = [hh.X[2:end,:,k]; zeros(p.nJ)']
#     hh.Z[:,:,k]     = [zeros(p.nJ)'; hh.Z[1:(end-1),:,k]]
#     hh.B[(p.nI*p.nJ*(k-1)+1):(p.nI*p.nJ*k),(p.nI*p.nJ*(k-1)+1):(p.nI*p.nJ*k)] =
#                         spdiagm(-1  => vcat(hh.X[:,:,k]...)[1:(end-1)],
#                                 0   => vcat(hh.Y[:,:,k]...),
#                                 1   => vcat(hh.Z[:,:,k]...)[2:end])
#     return nothing
# end

function D_hjb!(p::Param, hh::Household, ind::Tuple{Int64,Int64,Int64})
    # impose negative drift / use of backward finite difference at the upper bound
    # TODO: understand why...
    I  = p.lindex[ind...]
    d  = hh.d[ind...]
    if ind[2] < p.nJ
        hh.D[I,I]       = (min(d, 0) - max.(d, 0) - p.gArA[ind[2]] - p.gIncA[ind[3]]) / p.da
        hh.D[I,I+p.nI]  = (max(d, 0) + p.gArA[ind[2]] + p.gIncA[ind[3]]) / p.da
        if ind[2] > 1
            hh.D[I,I-p.nI]  =- min(d, 0) / p.da
        end
    else
        hh.D[I,I-p.nI]  = -(min(d,0) + p.gArA[end] + p.gIncA[ind[3]]) / p.da
        hh.D[I,I]       = (min(d,0) + p.gArA[ind[2]] + p.gIncA[ind[3]]) / p.da
    end
    # if ind[2] == p.nJ
    #     hh.Z[ind...] = 0.0
    #     hh.Y[ind...] = (min(hh.d[ind...],0) + p.gArA[ind[2]] + p.gIncA[ind[3]]) / p.da
    #     hh.X[ind...] = -(min(hh.d[ind...],0) + p.gArA[end] + p.gIncA[ind[3]]) / p.da
    # else
    #     hh.X[ind...] = - min(hh.d[ind...], 0) / p.da
    #     hh.Y[ind...] = (min(hh.d[ind...], 0) - max.(hh.d[ind...], 0) -
    #                     p.gArA[ind[2]] - p.gIncA[ind[3]]) / p.da
    #     hh.Z[ind...] = (max(hh.d[ind...], 0) + p.gArA[ind[2]] + pp.gIncA[ind[3]]) / p.da
    # end
    # cast together in sparse matrix
    # concD!(p,hh,k)
    return nothing
end

function D_kde!(p::Param, hh::Household, ind::Tuple{Int64,Int64,Int64})
    I           = p.lindex[ind...]
    adot        = hh.adot[ind...]
    hh.D[I,I]   = (min(adot, 0) - max(adot, 0)) / p.da
    if ind[2] == 1
        hh.D[I,I+p.nI] = max(adot, 0) / p.da
    elseif ind[2] == p.nJ
        hh.D[I,I-p.nI] = - min(adot, 0) / p.da
    else
        hh.D[I,I-p.nI] = - min(adot, 0) / p.da
        hh.D[I,I+p.nI] = max(adot, 0) / p.da
    end
    return nothing
end

# function D_kde!(p::Param, hh::Household, k::Integer)
#     upwind_mat!(hh, hh.adot[:,:,k], p.da, k)
#     concD!(p, hh, k)
#     return nothing
# end

# function concD!(p::Param, hh::Household, k::Integer)
#     hh.D[(p.nI*p.nJ*(k-1)+1):(p.nI*p.nJ*k),(p.nI*p.nJ*(k-1)+1):(p.nI*p.nJ*k)] =
#                         spdiagm(p.nI    => vcat(hh.Z[:,1:(end-1),k]...),
#                                 0       => vcat(hh.Y[:,:,k]...),
#                                 -p.nI   => vcat(hh.X[:,2:end,k]...))
#     return nothing
# end

function A!(p::Param, hh::Household)
    hh.A = hh.B + hh.D + hh.Λ
    return nothing
end

function inverse_hjb!(p::Param, hh::Household)
    hh.Vupdtvec[:]  = vcat(hh.V...)
    Pardiso.solve!(p.ps, hh.Vupdtvec, p.Δmat - hh.A, vcat(u(p, hh.c) + hh.V./ p.Δ...))
    hh.Vupdt[:]     = reshape(hh.Vupdtvec, p.nI, p.nJ, p.nK)
    return nothing
end
