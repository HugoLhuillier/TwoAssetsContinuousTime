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

function A!(p::Param, hh::Household)
    hh.A = hh.B + hh.D + hh.Λ
    return nothing
end
