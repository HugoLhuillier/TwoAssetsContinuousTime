function B!(p::Param, hh::Household, k::Integer)
    hh.X[:,:,k]     = - (min.(0, hh.scB[:,:,k]) .+ min.(0, hh.sdB[:,:,k])) ./ p.db
    hh.Y[:,:,k]     = (min.(hh.scB[:,:,k], 0.0) .+ min.(hh.sdB[:,:,k], 0.0)) ./ p.db .-
        (max.(hh.scF[:,:,k], 0.0) .+ max.(hh.sdF[:,:,k], 0.0)) ./ p.db
    hh.Z[:,:,k]     = (max.(hh.scF[:,:,k], 0.0) .+ max.(hh.sdF[:,:,k], 0.0)) ./ p.db
    # make sure that we use backward difference at bmax and forward difference at bmin for deposit
    hh.X[1,:,k]     = - min.(0, hh.scB[1,:,k]) ./ p.db
    hh.Y[1,:,k]     = min.(hh.scB[1,:,k], 0.0) ./ p.db .-
        (max.(hh.scF[1,:,k], 0.0) .+ max.(hh.sdF[1,:,k], 0.0)) ./ p.db
    hh.X[end,:,k]   = - (min.(0, hh.scB) .+ hh.sdB)[end,:,k] ./ p.db
    hh.Y[end,:,k]   = (min.(hh.scB, 0.0) .- max.(hh.scF, 0.0) .+ hh.sdB)[end,:,k] ./ p.db
    hh.Z[end,:,k]   = max.(hh.scF[end,:,k], 0.0)  ./ p.db
    # put the diagonals together in a tridiagonal matrix
    concB!(p,hh,k)
    return nothing
end

function B2!(p::Param, hh::Household, k::Integer)
    upwind_mat!(hh, hh.bdot[:,:,k], p.db, k)
    concB!(p,hh,k)
end

function concB!(p::Param, hh::Household, k::Integer)
    # do not need X[1,:,:] and Z[end,:,:]
    hh.X[:,:,k]     = [hh.X[2:end,:,k]; zeros(p.nJ)']
    hh.Z[:,:,k]     = [zeros(p.nJ)'; hh.Z[1:(end-1),:,k]]
    hh.B[(p.nI*p.nJ*(k-1)+1):(p.nI*p.nJ*k),(p.nI*p.nJ*(k-1)+1):(p.nI*p.nJ*k)] =
                        spdiagm(-1  => vcat(hh.X[:,:,k]...)[1:(end-1)],
                                0   => vcat(hh.Y[:,:,k]...),
                                1   => vcat(hh.Z[:,:,k]...)[2:end])
    dropzeros!(hh.B)
    return nothing
end

function D!(p::Param, hh::Household, k::Integer)
    hh.X[:,:,k]     = - min.(hh.d[:,:,k], 0) ./ p.da
    hh.Y[:,:,k]     = (min.(hh.d[:,:,k], 0) .- max.(hh.d[:,:,k], 0) .-
                    (ra(p, p.gA) .* p.gA)' .- p.gZ[k] * p.w * p.ξ) ./ p.da
    hh.Z[:,:,k]     = (max.(hh.d[:,:,k], 0) .+ (ra(p, p.gA) .* p.gA)' .+ p.gZ[k] * p.w * p.ξ) ./ p.da
    # impose negative drift / use of backward finite difference at the upper bound
    # TODO: understand why...
    hh.Z[:,end,k]   = zeros(p.nI)
    hh.Y[:,end,k]   = (min.(hh.d[:,end,k],0) .+ ra(p, p.gA[end]) .* p.gA[end] .+ p.gZ[k] * p.w * p.ξ) ./ p.da
    hh.X[:,end,k]   = -(min.(hh.d[:,end,k],0) .+ ra(p, p.gA[end]) .* p.gA[end] .+ p.gZ[k] * p.w * p.ξ) ./ p.da
    # cast together in sparse matrix
    concD!(p,hh,k)
    return nothing
end

function D2!(p::Param, hh::Household, k::Integer)
    upwind_mat!(hh, hh.adot[:,:,k], p.da, k)
    concD!(p, hh, k)
    return nothing
end

function concD!(p::Param, hh::Household, k::Integer)
    hh.D[(p.nI*p.nJ*(k-1)+1):(p.nI*p.nJ*k),(p.nI*p.nJ*(k-1)+1):(p.nI*p.nJ*k)] =
                        spdiagm(p.nI    => vcat(hh.Z[:,1:(end-1),k]...),
                                0       => vcat(hh.Y[:,:,k]...),
                                -p.nI   => vcat(hh.X[:,2:end,k]...))
    dropzeros!(hh.D)
    return nothing
end

function A!(p::Param, hh::Household)
    hh.A = hh.B + hh.D + hh.Λ
    return nothing
end

function solve_hjb!(p::Param, hh::Household)
    hh.Vupdtvec[:]  = vcat(hh.V...)
    Pardiso.solve!(p.ps, hh.Vupdtvec, p.Δmat - hh.A, vcat(u(p, hh.c) .+ hh.V./ p.Δ...))
    hh.Vupdt[:]     = reshape(hh.Vupdtvec, p.nI, p.nJ, p.nK)
    return nothing
end
