function B!(p, hh, k::Int)
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
    # put the diagonals together in a tridiagonal matrix. do not need X[1,:,:] and Z[end,:,:]
    hh.X[:,:,k]     = [hh.X[2:end,:,k]; zeros(p.nJ)']
    hh.Z[:,:,k]     = [zeros(p.nJ)'; hh.Z[1:(end-1),:,k]]
    hh.B[(p.nI*p.nJ*(k-1)+1):(p.nI*p.nJ*k),(p.nI*p.nJ*(k-1)+1):(p.nI*p.nJ*k)] =
                        spdiagm(-1  => vcat(hh.X[:,:,k]...)[1:(end-1)],
                                0   => vcat(hh.Y[:,:,k]...),
                                1   => vcat(hh.Z[:,:,k]...)[2:end])
    return nothing
end

function D!(p, hh, k::Int)
    # NOTE: last column is different from theirs
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
    hh.D[(p.nI*p.nJ*(k-1)+1):(p.nI*p.nJ*k),(p.nI*p.nJ*(k-1)+1):(p.nI*p.nJ*k)] =
                        spdiagm(p.nI    => vcat(hh.Z[:,1:(end-1),k]...),
                                0       => vcat(hh.Y[:,:,k]...),
                                -p.nI   => vcat(hh.X[:,2:end,k]...))
    return nothing
end

function A!(p,hh)
    hh.A[:] = hh.B .+ hh.D .+ hh.Λ
    return nothing
end

function solve!(p,hh)
    hh.Vupdtvec[:]  = vcat(hh.V...)
    idrs!(hh.Vupdtvec, (1 / p.Δ + p.ρ) .* spdiagm(0 => ones(p.nI*p.nJ*p.nK)) - hh.A, vcat(u(p, hh.c) .+ hh.V./ p.Δ...))
    hh.Vupdt[:]     = reshape(hh.Vupdtvec, p.nI, p.nJ, p.nK)
    return nothing
end
