function hjb_optimal!(p::Param, hh::Household)
    @inbounds for ki in eachindex(p.gZ)
        TwoAssetsContinuousTime.∂V!(p,hh,ki)
        @inbounds for bi in eachindex(p.gB), ai in eachindex(p.gA)
            TwoAssetsContinuousTime.∂V_boundary!(p,hh,(bi,ai,ki))
            TwoAssetsContinuousTime.consumption!(p,hh,(bi,ai,ki))
            TwoAssetsContinuousTime.deposit!(p,hh,(bi,ai,ki))
            TwoAssetsContinuousTime.B_hjb!(p,hh,(bi,ai,ki))
            TwoAssetsContinuousTime.D_hjb!(p,hh,(bi,ai,ki))
        end
    end
    check(hh)
    TwoAssetsContinuousTime.A!(p,hh)
end

function hjb!(p::Param, hh::Household, maxIter)
    for i in 1:maxIter
        hjb_optimal!(p, hh)
        # transition matrix should sum up to zero
        s       = maximum(abs.(sum(hh.A, dims = 2)))
        if s > p.ε; @warn "improper transition matrix, ∑ = $s"; end
        TwoAssetsContinuousTime.inverse_hjb!(p,hh)
        e       = maximum(hh.V .- hh.Vupdt)
        hh.V[:] = hh.Vupdt[:]
        if abs(e) <= p.ε
            println("converged in $i steps w/ residuals $e")
            loms!(p,hh)
            break
        elseif i == maxIter
            @warn "convergence failed w/ residuals $e"
            break
        else
            @info "iteration $i w/ residual = $e"
        end
    end
    return nothing
end

function kde!(p::Param, hh::Household)
    @inbounds for ki in eachindex(p.gZ)
        @inbounds for bi in eachindex(p.gB), ai in eachindex(p.gA)
            B_kde!(p,hh,(bi,ai,ki))
            D_kde!(p,hh,(bi,ai,ki))
        end
    end
    A!(p,hh)
    # @time val, vec = Arpack.eigs(hh.A', which = :LR, nev = 1)
    hh.A      = hh.A'
    hh.A[1,:] = [1; zeros(size(hh.A)[1] - 1)]
    Pardiso.solve!(p.ps, hh.μvec, SparseMatrixCSC(hh.A), [.01;hh.μvec[2:end]])
    # vec     = Real.(vec)
    hh.μvec ./= sum(hh.μvec)
    hh.μ    = reshape(hh.μvec, (p.nI, p.nJ, p.nK))
    return nothing
end

function inverse_hjb!(p::Param, hh::Household)
    hh.Vupdtvec[:]  = vcat(hh.V...)
    Pardiso.solve!(p.ps, hh.Vupdtvec, p.Δmat - hh.A, vcat(u(p, hh.c) + hh.V./ p.Δ...))
    hh.Vupdt[:]     = reshape(hh.Vupdtvec, p.nI, p.nJ, p.nK)
    return nothing
end
