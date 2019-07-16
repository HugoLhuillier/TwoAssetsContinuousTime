function hjb!(p::Param, hh::Household, maxIter)
    for i in 1:maxIter
        @inbounds for k in eachindex(p.gZ)
            ∂V!(p,hh,k)
            consumption!(p,hh,k)
            deposit!(p,hh,k)
            B!(p,hh,k)
            D!(p,hh,k)
        end
        check(hh)
        A!(p,hh)
        # transition matrix should sum up to zero
        s       = maximum(abs.(sum(hh.A, dims = 2)))
        if s > p.ε; @warn "improper transition matrix, ∑ = $s"; end
        solve_hjb!(p,hh)
        e       = maximum(hh.V .- hh.Vupdt)
        # maybe replace this with copyto!(hh.V, hh.Vupdt)
        copyto!(hh.V, hh.Vupdt)
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
    @inbounds for k in eachindex(p.gZ)
        B2!(p,hh,k)
        D2!(p,hh,k)
    end
    TwoAssetsContinuousTime.A!(p,hh)
    hh.A      = hh.A'
    hh.A[1,:] = [1; zeros(size(hh.A)[1] - 1)]
    Pardiso.solve!(p.ps, hh.μvec, SparseMatrixCSC(hh.A), [.01;hh.μvec[2:end]])
    # vec     = Real.(vec)
    hh.μvec ./= sum(hh.μvec)
    hh.μ    = reshape(hh.μvec, (p.nI, p.nJ, p.nK))
    return nothing
end

function solve_hjb!(p::Param, hh::Household)
    hh.Vupdtvec[:]  = vcat(hh.V...)
    Pardiso.solve!(p.ps, hh.Vupdtvec, p.Δmat - hh.A, vcat(u(p, hh.c) .+ hh.V./ p.Δ...))
    hh.Vupdt[:]     = reshape(hh.Vupdtvec, p.nI, p.nJ, p.nK)
    return nothing
end
