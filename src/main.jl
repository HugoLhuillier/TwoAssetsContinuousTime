module TwoAssetsContinuousTime

    using JSON
    using SparseArrays, LinearAlgebra, IterativeSolvers, SharedArrays, Arpack
    include(joinpath(dirname(@__FILE__),"Markov/src/Markov.jl"))
    include("param.jl")
    include("household.jl")
    include("diff.jl")
    include("finance.jl")
    include("foc.jl")
    include("mat_const.jl")

    function hjb!(p::Param, hh::Household, maxIter)
        for i in 1:maxIter
            for k in eachindex(p.gZ)
                ∂V!(p,hh,k)
                consumption!(p,hh,k)
                deposit!(p,hh,k)
                B!(p,hh,k)
                D!(p,hh,k)
            end
            A!(p,hh)
            # transition matrix should sum up to zero
            s       = maximum(abs.(sum(hh.A, dims = 2)))
            if s > p.ε; @warn "improper transition matrix, ∑ = $s"; end
            solve_hjb!(p,hh)
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
        for k in eachindex(p.gZ)
            B2!(p,hh,k)
            D2!(p,hh,k)
        end
        TwoAssetsContinuousTime.A!(p,hh)
        # NOTE: compared to ben's code, take the eigenvector associated to
        # the largest eigenvalue (= zero)
        val, vec = Arpack.eigs(hh.A', which = :LR, nev = 1)
        vec     = Real.(vec)
        vec     ./= sum(vec)
        hh.μ    = reshape(vec, (p.nI, p.nJ, p.nK))
        return nothing
    end

    function solution(; doParallel::Bool = false, doAR1::Bool = false, maxIter::Integer = 30)
        p   = Param(doAR1)
        hh  = Household(p, doParallel)
        hjb!(p,hh,maxIter)
        kde!(p, hh)
        return p,hh
    end

end
