module Two_Asset

    using ProgressMeter
    using JSON
    using SparseArrays, LinearAlgebra, IterativeSolvers, SharedArrays
    # include(joinpath(dirname(@__FILE__),"Markov/src/AR1.jl"))
    include("param.jl")
    include("household.jl")
    include("diff.jl")
    include("finance.jl")
    include("foc.jl")
    include("mat_const.jl")

    # NOTE: working but slow...

    function hjb!(p,hh,maxIter)
        # prog = ProgressThresh(p.ε, "residual:")
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
            # hh.A[4081,:]
            if s > p.ε; @warn "improper transition matrix, ∑ = $s"; end
            solve!(p,hh)
            e       = maximum(hh.V .- hh.Vupdt)
            hh.V[:] = hh.Vupdt[:]
            if abs(e) <= p.ε
                println("converged in $i steps w/ residuals $e")
                loms!(p,hh)
                break
            elseif i == maxIter
                warn("convergence failed w/ residuals $e")
                break
            else
                # ProgressMeter.update!(prog, e)
                @info "iteration $i w/ residual = $e"
            end
        end
        return nothing
    end

    function hjb(; maxIter::Int = 30)
        p   = Param();
        hh  = Household(p);
        hjb!(p,hh,maxIter)
        return p,hh
    end

end
