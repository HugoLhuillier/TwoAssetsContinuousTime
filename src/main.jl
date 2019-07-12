module TwoAssetsContinuousTime

    using JSON
    using SparseArrays, LinearAlgebra, SharedArrays
    using Pardiso
    using DataStructures
    include(joinpath(dirname(@__FILE__),"Markov/src/Markov.jl"))
    include("param.jl")
    include("household.jl")
    include("diff.jl")
    include("finance.jl")
    include("foc.jl")
    include("mat_const.jl")
    include("solver.jl")

    function solution(; doParallel::Bool = false, doAR1::Bool = false, maxIter::Integer = 30)
        p   = TwoAssetsContinuousTime.Param(doAR1)
        hh  = TwoAssetsContinuousTime.Household(p, doParallel)
        hjb!(p,hh,maxIter)
        kde!(p, hh)
        return p,hh
    end

end
