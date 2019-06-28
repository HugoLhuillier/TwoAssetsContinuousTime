mutable struct Param{T <: Float64}
    ρ       ::T
    σ       ::T
    Δ       ::T
    χ0      ::T
    χ1      ::T
    ξ       ::T
    ε       ::T

    rb      ::T
    κ       ::T
    τ       ::T
    ra      ::T
    w       ::T

    gA      ::Array{T}
    gB      ::Array{T}
    gZ      ::Array{T}
    da      ::T
    db      ::T
    dz      ::T
    nI      ::Int
    nJ      ::Int
    nK      ::Int
    λ       ::Array{T}

    ps      ::MKLPardisoSolver

    function Param(doAR1::Bool; par = Dict())
        f = open(joinpath(dirname(@__FILE__),"params.json"))
        # f = open("./params.json")
        j = JSON.parse(f)
        close(f)

        this = new{Float64}()
        for (k,v) in j
            if ∈(Symbol(k), fieldnames(Param))
                setfield!(this,Symbol(k),v["value"])
            end
        end

        if length(par) > 0
            for (k,v) in par
                setfield!(this,k,v)
            end
        end

        # NOTE: equidistant grid points
        this.gA = range(0.0, stop = j["amax"]["value"], length = this.nJ)
        this.gB = range(j["bmin"]["value"], stop = j["bmax"]["value"], length = this.nI)
        # discretize the wage process with rowenhorst
        if(doAR1)
            ω       = - j["σ2"]["value"] * (1 - j["ρ_z"]["value"]) / (2 * (1 - j["ρ_z"]["value"]^2))
            prod    = Markov.rowenhorst(ω, j["ρ_z"]["value"], j["σ2"]["value"], j["nK"]["value"])
            this.gZ = exp.(prod.grid)
            this.λ  = log(prod.π)
        else
            if this.nK > 2; @error "nK must be equal to 2. set doAR1 = true to keep nK > 2"; end
            this.gZ = [.8; 1.3]
            this.λ  = [-1/3 1/3; 1/3 -1/3]
        end
        this.da = this.gA[2] - this.gA[1]
        this.db = this.gB[2] - this.gB[1]
        this.dz = this.gZ[2] - this.gZ[1]

        if this.ra >= 1 / this.χ1; error("ra >= 1 / χ1: ∞ accumulation of illiquid wealth"); end

        this.ps = Pardiso.MKLPardisoSolver()

        return this
    end
end
