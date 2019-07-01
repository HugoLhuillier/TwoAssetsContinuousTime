mutable struct Param
    ρ       ::Float64
    σ       ::Float64
    Δ       ::Float64
    χ0      ::Float64
    χ1      ::Float64
    ξ       ::Float64
    ε       ::Float64

    rb      ::Float64
    κ       ::Float64
    τ       ::Float64
    ra      ::Float64
    w       ::Float64

    gA      ::Array{Float64}
    gB      ::Array{Float64}
    gZ      ::Array{Float64}
    da      ::Float64
    db      ::Float64
    dz      ::Float64
    nI      ::Int
    nJ      ::Int
    nK      ::Int
    λ       ::Array{Float64}


    function Param(doAR1::Bool; par = Dict())
        f = open(joinpath(dirname(@__FILE__),"params.json"))
        # f = open("./params.json")
        j = JSON.parse(f)
        close(f)

        this = new()
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

        return this
    end
end

function Base.show(io::IO, p::Param)
    println("parameter values")
    println("-----------------------")
    for (k,v) in OrderedDict(:ρ => "discount factor",
                  :σ     => "risk aversion",
                  :χ0    => "cost function, χ0",
                  :χ1    => "cost function, χ1",
                  :ξ     => "% invested in illiquid asset",
                  :ε     => "tolerance")
        println("$v: $(getfield(p,k))")
    end
    println("-----------------------")
    println("prices")
    println("-----------------------")
    for (k,v) in OrderedDict(:rb => "liquid interest rate",
                  :κ      => "liquid borrowing rate (net)",
                  :ra     => "illiquid interest rate",
                  :w      => "wage rate",
                  :τ      => "tax on liquid rate")
        println("$v: $(getfield(p,k))")
    end
    println("-----------------------")
    println("grids")
    println("-----------------------")
    println("illquid asset: [$(minimum(p.gB)), $(maximum(p.gB))] (# $(p.nI))")
    println("liquid asset: [$(minimum(p.gA)), $(maximum(p.gA))] (# $(p.nJ))")
    println("productivity levels: $(p.gZ) (# $(p.nK))")
    println("-----------------------")
    return nothing
end
