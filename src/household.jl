mutable struct Household{T <: Union{SharedArray{Float64}, Array{Float64}}, U <: SparseMatrixCSC{Float64, Int64}}

    V       ::T
    Vupdt   ::T
    Vupdtvec::T
    c       ::T
    d       ::T
    adot    ::T
    bdot    ::T
    μ       ::T

    VaB     ::T
    VaF     ::T
    VbB     ::T
    VbF     ::T
    cB      ::T
    cF      ::T
    dBB     ::T
    dBF     ::T
    dFB     ::T
    dFF     ::T
    dB      ::T
    dF      ::T
    scB     ::T
    scF     ::T
    sdB     ::T
    sdF     ::T

    X       ::T
    Y       ::T
    Z       ::T
    Λ       ::U
    B       ::U
    D       ::U
    A       ::U

    function Household(p::Param, doParallel::Bool)
        # state space:
        # first dimension is liquid asset
        # second dimension is illiquid asset
        # third dimension is productivity
        ms       = spzeros(p.nI * p.nJ * p.nK, p.nI * p.nJ * p.nK)
        if doParallel
            this = new{SharedArray{Float64}, SparseMatrixCSC{Float64, Int64}}()
            m    = SharedArray{Float64,3}((p.nI, p.nJ, p.nK))
        else
            this = new{Array{Float64}, SparseMatrixCSC{Float64, Int64}}()
            m    = zeros(p.nI, p.nJ, p.nK)
        end
        for e in fieldnames(Household)
            if !∈(e, [:Λ; :B; :D; :A; :Vupdtvec])
                setfield!(this, e, copy(m))
            elseif e == :Vupdtvec
                this.Vupdtvec = zeros(p.nI * p.nJ * p.nK)
            else
                setfield!(this, e, copy(ms))
            end
        end
        # TODO: so far valid only for 2 state income process
        this.Λ = spdiagm(0           => [-p.λ[1,2] * ones(p.nI * p.nJ); -p.λ[2,1] * ones(p.nI * p.nJ)],
                            -p.nI*p.nJ  => p.λ[2,1] * ones(p.nI * p.nJ),
                            p.nI*p.nJ   => p.λ[1,2] * ones(p.nI * p.nJ))

        for (k,z) in enumerate(p.gZ)
            this.V[:,:,k] = u(p, (1 - p.ξ) * z * p.w .+ (p.ra * p.gA)' .+ (p.rb + p.κ).* p.gB) ./ p.ρ
        end

        return this
    end
end

u(p::Param, c)     = c.^(1 - p.σ) ./ (1 - p.σ)
∂u(p::Param, c)    = c.^(-p.σ)
inv_∂u(p::Param, x) = (x).^(-1 / p.σ)
