mutable struct Household{T <: Union{SharedArray{Float64}, Array{Float64}}, U <: SparseMatrixCSC{Float64, Int64}}
    V       ::T
    Vupdt   ::T
    Vupdtvec::T
    c       ::T
    d       ::T
    adot    ::T
    bdot    ::T
    μ       ::T
    μvec    ::T

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
        else
            this = new{Array{Float64}, SparseMatrixCSC{Float64, Int64}}()
            m    = zeros(p.nI, p.nJ, p.nK)
        end
        for e in fieldnames(Household)
            if !∈(e, [:Λ; :B; :D; :A; :Vupdtvec; :μvec])
                setfield!(this, e, doParallel ? SharedArray{Float64,3}((p.nI, p.nJ, p.nK)) : copy(m))
            elseif ∈(e, [:Vupdtvec; :μvec])
                setfield!(this, e, zeros(p.nI * p.nJ * p.nK))
            elseif ∈(e, [:Λ; :B; :D; :A])
                setfield!(this, e, copy(ms))
            end
        end
        Λ = p.λ
        for (k,z) in enumerate(p.gZ)
            this.V[:,:,k] = u(p, (1 - p.ξ) * z * p.w .+ (p.ra * p.gA)' .+ (p.rb + p.κ).* p.gB) ./ p.ρ
            Λ[k,k]        = -sum(p.λ[k,:]) + p.λ[k,k]
        end
        this.Λ = kron(Λ, spdiagm(0 => ones(p.nI * p.nJ)))

        return this
    end
end

function Base.show(io::IO, hh::Household)
    nI, nJ, nK = size(hh.V)[1], size(hh.V)[2], size(hh.V)[3]
    println("household type with ($nI, $nJ, $nK) grid points.")
    return nothing
end

u(p::Param, c)     = c.^(1 - p.σ) ./ (1 - p.σ)
∂u(p::Param, c)    = c.^(-p.σ)
inv_∂u(p::Param, x) = (x).^(-1 / p.σ)
