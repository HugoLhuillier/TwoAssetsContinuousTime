"stores parameters of AR(1) and discretized representation"
mutable struct AR1
    ω       ::Float64           # intercept
    ρ       ::Float64           # persistence parameter
    σ2      ::Float64           # variance of (normal) shock
    grid    ::Array{Float64}    # state space
    π       ::Array{Float64}    # transition matrix

    function AR1(ω::Float64, ρ::Float64, σ2::Float64, N::Int)
        this            = new()
        this.ω          = ω
        this.ρ          = ρ
        this.σ2         = σ2
        this.grid       = zeros(N)
        this.π          = zeros(N,N)
        return this
    end
end
