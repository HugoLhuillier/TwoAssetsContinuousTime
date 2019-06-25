"computes the state space and the transition matrix of AR(1) based on the
rowenhorst algorithm"
function rowenhorst_estimation!(t)
    # NOTE: treat the process as if had mean zero and add at the end the mean
    # this way, the moment equations are simpler to solve
    # ψ defines the state space. (p,q) defines the transition matrix# 
    p = (1 + t.ρ) / 2
    q = p
    N = length(t.grid)
    ψ = sqrt(t.σ2 / (1 - t.ρ^2)) * sqrt(N - 1)
    μ = t.ω / (1 - t.ρ)
    Θ = Dict(2 => [p 1 - p; 1 - q q])

    # defines the transition matrix as sum of binary markov processes, with
    # transition matrix Θ
    if N > 2
        for n in 3:N
            z    = zeros(n-1)
            Θ[n] = p       * [Θ[n-1] z; z' 0] +
                   (1 - p) * [z Θ[n-1]; 0 z'] +
                   (1 - q) * [z' 0; Θ[n-1] z] +
                   q       * [0 z'; z Θ[n-1]]
            Θ[n][2:(end-1),:] ./= 2
        end
    end
    t.π[:]    = Θ[N]
    t.grid[:] = range(μ - ψ, stop = μ + ψ, length = N)
    return nothing
end

"rowenhorst wrapper"
function rowenhorst(ω::Float64, ρ::Float64, σ2::Float64, N::Int)
    r = AR1(ω, ρ, σ2, N)
    rowenhorst_estimation!(r)
    return r
end
