rb(p::Param, b)        = ifelse.(b .>= 0, p.rb, p.rb + p.κ)
# NOTE: in their code, precisse that need decreasing returns for very large level of illiquid wealth
# if ra >> rb
ra(p::Param, a)        = p.ra .* (1 .- (1.33 .* p.gA[end] ./ a).^(1 - p.τ))
χ(p::Param, d, a)      = p.χ0 .* abs.(d) .+ (p.χ1 / 2) .* (d ./ max.(a,p.ε)).^2 .* a
