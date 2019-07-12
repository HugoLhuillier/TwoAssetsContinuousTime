rb(p::Param, b)        = ifelse.(b .>= 0, p.rb, p.rb + p.κ)
ra(p::Param, a)        = p.ra .* (1 .- (1.33 .* p.gA[end] ./ a).^(1 - p.τ))
χ(p::Param, d::Float64, a::Float64) = p.χ0 * abs(d) + (p.χ1 / 2) * (d / max(a,p.ε))^2 * a
