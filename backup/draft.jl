# function get_slice!(ψ::Matrix{T}, x::Wsheet, τ::f64 = 0.) where {T<:Real}
#     if iszero(τ)
#         for i ∈ eachindex(ψ, x.wl)
#             ψ[i] = x[i][end].n_L
#         end
#     else
#         for i ∈ eachindex(ψ, x.wl)
#             ψ[i] = element_around(x[i], τ, +1).n_L
#         end
#     end
#     return nothing
# end
# get_slice(li::Wline)::StateType = li[end].n_L
# get_slice(li::Wline, τ::f64)::StateType = element_around(li, τ, +1).n_L