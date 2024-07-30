using Base: @ntuple, @nexprs
@generated function merge_sorted_oN!(dst::Vector{T}, vs::NTuple{N, Vector{T}})::Nothing where {N,T}
    if N == 1
        quote
            empty!(dst)
            append!(dst, first(vs))
            return nothing
        end
    else
        quote
            L0 = 0
            @nexprs $(N) i->begin
                v_i = vs[i]
                L_i = length(v_i)
                t_i = 1
                L0 += L_i
            end
            resize!(dst, L0)
            @inbounds for t_dst ∈ 1:L0
                min_idx = 0
                local min_val::T
                @nexprs $(N) i->begin
                    if t_i ≤ L_i && (min_idx == 0 || v_i[t_i] < min_val)
                        min_val = v_i[t_i]
                        min_idx = i
                    end
                end
                dst[t_dst] = min_val
                @nexprs $(N) i -> if i == min_idx
                    t_i += 1
                end
            end
            return nothing
        end
    end
end
