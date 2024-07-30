
using Base: @ntuple, @nexprs

@generated function merge_sorted!(dst::Vector{T}, vs::Vararg{Vector{T}, N}) where {T,N}
    if N == 1
        quote
            empty!(dst)
            append!(dst, first(vs))
            return nothing
        end
    end
    quote
        empty!(dst)
        L0 = 0
        @nexprs $(N) i->begin
            v_i = vs[i]
            L_i = length(v_i)
            t_i = 1
            L0 += L_i
        end
        min_val = nothing
        min_idx = 0
        r = 0
        updated = false
        while true
            r += 1
            @nexprs $(N) i->begin
                if t_i ≤ L_i && (isnothing(min_val) || v_i[t_i] ≤ min_val)
                    min_val = v_i[t_i]
                    min_idx = i
                    updated = true
                end
            end
            println("round - $(r), select $(min_idx), val = $min_val, updated=$(updated)")
            println(updated)
            if updated
                push!(dst, min_val)
                @nexprs $(N) i->begin
                    if i == min_idx
                        t_i += 1
                    end
                end
            else
                @assert length(dst) == L0
                return nothing
            end
        end
        # merged = T[]
        # indices = [1 for _ in 1:N]
        # lengths = [length(vectors[i]) for i in 1:N]
        
        # while true
        #     min_val, min_idx = nothing, nothing            
        #     for i in 1:N
        #         if indices[i] <= lengths[i]
        #             if min_val === nothing || vectors[i][indices[i]] < min_val
        #                 min_val = vectors[i][indices[i]]
        #                 min_idx = i
        #             end
        #         end
        #     end
            
        #     if min_val === nothing
        #         break
        #     else
        #         push!(merged, min_val)
        #         indices[min_idx] += 1
        #     end
        # end
        
        # return merged
    end
end


# aa = rand(10)
# using BenchmarkTools
# @btime partialsort!($aa,1:2)

v0 = Int[]
v1 = [1, 4, 7]
v2 = [2, 5, 8]
v3 = [3, 6, 9]
merge_sorted!(v0, v2, v1)
v0

merged_multiple = merge_multiple_sorted_vectors(vec1, vec2, vec3)
println("Merged multiple vectors: ", merged_multiple)
