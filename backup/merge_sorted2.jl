using Base: @ntuple, @nexprs
using StaticArrays
MVector()
SVector
@ntuple 3 i->0
@generated function merge_sorted_oN!(dst::Vector{T}, vs::NTuple{N, Vector{T}})::Nothing where {N,T}
    if N == 1
        quote
            empty!(dst)
            append!(dst, first(vs))
            # @assert issorted(v0)
            return nothing
        end
    else
        quote
            # empty!(dst)
            L0::Int = 0
            t = MVector{$(N),Int}(@ntuple $(N) i -> 1)
            L = SVector{$(N),Int}(@ntuple $(N) i -> length(vs[i]))
            L0 = sum(L)
            resize!(dst, L0)
            local min_idx::Int
            local min_val::T
            @inbounds for t_dst ∈ 1:L0
                min_idx = 0
                @nexprs $(N) i->begin
                    if t[i] ≤ L[i] && (min_idx == 0 || vs[i][t[i]] < min_val)
                        min_val = vs[i][t[i]]
                        min_idx = i
                    end
                end
                dst[t_dst] = min_val
                t[min_idx] += 1
            end
            return nothing
        end
    end
end

@generated function merge_sorted_cache!(dst::Vector{T}, vs::NTuple{N,Vector{T}}) where {N,T}
    if N == 1
        println("N==1 special case")
        quote
            empty!(dst)
            append!(dst, first(vs))
            return nothing
        end
    else
        quote
            # empty!(dst)
            L0::Int = 0
            @nexprs $(N) i -> begin
                L_i::Int = length(vs[i])
                t_i::Int = 1
                local e_i::T
                L0 += L_i
            end
            resize!(dst, L0)
            local min_idx::Int
            local min_val::T
            @nexprs $(N) i -> begin
                if t_i ≤ L_i
                    e_i = first(vs[i])
                end
            end
            @inbounds for t_dst ∈ 1:L0
                min_idx = 0
                @nexprs $(N) i -> begin
                    if t_i ≤ L_i
                        if min_idx == 0 || e_i < min_val
                            min_val = e_i
                            min_idx = i
                        end
                    end
                end
                dst[t_dst] = min_val
                @nexprs $(N) i -> if i == min_idx
                    t_i += 1
                    if t_i ≤ L_i
                        e_i = vs[i][t_i]
                    end
                end
            end
        end
    end
end

using WormQMC, Accessors
v0 = Element[]
# v0 = Float64[]
v_mat = [ [Element() << t for t ∈ sort(rand(50))] for i ∈ 1:100, j ∈ 1:100]
vs = (rand(v_mat, 3)...,)

# vs = ntuple(i -> sort(rand(10000)), 8)
using BenchmarkTools
merge_sorted_oN!(v0, vs)
issorted(v0)
@btime merge_sorted_oN!($v0, $(vs))
@code_warntype merge_sorted_oN!(v0, vs)
@btime sort!(vcat($(vs)...))
@btime merge_sorted_cache!($v0, $(vs))



# issorted(v0)



# @btime merge_sorted!($v0, $(vs[1]), $(vs[2]))
# issorted(v0)



# @btime merge_sorted_oN!($v0, $(vs[1]), $(vs[2]))
# @btime merge_sorted!($v00, $(vs[1]), $(vs[2]))


# issorted(v0)
# # @code_warntype merge_sorted_oN!(v0, vs)
# @btime naive_merge_sorted!($v0, $(vs))


# v0 |> issorted
# v0
# v0
# v0 |> issorted

# merged_multiple = merge_multiple_sorted_vectors(vec1, vec2, vec3)
# println("Merged multiple vectors: ", merged_multiple)
using WormQMC, Accessors
v0 = Element[]
# v0 = Float64[]
v_mat = [[Element() << t for t ∈ sort(rand(50))] for i ∈ 1:100, j ∈ 1:100]
vs = (rand(v_mat, 3)...,)

# vs = ntuple(i -> sort(rand(10000)), 8)
using BenchmarkTools
@btime merge_sorted_oN!($v0, $(vs))
@btime sort!(vcat($(vs)...))
@btime merge_sorted_cache!($v0, $(vs))



# issorted(v0)



# @btime merge_sorted!($v0, $(vs[1]), $(vs[2]))
# issorted(v0)



# @btime merge_sorted_oN!($v0, $(vs[1]), $(vs[2]))
# @btime merge_sorted!($v00, $(vs[1]), $(vs[2]))


# issorted(v0)
# # @code_warntype merge_sorted_oN!(v0, vs)
# @btime naive_merge_sorted!($v0, $(vs))


# v0 |> issorted
# v0
# v0
# v0 |> issorted

# merged_multiple = merge_multiple_sorted_vectors(vec1, vec2, vec3)
# println("Merged multiple vectors: ", merged_multiple)
@generated function merge_sorted_cache!(dst::Vector{T}, vs::NTuple{N,Vector{T}}) where {N,T}
    if N == 1
        println("N==1 special case")
        quote
            empty!(dst)
            append!(dst, first(vs))
            return nothing
        end
    else
        quote
            # empty!(dst)
            L0::Int = 0
            @nexprs $(N) i -> begin
                L_i::Int = length(vs[i])
                t_i::Int = 1
                local e_i::T
                L0 += L_i
            end
            resize!(dst, L0)
            local min_idx::Int
            local min_val::T
            @nexprs $(N) i -> begin
                if t_i ≤ L_i
                    e_i = first(vs[i])
                end
            end
            @inbounds for t_dst ∈ 1:L0
                min_idx = 0
                @nexprs $(N) i -> begin
                    if t_i ≤ L_i
                        if min_idx == 0 || e_i < min_val
                            min_val = e_i
                            min_idx = i
                        end
                    end
                end
                dst[t_dst] = min_val
                @nexprs $(N) i -> if i == min_idx
                    t_i += 1
                    if t_i ≤ L_i
                        e_i = vs[i][t_i]
                    end
                end
            end
        end
    end
end