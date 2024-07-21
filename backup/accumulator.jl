import Statistics: mean
mutable struct Accum{T}
    val::T
    num::Int
end
function Accum(x)
    return Accum(zero(x), 0)
end
import Base: push!, append!, empty!
function push!(A::Accum{T}, x::T) where {T<:Number}
    A.val += x
    A.num += 1
    A
end
function push!(A::Accum{T1}, x::T2) where {T1<:AbstractArray,T2<:AbstractArray}
    A.val .+= x
    A.num += 1
    A
end
function append!(A::Accum{T}, xs::Vector{T}) where {T}
    for x ∈ xs
        push!(A, x)
    end
    A
end
function empty!(A::Accum{T}) where {T<:Number}
    A.val = zero(T)
    A.num = 0
    A
end
function empty!(A::Accum{T}) where {T<:AbstractArray}
    A.val .= zero(eltype(T))
    A.num = 0
    A
end
function mean(A::Accum{T}) where {T}
    return A.val ./ A.num
end
import Base:merge
function merge(x1::Accum{T}, x2::Accum{T}) where T
    return Accum(x1.val + x2.val, x1.num + x2.num)
end
import Base.show
function Base.show(io::IO, a::Accum{T}) where T
    if a.num == 0
        print(io, typeof(a), " empty")
    else
        print(io, typeof(a), " mean=", a.val ./ a.num)
    end
end
