####### Type alias ######
const f64::Type = Float64
const f32::Type = Float32
const i32::Type = Int32
const i16::Type = Int32
const i8::Type = Int8
const StateType::Type = Int8
const IndexType::Type = Int16
@inline metro(p)::Bool = rand() < p
@inline randsign()::Int8 = rand(Bool) ? Int8(1) : Int8(-1)

function dumptoNamedTuple(x)
    keys = fieldnames(typeof(x))
    vals = [getfield(x, k) for k ∈ keys]
    return (; zip(keys, vals)...)
end

####### Element ######
@enum OpType::Int8 begin
    b_ = -1 # annihilation
    I_ = 0  # Identity (e.g. for dummy)
    b̂_ = +1 # creation
end
struct Element
    t::f64          # time of the element
    i::IndexType    # the site of element it self.
    j::IndexType    # the other site of a bond, 0 for the worm head or tail.
    n_L::StateType  # state to the left  of the element.
    n_R::StateType  # state to the right of the element.
    op::OpType      # the operator type
end
Element()::Element = Element(0.,IndexType(0), IndexType(0), i8(0), i8(0), I_)
function Element(op::OpType, t::f64,
    i::IndexType, n_L::StateType, n_R::StateType,
    j::IndexType=IndexType(0))::Element
    return Element(t, i, j, n_L, n_R, op)
end
const ElePair::Type = Pair{Element,Element} # two related Elements
const Wline::Type = Vector{Element}
# Note that : For each wl the last element must be I_ at β (dummy element)
empty_wline(i0::Integer, β::f64, n0::Integer) = Element[Element(I_, β, IndexType(i0),  StateType(n0), StateType(n0))]
is_bond(x::Element)::Bool = iszero(e.j)

function display_Wline(l::Wline)
    println("Wline : i = $(l[end].i), β = $(l[end].t)")
    for e ∈ l
        println("ψ $(e.n_L)→$(e.n_R), op:$(e.op), t=$(e.t)")
    end
    println(" ")
end

import Base: isless, isequal, ==
isless(a::Element, b::Element)::Bool = isless(a.t, b.t)
isless(a::Element, t::f64)::Bool = isless(a.t, t)
isless(t::f64, b::Element)::Bool = isless(t, b.t)
vindex(l::Wline, t::f64)::Int = searchsortedfirst(l, t)

# vertex id around t. if at = true, then exlude the vertex at t (it must exist).
@inline function vindex_around(l::Wline, t::f64, at::Bool = false)::Pair{i32, i32}
    # @static if worm_debug
    #     @assert 0 < t ≤ last(l).t
    # end
    len::i32 = lastindex(l)
    i::i32 = vindex(l, t)
    i_R::i32 = mod1(i+i32(at), len)
    i_L::i32 = mod1(i-i32(1), len)
    return i_L=>i_R
end
@inline function vindex_around(l::Wline, t::f64, D::Integer, at::Bool=false)::i32
    # @static if worm_debug
    #     @assert 0 < t ≤ last(l).t
    # end
    len::i32 = lastindex(l)
    i::i32 = vindex(l, t)
    if D > 0
        i = mod1(i + i32(at), len)
    else
        i = mod1(i - i32(1), len)
    end
    return i
end
@inline function element_around(l::Wline, t::f64, at::Bool=false)::ElePair
    # @static if worm_debug
    #     @assert 0 < t ≤ last(l).t
    # end
    len::i32 = lastindex(l)
    i::i32 = vindex(l, t)
    i_R::i32 = mod1(i + i32(at), len)
    i_L::i32 = mod1(i - i32(1), len)
    return l[i_L] => l[i_R]
end
@inline function element_around(l::Wline, t::f64, D::Integer, at::Bool=false)::Element
    # @static if worm_debug
    #     @assert 0 < t ≤ last(l).t
    # end
    len::i32 = lastindex(l)
    i::i32 = vindex(l, t)
    if D > 0
        i = mod1(i + i32(at), len)
    else
        i = mod1(i - i32(1), len)
    end
    return l[i]
end
element_le(l::Wline, t::f64)::Element = getindex(l, vindex_around(l,t).first)
element_geq(l::Wline, t::f64)::Element = getindex(l, vindex_around(l,t).second)
function erase_vertex!(l::Wline, e::Element)::Nothing
    # @static if worm_debug
    #     @assert e.t ≠ last(l).t "Dummy element cannot be erased"
    # end
    deleteat!(l, vindex(l, e.t))
    nothing
end

## Here goes definitions of world sheet and worm!
## The Wsheet and worm decouples!
struct Wsheet{N}
    β::f64
    wl::Array{Wline,N}
end
function Wsheet(β::f64, ψ0::AbstractArray{<:Integer,N})::Wsheet where {N}
    @assert β > 0
    return Wsheet(β,[empty_wline(i, β, ψ0[i]) for i ∈ LinearIndices(ψ0)])
end
import Base: getindex, eachindex, size, LinearIndices, CartesianIndices
function getindex(X::Wsheet{N}, i::Integer)::Wline where {N}
    return getindex(X.wl, i)
end
function getindex(X::Wsheet{N}, inds...)::Wline where {N}
    return getindex(X.wl, inds...)
end
eachindex(X::Wsheet) = eachindex(X.wl)
LinearIndices(X::Wsheet) = LinearIndices(X.wl)
CartesianIndices(X::Wsheet) = CartesianIndices(X.wl)
size(X::Wsheet) = size(X.wl)

@enum WormLocation::Int8 begin
    _at_null # null worm tag
    _at_free # at free time on a world line
    _at_stop # head just before or after (near) tail
    _at_kink # head is near a kink
    _at_nbkink # head is near a neighboring kink. Usually can pass
    _at_dummy # head is before β or after 0. near the dummy element
    _at_green # head is of the same time as tail, but not at the same site.
end
struct Worm
    tail::Element
    head::Element
    δ::StateType # head is before(left,δ=-1) or after(right,δ=+1) another element
    loc::WormLocation # see above WormLocation
end
Worm()::Worm = Worm(Element(),Element(),StateType(0),_at_null)