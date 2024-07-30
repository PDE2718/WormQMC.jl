abstract type Node_Z end
for N ∈ 1:2
    local x = Meta.parse("Node_$(N)")
    Meta.parse(
        "@kwdef mutable struct $(x) <: Node_Z \n" *
        "prev::Ptr{$(x)} \n" * 
        prod(ntuple(i->"    s_$(i)::Ptr{$(x)} \n", N)) *
        "t::Float64 \n" *
        "end"
    ) |> eval
end


mutable struct NodeZ{Z}
    nbs::NTuple{Z, Ptr{Nothing}}
    self::Ptr{Nothing}
    prev::Ptr{Nothing}
    next::Ptr{Nothing}
    t::Float64
    function NodeZ{Z}() where {Z}
        x = new{Z}()
        x.self = pointer_from_objref(x)
        return x
    end
    function NodeZ{Z}(t::Float64) where {Z}
        x = new{Z}()
        x.t = t
        return x
    end
end
aa = NodeZ{3}()
unsafe_load(aa |> pointer_from_objref |> Ptr{NodeZ{3}})
import Base:getindex, setindex!
function getindex(x::NodeZ{Z}, i::R)::NodeZ{Z} where {Z, R<:Integer}
    return unsafe_load(x.nbs[i] |> Ptr{NodeZ{Z}})
end
function setindex!(x::NodeZ{Z}, nb::NodeZ{Z}, i::R) where {Z, R<:Integer}
    unsafe_store!(Ptr{Ptr{Nothing}}(x.self), nb |> pointer_from_objref, i)
    return nb
end
()
aa[1] = aa


mutable struct MyFoo{N}
    t::Float64
    prev::Ptr{MyFoo{N}}
    next::Ptr{MyFoo{N}}
    const nbs::NTuple{N, Ptr{MyFoo{N}}}
    function MyFoo{N}() where {N}
        x = new{N}()
        return x
    end
    function MyFoo{N}(tag) where {N}
        x = new{N}()
        x.tag = d
        return x.tag
    end
end



Foo1(1,1.0)
Foo2(1,2.0,2.0)

1

for NN ∈ 1:3
    eval(
end
MyFoo1(1)


:(x_ = 3)

mutable struct MyFoo{N}
    tag::Int
    prev::RefValue{MyFoo{N}}
    next::RefValue{MyFoo{N}}
    nbs::NTuple{N,RefValue{MyFoo{N}}}
    function MyFoo{N}() where {N}
        x = new{N}()
        return x
    end
    function MyFoo{N}(tag) where {N}
        x = new{N}()
        x.tag = d
        return x.tag
    end
end

function set_nb(s::MyFoo{N}, i::Int)
    
end
function get_nb()
    
end
pointer_from_objref
sizeof(MyFoo{4}())
MyFoo{10}()
aa = MyFoo{3}()
aa.tag = 114
aa.prev = pointer_from_objref(aa)

pointer_from_objref(aa)
unsafe_load(Ptr{typeof(aa)}(pointer_from_objref(aa)+100))

aa.prev = aa |> pointer_from_objref
unsafe_load(aa.nbs[1])
load


Base.RefValue{Int}(1)[]



sizeof(MyFoo{2}())

sizeof(MyFoo{4}())

MyFoo{4}()

MVector{2, MyFoo}()

(MyFoo(), MyFoo())

SVector()
aa = MyFoo()
aa.m = pointer_from_objref(aa)
sizeof(aa)

using BenchmarkTools
@btime unsafe_load($(aa.m))
@btime aa.m

finalize(aa)

unsafe_load(Ptr{Int32}(pointer_from_objref(aa)) + 8)
unsafe_load(Ptr{Int32}(pointer_from_objref(aa)) + 8)



unsafe_load(Ptr{MyFoo}(Ptr{Nothing}()))

function test_foo(a, b)
    f = MyFoo(0, 0.)
    f.a = a
    f.b = b
    f.b += rand()
    return f.a + f.b
end
function test_small_array(a, b)
    m = [a,b]
    return sum(m)
end

using BenchmarkTools
@benchmark test_foo($(1), $(2.0))
@benchmark test_small_array(1,2)


eval(Expr(:struct, true, :(MNT200{$((Symbol(:T, i) for i in 1:200)...)}), Expr(:block, ((:($(Symbol('a' + i -1)) :: $(Symbol(:T,i))) for i in 1:200 ))...)))

@generated function Base.NamedTuple(mnt::MNT200)
           Expr(:call, NamedTuple{fieldnames(MNT200)}, Expr(:tuple, (:(getfield(mnt, $i)) for i ∈ 1:200)...))
       end;

@generated function MNT200(nt::NamedTuple)
           Expr(:call, MNT200, (:(nt[$i]) for i ∈ 1:200)...)
       end;

       N = 200
nt = NamedTuple{ntuple(i -> Symbol('a' + i - 1), N)}(ntuple(i -> i == 2 ? rand(Int) : rand(("hi", 1 + im, [1, 2], Ref{Any}(1))), N))
Base.RefValue
let N = 200
           mnt = MNT200(nt)
           @btime $mnt.b = 10
end;
