using BenchmarkTools
using StaticArrays
using Distributions

struct myfoo
    a::Int16
end
@btime myfoo(Int32(10))

itype = UInt128
A = rand(Float64, 64)
B = rand(itype, 100)
randu128() = rand(UInt128)


function test_plus(A, B)
    A .-= B
end

@btime test_plus($A,$B)

@btime UInt($1) + UInt(1)
@btime UInt128($1) + UInt128(1)
using Quadmath
@btime rand(UInt64) + 1
@btime rand(Float128) + 1.0
@btime rand(Float64) + 1.0
nextfloat(45.0,1)-45.0
nextfloat(5.0, 1)-5.0
using Random
@btime Random.rand!($A)
@btime for i ∈ 1:1000
    randu128()
end

randexp!(A)
rand!(A)
A

using Distributions
Distributions.
@btime for i ∈ 1:1000
    rand(Exponential(2.0))
end
@btime for i ∈ 1:1000
    randexp()
end

using Base.Cartesian.Base:@nexprs
@nexprs 4 i->begin
    println(i)
end

@btime rand(Uniform(nextfloat(0.), prevfloat(10.)))
@btime rand()
myrng = Random.MersenneTwister(123)
myrng = Random.XoshiroSimd
@btime rand($myrng)

using Base.Cartesian.Base:@nextract, @nexprs

qq = 2 .+ (1:6)
@nextract 6 xx i -> cos(i)
@nextract 6 xx qq
@nexprs 6 begin
    
end
@macroexpand 

Base.Cartesian.@nexprs 4 i -> y_i = i
Base.@nexprs 3 i -> begin 
    assoc_i = Element()
end

@macroexpand Base.@ntuple 2 i->assoc_i.n_L
sizeof(:n_L)

Base.@ntuple 2 i->2
-log(prevfloat(1.))
nextfloat(0.)


@macroexpand @nexprs 4 inb -> begin
    if assoc_inb.t < v_near.t
        v_near = assoc_i
    end
end

xx = if false
    1
else
    2
end::Int

assoc_1


y_1
y_2

xx_1


xx_7
xx_2
xx_3
using Base:@nexprs
@generated function myfoo(a::NTuple{z,Int})::Float64 where z
    println("z=$(z)")
    println("a=$(a)")
    return quote
        s = 0
        Base.@nexprs $z i->begin
            s += a[i]
        end
        return s
    end
end
@generated function mybar(u::Real, a::NTuple{z,Int}) where z
    println("z=$(z)")
    println("a=$(a)")
    return quote
        s = 0
        Base.@nexprs $z i->begin
            s += a[i]
        end
        return s + u
    end
end
using BenchmarkTools
qq = Base.@ntuple 20 i-> 2i
myfoo(qq) |> typeof
mybar(1 |> Int32,qq)
mybar
methods(mybar)

@code_warntype myfoo(qq)

using BenchmarkTools
@btime myfoo($qq)
@btime sum($qq)

myfoo(::Complex) = 10
myfoo(::Int) = 20
myfoo(1+1im)


AA = rand(1000)
function bench_deleteat!(AA)
    deleteat!(AA,1)
    push!(AA,1.)
end
function bench_pop!(AA)
    a = pop!(AA)
    pushfirst!(AA)
end
using BenchmarkTools
@benchmark bench_deleteat!($AA)
@benchmark bench_deleteat!($AA)


xx = Element()
xx <<= 10.
xx


@generated function mysum(a::NTuple{z,T}) where {z, T<:Real}
    println("compiling")
    return quote
        s = 0
        Base.@nexprs $z i->begin
            s += a[i]
        end
        return s
    end
end
using Base:@ntuple
A = @ntuple 1000 i -> Float64(i)
using BenchmarkTools

@benchmark sum($(A))
@benchmark mysum($(A))
Base.JLOptions().opt_level

Base.@nexprs 5 i->begin
    xx_i = i
end


using Random
using RandomNumbers
using BenchmarkTools
# const 
const myrng::Xoshiro = Random.Xoshiro()
myrng2 = RandomNumbers.()
Random.default_rng() |> typeof
@btime Random.default_rng()
function randfill1(A, rng=Random.default_rng())
    for i ∈ eachindex(A)
        A[i] = rand(rng)
    end
end

A = rand(1000)
@btime randfill1($A)
@btime randfill1($A, $myrng2)


using LegendrePolynomials
using BenchmarkTools
aaa = rand(100)

@btime collectPl!($aaa, $0.44, norm=Val(:schmidt))

xs = [0.1,0.1,0.2,0.2,0.3,0.4]

aaa = rand(100)
@btime Pl_schmidt_lazy!($aaa, $1.0)

0.44 * √3
