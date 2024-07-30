pwd()
@static if false
    include("./src/WormQMC.jl")
end
using Revise, WormQMC
using EDTools, LinearAlgebra, Statistics, SparseArrays, Dates, DelimitedFiles
using Logging

H = BH_Pyroch(nmax=1,
    Lx=16, Ly=16, U=0.0,
    J1=1.0, J2=2.0,
    V=10.0, μ=0.5
)

# H = BH_Square(nmax=1,
#     Lx=8, Ly=8,
#     U=+0.0, V=+1.0,
#     μ=+0.5, J=1.0,
# )

β = 16.0
update_const = UpdateConsts(0.5, 2.0, 1.0)
cycle_prob = CycleProb(1, 1, 1, 1)
time_ther = Second(10)
time_simu = Second(10)
x = Wsheet(β, H)
m = WormMeasure(x, update_const)
onesimu!(x::Wsheet{3}, H::BH_Pyroch, m::WormMeasure,
    update_const, cycle_prob, time_ther, time_simu)

m.Gfunc.insertion_trial


m.Gfunc |> normalize_density_matrix
m.Gfunc
function batched_cycles(x, H, m, update_const, cycle_prob, nbatch)
    total_size = 0
    for t ∈ 1:nbatch
        total_size += worm_cycle_!(x, H,
            update_const, cycle_prob, m.Gfunc
        )
    end
    return total_size
end
@time batched_cycles(x, H, m, update_const, cycle_prob, 100)
@time batched_cycles(x, H, m, update_const, cycle_prob, 100000)


mean(length.(x.wl)) ./ β

m.simple |> display
m.winding |> display
Cab(m.Sfact, 1, 1)
using Plots
heatmap(Cab(m.Sfact, 1, 1))

0.2029
mean(m.simple.N) / 64
m.Gfunc |> normalize_density_matrix
m.Gfunc

m.Gfunc.G0
m.simple.E.num
m.winding
m.Gfunc
res.winding


aaa = 1
quote
    $(
    aaa == 1 ? quote
        2 + yy
    end : quote
        3 + xx
    end
    )
    $(
        quote
            nothing
        end
    )
end |> typeof



Cij = Cab(res.Sfact, 1, 1)
using Plots
Cij |> heatmap

res.Sfact.Sk[1] |> full_Sk
function full_Sk(Sk)
    return vcat(Sk, Sk[end-1:-1:2, :])
end

heatmap(
    real.(res.Sfact.Sk[1]),
)

Gij = res.Gfunc |> normalize_density_matrix

heatmap(Gij[:, :, 2, 1])

inv.((1:16) .^ 2 .* real.(fft(Gij[:, :, 1, 1])[1, :]))

71911
res.simple
mean(res.simple.N) / (H.Lx * H.Ly) / 2
1 / 6
cal_Gτ(res.Gfunc, [0.0, 1.0])
(res |> Base.summarysize) / 8
(rand(100, 100) |> Base.summarysize) / 8

isbitstype(Element)
qq = [Element() for i ∈ 1:3]
sizeof(Element())
sizeof(qq)

reinterpret(Tuple{Int64,Int64}, qq)
reinterpret(Float64,reinterpret(NTuple{15, UInt8}, Element() << -1000.)[1:8])

using FieldFlags
@bitfield mutable struct foo
    a:64
    b:8
    c:16
    d:1
end
@bitfield struct bar
    a:64
    b:8
    c:16
    d:1
end

sizeof()
xx = foo(10, 3, 2, 1)
@btime $(xx).a = rand(UInt64)

@btime Float64(MyBits(1,2).a)
using BenchmarkTools

