using Worm
# H = BH_Pyroch(nmax=1,
#     Lx=8, Ly=8, U=0.0,
#     J1=1.0, J2=2.0,
#     V=20.0, μ=0.5
# )

H = BH_Square(nmax=1,
    Lx=8, Ly=8,
    U=+0.0, V=+20.0,
    μ=+0.5, J=1.0,
)

β = 20.0
update_const = UpdateConsts(0.5, 4.0, 1.0)
cycle_prob = CycleProb(1, 1, 1, 1)
time_ther = Second(20)
time_simu = Second(60)
x = Wsheet(β, H)
m = WormMeasure(x, update_const)
onesimu!(x, H, m, update_const, cycle_prob, time_ther, time_simu)

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