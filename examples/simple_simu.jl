using Dates
function simple_simu(H::BH_Parameters, β::f64;
    update_const::UpdateConsts = UpdateConsts(0.5, 4.0, 1.0),
    cycle_prob::CycleAccumProb = CycleProb(1, 1, 1, 1),
    time_simu::TimePeriod=Second(5),
    time_ther::TimePeriod=Second(5),
    mes_every::Int = 2,
    green_lmax::Int = 100
)
    x::Wsheet = Wsheet(β, H)
    # prepare measurements
    bond_buffer = Element[]
    mes_simple = SimpleMeasure()
    mes_Gfunc = GreenFuncBin(x, green_lmax, update_const.Cw)
    mes_Sfact = StructureFactor2D(H)
    Sk_rFFTWPlan = mes_Sfact |> rfftplan
    mes_winding = (Accum(0), Accum(0), Accum(0))
    average_cycle_size::f64 = 0.
    cycle_size::Int = insertion_trial::Int = insertion_success::Int = 0
    N_cycle::Int = N_measure_simple::Int = 0

    println("thermalizing for $(time_ther)")
    time_limit::f64 = time() + (time_ther |> Second).value    
    total_cycle_size::Int = 0
    while time() < time_limit
        for _ ∈ 1:mes_every
            total_cycle_size += worm_cycle!(x, H, update_const, cycle_prob)
            N_cycle += 1
        end
    end
    average_cycle_size = total_cycle_size / N_cycle
    println("Thermalization info:")
    println("average_cycle_size = $(average_cycle_size)")
    println("total_cycle = $(N_cycle)")
    mes_every = clamp(
        ceil(Int, (H.Lx * H.Ly * β) / average_cycle_size),
        mes_every, 10mes_every
    )

    N_cycle = total_cycle_size = 0
    println("simulating for $(time_simu)")
    time_limit = time() + (time_simu |> Second).value

    while time() < time_limit
        for _ ∈ 1:mes_every
            cycle_size = worm_cycle!(x, H,
                update_const, cycle_prob, mes_Gfunc)
            total_cycle_size += cycle_size
            insertion_trial += 1
            insertion_success += (cycle_size>0)
            N_cycle += 1
        end
        measure!(mes_simple, simple_measure(x, H, bond_buffer))
        measure_Sk2D!(mes_Sfact, x, Sk_rFFTWPlan)
        Wx, Wy=winding_number(x, H)
        push!.(mes_winding, (abs2(Wx), abs2(Wy), abs(Wx * Wy)))
        N_measure_simple += 1
    end
    println("[finished]")
    println("average WL length / β =", mean(length, x.wl)/β)
    @assert all(issorted, x.wl)
    return (
        mes_simple = mes_simple,
        mes_Gfunc = mes_Gfunc,
        mes_Sfact = mes_Sfact,
        mes_winding = mes_winding,
        average_cycle_size = average_cycle_size,
        mes_every = mes_every,
        N_measure_simple = N_measure_simple,
        insertion_trial = insertion_trial,
        insertion_success = insertion_success,
    )
    # update_statistics = (
    #     average_cycle_size=average_cycle_size,
    #     mes_every=mes_every,
    #     N_measure_simple=N_measure_simple,
    # )
    # return (
    #     mes_simple |> dumptoNamedTuple,
    #     mes_Gfunc,
    #     mes_Sfact,
    #     mes_winding,
    #     update_statistics,
    # )
end
