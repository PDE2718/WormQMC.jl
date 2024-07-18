function onesimu!(x::Wsheet{N}, H::Ham, m::WormMeasure,
    update_const::UpdateConsts,
    cycle_prob::CycleAccumProb,
    time_ther::Second,
    time_simu::Second,
    sweep_size_limit::Tuple{Int,Int} = (1,100),
    )::Nothing where {N, Ham <: BH_Parameters}
    # prepare buffers measurements
    @assert m.Gfunc.Cw == update_const.Cw "Gfunc Cw is not consistent with the update_const!"
    bond_buffer = Element[]
    Sk_plan = m.Sfact |> rfftplan
    N_cycle::Int = total_cycle_size::Int = 0
    sweep_size::Int = 100
    
    @info "thermalizing for $(time_ther)"
    @assert iszero(time_ns() >> 62) "current time_ns is close to wrap and can produce bugs!"
    T_ther::UInt64 = Nanosecond(time_ther).value
    T_simu::UInt64 = Nanosecond(time_simu).value
    t_limit::UInt64 = time_ns() + T_ther
    while time_ns() < t_limit
        for _ ∈ 1:sweep_size
            total_cycle_size += worm_cycle!(x, H, update_const, cycle_prob)
            N_cycle += 1
        end
    end
    average_cycle_size::f64 = total_cycle_size / N_cycle
    sweep_size = ceil(Int, sum(length, x.wl) / average_cycle_size)
    sweep_size = clamp(sweep_size, sweep_size_limit...)

    @info """Thermalization Statistics
    [total_cycle_size] = $(total_cycle_size)
    [N_cycles        ] = $(N_cycle)
    [average_size    ] = $(total_cycle_size/N_cycle)
    [wl length / β   ] = $((mean(length, x.wl)-1) / x.β)
    [N cycle per mes ] = $(sweep_size)
    """
    @assert all(issorted, x.wl)
    
    N_cycle = total_cycle_size = 0
    @info "simulating for $(time_simu)"
    T_measure::UInt64 = 0
    t_limit = time_ns() + T_simu
    while time_ns() < t_limit
        for _ ∈ 1:sweep_size
            total_cycle_size += worm_cycle!(x, H, update_const, cycle_prob, m.Gfunc)
            N_cycle += 1
        end
        T_measure += measure!(m, x, H, bond_buffer, Sk_plan)
    end

    @assert all(issorted, x.wl)
    average_cycle_size = total_cycle_size / N_cycle
    @info """[FINISHED] Simulation Statistics
    [total_cycle_size ] = $(total_cycle_size)
    [N_cycles         ] = $(N_cycle)
    [average_size     ] = $(total_cycle_size/N_cycle)
    [wl length / β    ] = $((mean(length, x.wl)-1) / x.β)
    [total measure    ] = $(m.simple.n_measure)
    [total   time (ns)] = $(T_simu)
    [measure time (ns)] = $(T_measure)
    [measurement]: $(m.simple)
    """
    return nothing
end

function onesimu!(x::Wsheet{N}, H::Ham, m::WormMeasure, ψsnaps::Vector{Array{i8,N}},
    update_const::UpdateConsts,
    cycle_prob::CycleAccumProb,
    time_ther::Second,
    time_simu::Second,
    sweep_size_limit::Tuple{Int,Int}=(1, 256),
    save_snap_every::Int=128,
    max_snap::Int=1024,
)::Nothing where {N, Ham<:BH_Parameters}
    # prepare buffers measurements
    @assert m.Gfunc.Cw == update_const.Cw "Gfunc Cw is not consistent with the update_const!"
    bond_buffer = Element[]
    Sk_plan = m.Sfact |> rfftplan
    N_cycle::Int = total_cycle_size::Int = 0
    sweep_size::Int = 100

    @info "thermalizing for $(time_ther)"
    @assert iszero(time_ns() >> 62) "current time_ns is close to wrap and can produce bugs!"
    T_ther::UInt64 = Nanosecond(time_ther).value
    T_simu::UInt64 = Nanosecond(time_simu).value
    t_limit::UInt64 = time_ns() + T_ther
    while time_ns() < t_limit
        for _ ∈ 1:sweep_size
            total_cycle_size += worm_cycle!(x, H, update_const, cycle_prob)
            N_cycle += 1
        end
    end
    average_cycle_size::f64 = total_cycle_size / N_cycle
    sweep_size = ceil(Int, sum(length, x.wl) / average_cycle_size)
    sweep_size = clamp(sweep_size, sweep_size_limit...)

    @info """Thermalization Statistics
    [total_cycle_size] = $(total_cycle_size)
    [N_cycles        ] = $(N_cycle)
    [average_size    ] = $(total_cycle_size/N_cycle)
    [wl length / β   ] = $((mean(length, x.wl)-1) / x.β)
    [N cycle per mes ] = $(sweep_size)
    """
    @assert all(issorted, x.wl)

    N_cycle = total_cycle_size = 0
    @info "simulating for $(time_simu)"
    T_measure::UInt64 = 0
    t_limit = time_ns() + T_simu
    while time_ns() < t_limit
        for _ ∈ 1:sweep_size
            total_cycle_size += worm_cycle!(x, H, update_const, cycle_prob, m.Gfunc)
            N_cycle += 1
            # For taking snaps
            if total_cycle_size % save_snap_every == 0
                if length(ψsnaps) < max_snap
                    push!(ψsnaps, map(l->first(l).n_L, x.wl))
                else
                    save_snap_every *= 2
                    deleteat!(ψsnaps, 1:2:length(ψsnaps))
                end
            end
        end
        T_measure += measure!(m, x, H, bond_buffer, Sk_plan)
    end
    @assert all(issorted, x.wl)
    average_cycle_size = total_cycle_size / N_cycle
    @info """[FINISHED] Simulation Statistics
    [total_cycle_size ] = $(total_cycle_size)
    [N_cycles         ] = $(N_cycle)
    [average_size     ] = $(total_cycle_size/N_cycle)
    [wl length / β    ] = $((mean(length, x.wl)-1) / x.β)
    [total measure    ] = $(m.simple.n_measure)
    [total   time (ns)] = $(T_simu)
    [measure time (ns)] = $(T_measure)
    [measurement]: $(m.simple)
    """
    return nothing
end