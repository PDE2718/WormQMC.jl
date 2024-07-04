@generated function worm_cycle!(x::Wsheet{Ndim}, H::Ham,
    Q::UpdateConsts, Y::CycleAccumProb, G::T_G=nothing)::Int where {Ndim,Ham<:BH_Parameters,T_G<:Union{Nothing,GreenFuncBin}}

    @assert Ndim == N_wldim(H)
    znbs::Int = N_nbs(H)
    mes_green = (T_G == GreenFuncBin)
    # println("compiling update for $(Ham) (Ndim = $(Ndim), znbs = $(znbs)), GF = $(mes_green)")
    quote
        $(if mes_green quote G.insertion_trial += 1 end end)
        tail::Element = Element()
        head::Element = Element()
        head_new::Element = Element()
        β::f64 = x.β
        dice::f64 = 0.0
        δ::StateType = D::StateType = randsign() # relative direction from b to b̂ / other direction choice
        loc::WormLocation = _at_stop
        cycle_size::Int = 0
        #############################################################################
        ############################# [INSERT_WORM] #################################
        #############################################################################
        i::IndexType = j::IndexType = rand(eachindex(x))
        t::f64 = clamp(rand(), nextfloat(0.0, 2), prevfloat(1.0, 2))
        li::Wline = lj::Wline = x[i]
        head_id::Int = del_id::Int = vindex(li, t)
        n0::StateType = li[head_id].n_L
        np::StateType = n0 - D
        P_acc::f64 = 2 * Q.Cw * max(np, n0)
        if np < StateType(0) || np > H.nmax || !(metro(P_acc)) # case insertion failed!
            return 0
        end # otherwise accepted
        nbs::NTuple{$znbs,Int} = get_nbs(H, i)
        if close_to_any(li, t)
            return 0
        end
        @nexprs $znbs k -> begin
            if close_to_any(x[nbs[k]], t)
                return 0
            end
        end
        tail = Element(t, i, IndexType(0),
            D == i8(+1) ? n0 : np,
            D == i8(+1) ? np : n0,
            b_
        )
        head = Element(nextfloat(t, D), i, StateType(0), tail.n_R, tail.n_L, b̂_)
        insert!(li, head_id, tail)
        insert!(li, head_id + (D == i8(+1)), head)
        if rand(Bool) # swap head and tail
            head_new = tail
            tail = head
            head = head_new
            δ = -δ
        end
        cycle_size += 1

        ################################### begin worm cycle
        # @label CYCLE_START🔁
        for cycle_iter ∈ 1:100_000_000_000_000
            $(if mes_green quote accum_green!(G, tail, head, loc, H) end end)
            dice = rand()
            if dice < Y.AP_move_worm # [MOVE_WORM]
                D = randsign()
                i = head.i
                li = x[head.i]
                if D * δ == i8(-1)
                    if loc == _at_stop || loc == _at_kink
                        continue
                    elseif loc == _at_nbkink || loc == _at_green
                        # passby nb interaction && continue to move
                        head_id = vindex(li, head.t)
                        head <<= nextfloat(head.t, 2D)
                        li[head_id] = head
                        δ = -δ
                    else #if loc == _at_dummy
                        dummy::Element = li[end]
                        if δ == i8(+1) && D == i8(-1)
                            # @assert head == li[1] && head.t == nextfloat(0.)
                            head <<= prevfloat(1.)
                            dummy = dummy_element(dummy.i, head.n_R)
                            li[end] = head
                            push!(li, dummy)
                            popfirst!(li)
                        else  # if δ == -1 && D == +1
                            # @assert head == li[end-1] && head.t == prevfloat(1.)
                            head <<= nextfloat(0.0)
                            dummy = dummy_element(dummy.i, head.n_L)
                            pushfirst!(li, head)
                            pop!(li)
                            li[end] = dummy
                        end
                        δ = -δ
                    end
                end
                t = head.t
                head_id = vindex(li, head.t)
                v_near::Element = li[mod1(head_id + D, length(li))]
                nbs = get_nbs(H, i)
                @nexprs $znbs k -> begin
                    assoc_k::Element = element_around(x[nbs[k]], t, D)
                end
                nb_states::NTuple{$znbs,StateType} = if D == StateType(+1)
                    @ntuple $znbs k -> assoc_k.n_L
                else
                    @ntuple $znbs k -> assoc_k.n_R
                end
                λ₊::f64 = λ₋::f64 = ΔE::f64 = Q.Eoff
                Eb::f64 = diagE(H, head.n_L, nb_states)
                Ef::f64 = diagE(H, head.n_R, nb_states)
                if D == StateType(-1)
                    ΔE = Ef
                    Ef = Eb
                    Eb = ΔE
                end
                # [TODO] add check for illegal configs?
                if Ef > Ecutoff || Eb > Ecutoff
                    # @assert !(Ef > Ecutoff && Eb > Ecutoff) "problemastic update"
                    # @assert loc ≠ _at_free
                    continue
                end
                ΔE = abs(Ef - Eb)
                Δt = randexp()
                ##################################### try to move
                if D == StateType(+1) # move forward -->
                    @nexprs $znbs k -> begin
                        if assoc_k.t < v_near.t
                            v_near = assoc_k
                        end
                    end
                    # however, the head cannot pass the tail
                    if head.t < tail.t < v_near.t
                        v_near = tail
                        # @assert tail.i ≠ head.i # otherwise, we have tail == v_near already
                    end

                    if Ef ≥ Eb
                        λ₋ += ΔE
                    else
                        λ₊ += ΔE
                    end
                    Δt /= (β * λ₊)
                    if Δt < t_eps
                        continue
                    end
                    t_new = t + Δt
                    t_bound = v_near.t - t_eps
                    if (t_new < t_bound) && metro(δ == 0 ? (λ₋/λ₊) : inv(λ₊))
                        # no interaction
                        head <<= t_new
                        li[head_id] = head
                        δ = i8(0)
                        loc = _at_free
                        cycle_size += 1
                    # elseif (t_new ≥ t_bound) && metro(δ == 0 ? λ₋ : 1.)
                    elseif (t_new ≥ t_bound) && (δ ≠ 0 || metro(λ₋))
                        # interaction encountered
                        head <<= prevfloat(v_near.t)
                        li[head_id] = head
                        δ = i8(-1)
                        loc = if v_near.op == I_
                            _at_dummy
                        elseif v_near.j == IndexType(0) # v_near is the tail
                            v_near.i == head.i ? _at_stop : _at_green
                        else # v_near is a kink
                            v_near.i == head.i ? _at_kink : _at_nbkink
                        end::WormLocation
                        cycle_size += 1
                    end
                else # if D == -1 # move backward <--
                    if v_near.t == 1.0 # dummy element is at both 0 and β, here we regard it at 0.
                        v_near <<= 0.0
                    end
                    @nexprs $znbs k -> begin
                        if v_near.t < assoc_k.t < 1.0
                            v_near = assoc_k
                        end
                    end
                    # however, the head cannot pass the tail
                    if v_near.t < tail.t < head.t
                        v_near = tail
                        # @assert tail.i ≠ head.i # otherwise, we have tail == v_near already
                    end
                    if Ef ≥ Eb
                        λ₋ += ΔE
                    else
                        λ₊ += ΔE
                    end
                    Δt /= (β * λ₊)
                    if Δt < t_eps
                        continue
                    end
                    t_new = t - Δt
                    t_bound = v_near.t + t_eps
                    if (t_new > t_bound) && metro(δ == 0 ? (λ₋ / λ₊) : (1.0 / λ₊))
                        # no interaction
                        head <<= t_new
                        li[head_id] = head
                        δ = i8(0)
                        loc = _at_free
                        cycle_size += 1
                    # elseif (t_new ≤ t_bound) && metro(δ == 0 ? λ₋ : 1.)
                    elseif (t_new ≤ t_bound) && (δ ≠ 0 || metro(λ₋))
                        head <<= nextfloat(v_near.t)
                        li[head_id] = head
                        δ = i8(+1)
                        loc = if v_near.op == I_
                            _at_dummy
                        elseif v_near.j == IndexType(0) # v_near is the tail
                            v_near.i == head.i ? _at_stop : _at_green
                        else # v_near is a kink
                            v_near.i == head.i ? _at_kink : _at_nbkink
                        end::WormLocation
                        cycle_size += 1
                    end
                end
            elseif dice < Y.AP_insert_kink # [INSERT_KINK]
                if loc ≠ _at_free
                    continue
                end
                D = randsign()
                i = head.i
                nbs = get_nbs(H, i)
                j = rand(nbs) |> IndexType
                B = StateType(head.op)
                li = x[i]
                lj = x[j]
                jqL::Int, jqR::Int = vindex_around(lj, head.t)
                nj = lj[jqR].n_L # @assert nj == lj[jqL].n_R
                nmid = nj - B * D
                if nmid > H.nmax || nmid < i8(0)
                    continue
                end
                Wk = max(nj, nmid)
                P_acc = 2 * $znbs * Wk * bond_weight(H, i, j) * Q.P_del2ins
                if metro(P_acc)
                    # [TODO] check for order
                    nbs = get_nbs(H, j)
                    @nexprs $znbs k -> begin
                        if nbs[k] ≠ i && close_to_any(x[nbs[k]], head.t) # then insertion fail
                            continue
                        end
                    end
                    iq::Int = vindex(li, head.t) # @assert li[iq] == head
                    # first set head as part of K
                    li[iq] = Element(head.t, head.i, j, head.n_L, head.n_R, head.op)
                    # @set head.j = j
                    Kj = Element(head.t, j, i,
                        D == i8(+1) ? nj : nmid,
                        D == i8(+1) ? nmid : nj,
                        OpType(-B)
                    )
                    head_new = Element(nextfloat(head.t, D), j, IndexType(0), Kj.n_R, Kj.n_L, head.op)
                    if D == i8(+1)
                        insert!(lj, jqR, Kj)
                        insert!(lj, jqR + 1, head_new)
                    else
                        insert!(lj, jqR, head_new)
                        insert!(lj, jqR + 1, Kj)
                    end
                    # update worm
                    head = head_new
                    δ = D
                    loc = _at_kink
                    cycle_size += 1
                end

            elseif dice < Y.AP_delete_kink # [DELETE_KINK]
                if loc ≠ _at_kink
                    continue
                end
                i = head.i
                li = x[i]
                head_id = vindex(li, head.t)
                Ki_id::Int = head_id - δ
                Ki::Element = li[Ki_id]
                # @assert head.t == nextfloat(Ki.t, δ)
                if head.op == Ki.op # case 1 : pass kink with prob 1
                    li[head_id] = Element(Ki.t,
                        Ki.i, Ki.j, head.n_L, head.n_R, Ki.op
                    ) # now head_id is Ki
                    head_new = Element(nextfloat(Ki.t, -δ),
                        head.i, head.j, Ki.n_L, Ki.n_R, head.op
                    )
                    li[Ki_id] = head_new

                    # update worm
                    head = head_new
                    δ = -δ
                    loc = _at_kink
                    cycle_size += 1
                else # case 2 : delete_kink
                    j = Ki.j
                    lj = x[j]
                    Kj_id::Int = vindex(lj, Ki.t)
                    Kj::Element = lj[Kj_id]
                    Wk = max(head.n_L, head.n_R)
                    P_acc = inv(2 * $znbs * Wk * bond_weight(H, i, j) * Q.P_del2ins)
                    if metro(P_acc)
                        head_new = Element(Kj.t, Kj.i, IndexType(0), Kj.n_L, Kj.n_R, head.op)
                        lj[Kj_id] = head_new
                        del_id = min(Ki_id, head_id)
                        deleteat!(li, range(del_id, length=2))
                        # update worm
                        head = head_new
                        δ = i8(0)
                        loc = _at_free
                        cycle_size += 1
                    end
                end
            else #if glue worm # [DELETE_KINK]
                if loc == _at_stop && metro(inv(2 * Q.Cw * max(head.n_L, head.n_R)))
                    @assert head.i == tail.i && head.t == nextfloat(tail.t, δ)
                    li = x[head.i]
                    head_id = vindex(li, head.t)
                    del_id = min(head_id, head_id - δ)
                    deleteat!(li, range(del_id, length=2))
                    return cycle_size
                end
            end
        end
        error("loop size exceeds limit, something unusual happened!")
        return cycle_size
    end
end