## Some update constants
@kwdef struct UpdateConsts
    Cw::f64 = 0.5
    Eoff::f64 = 1.0
    P_del2ins::f64 = 1.0 # ratio p_delete/p_insert
end
@kwdef struct CycleAccumProb
    AP_move_worm::f64   = 0.25
    AP_insert_kink::f64 = 0.50
    AP_delete_kink::f64 = 0.75
end
import Base.show
function Base.show(io::IO, a::CycleAccumProb)
    print(io, """Worm Cycle Probabilities
    P_move   = $(a.AP_move_worm - 0.0)
    P_insert = $(a.AP_insert_kink - a.AP_move_worm)
    P_delete = $(a.AP_delete_kink - a.AP_insert_kink)
    P_glue   = $(1.0 - a.AP_delete_kink)
    """
    )
end
function CycleProb(P_move_worm, P_insert_kink, P_delete_kink, P_glue_worm)::CycleAccumProb
    P_tot = P_move_worm + P_insert_kink + P_delete_kink + P_glue_worm
    p = CycleAccumProb(
        P_move_worm / P_tot,
        (P_move_worm + P_insert_kink) / P_tot,
        (P_move_worm + P_insert_kink + P_delete_kink) / P_tot,
    )
    @assert 0. < p.AP_move_worm < p.AP_insert_kink < p.AP_delete_kink < 1.
    return p
end

function insert_worm!(x::Wsheet, w::Worm,
    H::BH_Parameters{znbs}, Q::UpdateConsts)::Worm where {znbs}
    if w.loc ≠ _at_null
        return w
    end
    β::f64 = x.β
    i::IndexType = rand(eachindex(x))
    t::f64 = rand(Uniform(nextfloat(0., +2), nextfloat(β, -2)))
    l::Wline = x[i]
    i_after::Int = vindex(l, t)
    e_after::Element = l[i_after]
    n0::StateType = e_after.n_L
    D::StateType = randsign() # relative direction for b to b̂
    np::StateType = n0 - D
    P_acc::f64 = 2 * Q.Cw * max(np, n0)
    if np < StateType(0) || np > H.nmax || !(metro(P_acc)) # case insertion failed!
        return w
    end # otherwise accepted
    nbs::NTuple{znbs,Int} = get_nbs(H, i)
    if close_to_any(l, t)
        return w
    end
    for inb ∈ nbs
        if close_to_any(x[inb], t)
            return w
        end
    end
    b = Element(t, i, IndexType(0),
        D == i8(+1) ? n0 : np,
        D == i8(+1) ? np : n0,
        b_
    )
    b̂ = Element(nextfloat(t, D), i, StateType(0), b.n_R, b.n_L, b̂_)
    insert!(l, i_after, b)
    insert!(l, i_after + (D == i8(+1)), b̂)
    if rand(Bool)
        return Worm(b, b̂, +D, _at_stop)
    else
        return Worm(b̂, b, -D, _at_stop)
    end
end

function glue_worm!(x::Wsheet, w::Worm,
    H::BH_Parameters{znbs}, Q::UpdateConsts)::Worm where {znbs}
    if w.loc == _at_stop && metro(inv(2 * Q.Cw * max(w.head.n_L, w.head.n_R))) # accept ratio
        # @assert w.δ == +1 || w.δ == -1
        # @assert w.head.i == w.tail.i
        l::Wline = x[w.head.i]
        id_head::Int = vindex(l, w.head.t)
        id_del::Int = min(id_head, id_head - w.δ)
        # @assert l[id_head] == w.head
        # @assert l[id_tail] == w.tail
        deleteat!(l, range(id_del, length=2))
        return Worm()
    else
        return w
    end
end

function move_worm!(x::Wsheet, w::Worm,
    H::BH_Parameters{znbs}, Q::UpdateConsts)::Worm where {znbs}
    if w.loc == _at_null
        return w
    end
    tail::Element = w.tail
    head::Element = w.head
    D::StateType = randsign()
    δ::StateType = w.δ
    β::f64 = x.β
    loc::WormLocation = w.loc
    i::IndexType = head.i
    li::Wline = x[head.i]
    if (D*δ == -1) # 
        # @assert loc ≠ _at_free "δ≠0 but free???"
        if loc == _at_stop || loc == _at_kink
            # cannot pass either tail or nb kink
            return w
        elseif loc == _at_nbkink || loc == _at_green
            # passby nb interaction && continue to move
            h_id::Int = vindex(li, head.t)
            @reset head.t = nextfloat(head.t, 2D)
            li[h_id] = head
            δ = -δ
            w = Worm(tail, head, δ, loc)
        else # if loc == _at_dummy
            # now pass dummy
            dummy::Element = li[end]
            if δ == +1 && D == -1
                # @assert head == li[1]
                # @assert head.t == nextfloat(0.)
                @reset head.t = prevfloat(β)
                @reset dummy.n_L = head.n_R
                @reset dummy.n_R = head.n_R
                li[end] = head
                push!(li, dummy)
                popfirst!(li)
            elseif δ == -1 && D == +1
                # @assert head == li[end-1]
                # @assert head.t == prevfloat(β)
                @reset head.t = nextfloat(0.)
                @reset dummy.n_L = head.n_L
                @reset dummy.n_R = head.n_L
                pushfirst!(li, head)
                pop!(li)
                li[end] = dummy
            end
            δ = -δ
            w = Worm(tail, head, δ, loc)
            # return w
        end
    end

    # @assert D * δ == +1 || D*δ == 0
    nbs::NTuple{znbs, Int} = get_nbs(H, i)
    t0::f64 = head.t
    head_id::Int = vindex(li, head.t)
    vi_id::Int = mod1(head_id + D, length(li))
    v_near::Element = vi::Element = li[vi_id]
    assocs::NTuple{znbs, Element} = (nb->element_around(x[nb], t0, D, false)).(nbs)
    njs::NTuple{znbs, StateType} = getfield.(assocs, D == StateType(+1) ? :n_L : :n_R)
    λ₊::f64 = λ₋::f64 = Q.Eoff
    Eb::f64 = diagE(H, head.n_L, njs)
    Ef::f64 = diagE(H, head.n_R, njs)
    ΔE::f64 = abs(Ef - Eb)
    if D == StateType(-1)
        Emid = Ef; Ef = Eb ; Eb = Emid
    end
    if D == StateType(+1)
        for e ∈ assocs
            if e.t < v_near.t
                v_near = e
            end
        end
        # v_near = argmin(x->x.t, assocs)
        # v_near = argmin(x->x.t, (vi, v_near))
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
        Δτ = randexp() / λ₊
        t_new = t0 + Δτ
        if (t_new < v_near.t) && metro(δ==0 ? (λ₋ / λ₊) : (1. / λ₊))
            # IF not encounter interaction
            @reset head.t = t_new
            li[head_id] = head
            return Worm(tail, head, i8(0), _at_free)

        elseif (t_new ≥ v_near.t) && metro(δ==0 ? (λ₋ / 1.) : (1. / 1.))
            @reset head.t = prevfloat(v_near.t)
            li[head_id] = head
            
            if v_near.op == I_
                loc = _at_dummy
            elseif v_near.j == IndexType(0) # v_near is the tail
                loc = v_near.i == head.i ? _at_stop : _at_green
            else # v_near is a kink
                loc = v_near.i == head.i ? _at_kink : _at_nbkink
            end
            return Worm(tail, head, i8(-1), loc)
        else
            return w
        end
    else # if D == -1
        for e ∈ assocs
            if mod(e.t, β) > mod(v_near.t, β)
                v_near = e
            end
        end
        # v_near = argmax(x -> mod(x.t, β), assocs)
        # v_near = argmax(x -> mod(x.t, β), (vi, v_near))
        if v_near.t == β # dummy element is at both 0 and β, here we regard it at 0.
            @reset v_near.t = 0.
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
        Δτ = randexp() / λ₊
        t_new = t0 - Δτ

        if (t_new > v_near.t) && metro(δ==0 ? (λ₋ / λ₊) : (1. / λ₊))
            # IF not encounter interaction
            @reset head.t = t_new
            li[head_id] = head
            return Worm(tail, head, i8(0), _at_free)
        elseif (t_new ≤ v_near.t) && metro(δ==0 ? (λ₋ / 1.) : (1. / 1.))
            @reset head.t = nextfloat(v_near.t)
            li[head_id] = head

            if v_near.op == I_
                loc = _at_dummy
            elseif v_near.j == IndexType(0) # v_near is the tail
                loc = v_near.i == head.i ? _at_stop : _at_green
            else # v_near is a kink
                loc = v_near.i == head.i ? _at_kink : _at_nbkink
            end
            return Worm(tail, head, i8(+1), loc)
        else
            return w
        end
    end
end

function insert_kink!(x::Wsheet, w::Worm,
    H::BH_Parameters{znbs}, Q::UpdateConsts)::Worm where {znbs}

    if w.loc ≠ _at_free
        return w
    end # early return for the trivial case

    D::StateType = randsign()
    head::Element = w.head
    i::IndexType = head.i
    nbs::NTuple{znbs, Int} = get_nbs(H, i)
    j::IndexType = rand(nbs) |> IndexType
    B::StateType = StateType(head.op)
    li::Wline = x[i]
    lj::Wline = x[j]
    jqL::Int, jqR::Int = vindex_around(lj, head.t, false)
    nj::StateType = lj[jqR].n_L # @assert nj == lj[jqL].n_R
    nmid::StateType = nj - B * D
    if nmid > H.nmax || nmid < i8(0)
        return w
    end
    Wk::StateType = max(nj, nmid)
    P_acc::f64 = 2 * znbs * Wk * bond_weight(H,i,j) * Q.P_del2ins
    if metro(P_acc)
        # for inb ∈ nbs
        #     if close_to_any(x[inb], head.t)
        #         return w
        #     end
        # end
        iq = vindex(li, head.t) # @assert li[iq] == head

        # first set head as part of K
        li[iq] = Element(head.t, head.i, j, head.n_L, head.n_R, head.op)
        # @set head.j = j

        Kj = Element(head.t, j, i,
            D > 0 ? nj : nmid,
            D > 0 ? nmid : nj,
            OpType(-B)
        )
        head_new = Element(nextfloat(head.t, D), j, IndexType(0), Kj.n_R, Kj.n_L, head.op)
        if D > 0
            insert!(lj, jqR, Kj)
            insert!(lj, jqR + 1, head_new)
        else
            insert!(lj, jqR, head_new)
            insert!(lj, jqR + 1, Kj)
        end
        return Worm(w.tail, head_new, D, _at_kink)
    else
        return w
    end
end

function delete_kink!(x::Wsheet, w::Worm, H::BH_Parameters{znbs}, Q::UpdateConsts)::Worm where {znbs}
    if w.loc ≠ _at_kink
        return w
    end # early return for the trivial case

    head::Element = w.head
    i::IndexType = head.i
    δ::StateType = w.δ
    li::Wline = x[i]
    id_head::Int = vindex(li, head.t)
    id_Ki::Int = id_head - δ
    Ki::Element = li[id_Ki]

    # case pass_interaction / swap Ki and head
    if head.op == Ki.op # case : pass kink with prob 1
        li[id_head] = Element(Ki.t,
            Ki.i, Ki.j, head.n_L, head.n_R, Ki.op
        ) # now id_head is Ki
        head_new = Element(nextfloat(Ki.t, -δ),
            head.i, head.j, Ki.n_L, Ki.n_R, head.op
        )
        li[id_Ki] = head_new
        return Worm(w.tail, head_new, -δ, _at_kink)
    end
    
    # case delete_kink
    j::IndexType = Ki.j
    lj::Wline = x[j]
    id_Kj::Int = vindex(lj, Ki.t)
    Kj::Element = lj[id_Kj]
    Wk::StateType = max(head.n_L, head.n_R)
    P_acc::f64 = inv(2 * znbs * Wk * bond_weight(H, i, j) * Q.P_del2ins)
    if metro(P_acc)
        head_new::Element = Element(Kj.t, Kj.i, IndexType(0), Kj.n_L, Kj.n_R, head.op)
        lj[id_Kj] = head_new
        id_del::Int = min(id_Ki, id_head)
        deleteat!(li, range(id_del, length=2))
        return Worm(w.tail, head_new, StateType(0), _at_free)
    else
        return w
    end
end

function worm_cycle!(x::Wsheet, H::BH_Parameters,
    Q::UpdateConsts, Y::CycleAccumProb,
    )::Int
    cycle_size::Int = hit_dummy::Int = 0
    w0::Worm = w::Worm = insert_worm!(x, Worm(), H, Q)
    while w.loc ≠ _at_null
        dice::f64 = rand()
        if dice < Y.AP_move_worm
            w = move_worm!(x, w, H, Q)
        elseif dice < Y.AP_insert_kink
            w = insert_kink!(x, w, H, Q)
        elseif dice < Y.AP_delete_kink
            w = delete_kink!(x, w, H, Q)
        else#if dice < Y.AP_glue_worm
            w = glue_worm!(x, w, H, Q)
        end
        if w0 ≠ w
            cycle_size += 1
            w0 = w
        end
    end
    return cycle_size # cycle_size ≠ 0 for successful insertion
end
function worm_cycle!(x::Wsheet, H::BH_Parameters,
    Q::UpdateConsts, Y::CycleAccumProb,
    G::GreenFuncBin)::Int
    G.insertion_trial += 1
    cycle_size::Int = 0
    w0::Worm = w::Worm = insert_worm!(x, Worm(), H, Q)
    while w.loc ≠ _at_null
        accum_green!(G, w, H)
        dice::f64 = rand()
        if dice < Y.AP_move_worm
            w = move_worm!(x, w, H, Q)
        elseif dice < Y.AP_insert_kink
            w = insert_kink!(x, w, H, Q)
        elseif dice < Y.AP_delete_kink
            w = delete_kink!(x, w, H, Q)
        else#if dice < Y.AP_glue_worm
            w = glue_worm!(x, w, H, Q)
        end
        if w0 ≠ w
            cycle_size += 1
            w0 = w
        end
    end
    return cycle_size # cycle_size ≠ 0 for successful insertion
end