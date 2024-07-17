@kwdef struct BH_Pyroch <: BH_Parameters
    nmax::StateType = 1 # nmax must be implemented
    Lx::IndexType = 0 # the size should also be implemented
    Ly::IndexType = 0
    J1::f64 = 1.0 # other parameters
    J2::f64 = 1.0 # other parameters
    U::f64 = 0.0 # on site repulsion, can be zero if nmax==1
    V::f64 = 1.0 # off site repulsion, set to a large value for shaped hard core
    μ::f64 = 0.0 # chemical potential
end
N_sublatt(::Type{BH_Pyroch}) = 2
N_wldim(::Type{BH_Pyroch})::Int = 3
N_nbs(::Type{BH_Pyroch})::Int = 6
N_hops(::Type{BH_Pyroch})::Int = 6

function get_nbs(H::BH_Pyroch, i::Integer)::NTuple{6,Int}
    Lx = Int(H.Lx)
    Ly = Int(H.Ly)
    x0, y0, sT = @inbounds CartesianIndices((Lx, Ly, 2))[i] |> Tuple
    xm, xp = mod1.((x0 - 1, x0 + 1), Lx)
    ym, yp = mod1.((y0 - 1, y0 + 1), Ly)
    idmap = LinearIndices((Lx, Ly, 2))
    @inbounds begin
        if sT == 1 # i in sublattice A
            return (
                idmap[xm, y0, 2],
                idmap[x0, y0, 2],
                idmap[xm, yp, 2],
                idmap[x0, yp, 2],
                idmap[x0, ym, 1],
                idmap[x0, yp, 1],
            )
        else # if sT == 2 # i in sublattice B
            return (
                idmap[x0, ym, 1],
                idmap[xp, ym, 1],
                idmap[x0, y0, 1],
                idmap[xp, y0, 1],
                idmap[xm, y0, 2],
                idmap[xp, y0, 2],
            )
        end
    end
end
get_hops(H::BH_Pyroch, i::Integer)::NTuple{6,Int} = get_nbs(H, i)

function diagE(H::BH_Pyroch, ni::StateType, njs::NTuple{6,StateType})::f64
    return 0.5 * H.U * ni * (ni - 1) - H.μ * ni + H.V * ni * sum(njs)
end
function bond_weight(H::BH_Pyroch, i::Integer, j::Integer)::f64
    N = H.Lx * H.Ly
    return (i > N) ⊻ (j > N) ? H.J2 : H.J1
end
Wsheet(β::f64, H::BH_Pyroch) = Wsheet(β, zeros(StateType, H.Lx, H.Ly, 2))
function site_diff(H::BH_Pyroch, i, j)::CartesianIndex{4}
    lattice = CartesianIndices((H.Lx, H.Ly, 2))
    xi, yi, si = Tuple(lattice[i])
    xj, yj, sj = Tuple(lattice[j])
    dx = mod1(xj - xi + 1, H.Lx)
    dy = mod1(yj - yi + 1, H.Ly)
    return CartesianIndex(dx, dy, si, sj)
end

function winding_number(x::Wsheet{3}, H::BH_Pyroch)::NTuple{2,Int}
    Wx::Int = Wy::Int = 0
    N::Int = H.Lx * H.Ly
    for i ∈ eachindex(x)
        for e::Element ∈ x[i]::Wline
            if e.op == b_
                @assert e.i == i
                j = Int(e.j)
                hops = get_hops(H, i)
                # dir = findfirst(==(e.j), get_hops(H, i))
                sid = i ≤ N ? 1 : 2
                if j == hops[1]
                    Wx += (-1)
                    Wy += (-1)
                elseif j == hops[2]
                    Wx += (+1)
                    Wy += (-1)
                elseif j == hops[3]
                    Wx += (-1)
                    Wy += (+1)
                elseif j == hops[4]
                    Wx += (+1)
                    Wy += (+1)
                elseif j == hops[5]
                    if sid == 1
                        Wy -= 2
                    else
                        Wx -= 2
                    end
                elseif j == hops[6]
                    if sid == 1
                        Wy += 2
                    else
                        Wx += 2
                    end
                else
                    error("dir illegal")
                end
            end
        end
    end
    @assert Wx % (2H.Lx) == Wy % (2H.Ly) == 0
    WxWy = (Wx ÷ (2H.Lx), Wy ÷ (2H.Ly))
    return WxWy
end
