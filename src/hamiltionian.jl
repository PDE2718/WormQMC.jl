abstract type BH_Parameters{znbs} end

# Example: all the functions need to be implement for a 2D lattice:
# parameters for 2D Square lattice with only J1
@kwdef struct BH_Square <: BH_Parameters{4}
    nmax::StateType = 1 # nmax must be implemented
    Lx::IndexType = 0 # the size should also be implemented
    Ly::IndexType = 0 # 
    J::f64 = 1. # other parameters
    U::f64 = 0.
    V::f64 = 1.
    μ::f64 = 0.
end
function get_nbs(H::BH_Square, i::Int)::NTuple{4, Int}
    Lx = Int(H.Lx)
    Ly = Int(H.Ly)
    x0, y0 = @inbounds CartesianIndices((Lx, Ly))[i] |> Tuple
    xm, xp = mod1.((x0 - 1, x0 + 1), Lx)
    ym, yp = mod1.((y0 - 1, y0 + 1), Ly)
    idmap = LinearIndices((Lx, Ly))
    @inbounds begin
        return (
            idmap[x0, ym],
            idmap[x0, yp],
            idmap[xm, y0],
            idmap[xp, y0],
        )
    end
end
get_nbs(H::BH_Square, i::Integer) = get_nbs(H, i|>Int)
function diagE(H::BH_Square, ni::StateType, njs::NTuple{4, StateType})::f64
    return 0.5 * H.U * ni * (ni - 1) - H.μ * ni + H.V * ni * sum(njs)
end
bond_weight(H::BH_Square, i::Integer, j::Integer)::f64 = H.J
function get_half_nbs(H::BH_Square, i::Integer)::NTuple{2, Int}
    nbs = get_nbs(H, i)
    return (nbs[1], nbs[3])
end
Wsheet(β::f64, H::BH_Square) = Wsheet(β, zeros(StateType, H.Lx, H.Ly))
function site_diff(H::BH_Square, i, j)::CartesianIndex{4}
    lattice = CartesianIndices((H.Lx, H.Ly))
    xi, yi = Tuple(lattice[i])
    xj, yj = Tuple(lattice[j])
    dx = mod1(xj-xi+1, H.Lx)
    dy = mod1(yj-yi+1, H.Ly)
    return CartesianIndex(dx,dy,1,1)
end

@kwdef struct BH_Pyroch <: BH_Parameters{6}
    nmax::StateType = 1 # nmax must be implemented
    Lx::IndexType = 0 # the size should also be implemented
    Ly::IndexType = 0
    J1::f64 = 1.0 # other parameters
    J2::f64 = 1.0 # other parameters
    U::f64 = 0.0 # on site repulsion, can be zero if nmax==1
    V::f64 = 1.0 # off site repulsion, set to a large value for shaped hard core
    μ::f64 = 0.0 # chemical potential
end
function get_nbs(H::BH_Pyroch, i::Int)::NTuple{6,Int}
    Lx = Int(H.Lx)
    Ly = Int(H.Ly)
    x0, y0, sT = @inbounds CartesianIndices((Lx, Ly, 2))[i] |> Tuple
    xm, xp = mod1.((x0-1, x0+1), Lx)
    ym, yp = mod1.((y0-1, y0+1), Ly)
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
get_nbs(H::BH_Pyroch, i::Integer)=get_nbs(H,i|>Int)
function get_half_nbs(H::BH_Pyroch, i::Integer)::NTuple{3,Int}
    nbs = get_nbs(H, i)
    return (nbs[1], nbs[2], nbs[5])
end
function diagE(H::BH_Pyroch, ni::StateType, njs::NTuple{6,StateType})::f64
    return 0.5 * H.U * ni * (ni - 1) - H.μ * ni + H.V * ni * sum(njs)
end
function bond_weight(H::BH_Pyroch, i::Integer, j::Integer)::f64
    N = H.Lx * H.Ly
    return (i>N)⊻(j>N) ? H.J2 : H.J1
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
