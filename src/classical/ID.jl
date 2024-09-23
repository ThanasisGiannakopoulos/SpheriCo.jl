
"""
Functions that define the ID of ϕ, Π, Q at t=0, i.e. ϕ(t=0,x).
Needs signature ϕ_ID(x::T) where {T<:Real}
"""

# # classical scalar fields:
# function Φ_ID end # scalar; even
# function Π_ID end # Π := ∂_t Π; even
# function Ψ_ID end # Ψ := ∂_r Ψ; odd

function ID_gauss(r::Real, amp::Float64, width::Float64, rc::Float64)
    f =  amp*exp(-((r - rc) / width)^2 ) + amp*exp(-((r + rc) / width)^2 )
end
function ID_gauss_r(r::Real, amp::Float64, width::Float64, rc::Float64)
    f =  -2.0*((r - rc)/(width^2)) *amp*exp(-((r - rc) / width)^2 )-
    2.0*((r + rc)/(width^2)) *amp*exp(-((r + rc) / width)^2 )
end

# initial data
# classical scalar fields
Φ_ID(r::T, amp::T, width::T, rc::T) where {T<:Float64} = ID_gauss(r, amp, width, rc)
Π_ID(r::T, amp::T, width::T, rc::T) where {T<:Float64} = 0.0
Ψ_ID(r::T, amp::T, width::T, rc::T) where {T<:Float64} = ID_gauss_r(r, amp, width, rc)

# # geometry vars
# function A0_r0 end # metric func; even
# function B_ID end # metric func; even
# function DB_ID end # D_B := (∂_r B)/B ; odd
# function Utld_ID end # U^tilde := (∂_r A)/A - 2 (∂_r B)/B - 2 B λ/A ; odd
# function K_ID end # K := -(1/2α)( (∂_t A)/A + 2(∂_t B)/B) ;  even
# function KB_ID end # K_B := -(1/2α)*(∂_tB)/B ; even
# function λ_ID end # λ := (1-A/B)/r ; odd
# function α_ID end # lapse ; even
# function Dα_ID end # D_α:= (∂_r α)/α ; odd

A0_r0() = 1.0
B_ID(x::T)  where {T<:Float64}  = 1.0
DB_ID(x::T) where {T<:Float64}  = 0.0
Dα_ID(x::T) where {T<:Float64}  = 0.0
K_ID(x::T)  where {T<:Float64}  = 0.0
KB_ID(x::T) where {T<:Float64}  = 0.0
α_ID(x::T)  where {T<:Float64}  = 1.0

"""
Expand that system with two variables Θ and Zr that try to control
violation of Hamiltonian and Momentum constraints, respectively
"""
Θ_ID(x::T)   where {T<:Float64} = 0.0
Zr_ID(x::T)  where {T<:Float64} = 0.0 #anything non-zero is violation of the constraint

"""
We want to be able to build double null coordinates. We expand the system to include
U, V (from ds^2 = dU dV) as well as functions f, g that are used to find
U(t,r) and V(t,r)
"""
f_ID(x::T)  where {T<:Float64}  = 1.0 # = f(0,r) = 1.0
g_ID(x::T)  where {T<:Float64}  = 1.0 # = g(0,r) = 1.0
U0_r0() = 0.0 # U(0,0) = 0.0; dU/dr = -sqrt(A)
V0_r0() = 0.0 # V(0,0) = 0.0; dV/dr = sqrt(A)


# start A0 ID

# classical; real valued fields

# A_ID
# mutates dv

"""
Explore(?):
Maybe there is a better way to solve for the A_ID,
maybe it is better to give the boundary value at r_max and not r=0
"""

function A0_eq_rhs!(dv::Float64,
                    v::Float64,
                    Ψ0::Float64,
                    ri::Float64,
                    p::Param)

    #Λ = p.CC
    oMp2 = p.overMp2 # 1.0 #8.0*π
    
    if ri==0.0
        dv = 0.0
        # the above is true only for A0(r=0)=1
    else
        # always assume CC=0 for classical ID
        dv = v*( (1.0 - v)/ri + 0.5*ri*(Ψ0^2)*oMp2 )
        # the 1st one below works wit CC=|theoretical value|; for backreaction
        #v*( (1.0 - v)/ri + 0.5*ri*(Ψ0^2)*oMp2 - ri*v*Λ*oMp2 )
    end
    
end

"""
RK4 for a scalar (just a number)
to solve for the A_ID; it needs interpolations lib (see later)
"""
#mutate v
function RK4_A0!(r00::Float64, hr::Float64,
                 Ψ00::Float64, Ψ05::Float64, Ψ1::Float64,
                 # for intermediate steps
                 v::Float64,
                 v1::Float64,
                 v2::Float64,
                 v3::Float64,
                 v4::Float64,
                 p::Param)

    v1 = A0_eq_rhs!(v1, v, Ψ00, r00, p)
    v2 = A0_eq_rhs!(v2, v +0.5*hr*v1, Ψ05, r00+0.5hr, p)
    v3 = A0_eq_rhs!(v3, v +0.5*hr*v2, Ψ05, r00+0.5*hr, p)
    v4 = A0_eq_rhs!(v4, v +hr*v3, Ψ1, r00+hr, p)
    v = v + hr*(0.166666667*v1 + 0.333333333*v2 + 0.333333333*v3 + 0.166666667*v4)

    v
end
function RK4_A0(r00::Float64, hr::Float64,
                Ψ00::Float64, Ψ05::Float64, Ψ1::Float64, v::Float64,
                p::Param)

    v1 = 0.0
    v2 = 0.0
    v3 = 0.0
    v4 = 0.0

    v = RK4_A0!(r00, hr, Ψ00, Ψ05, Ψ1, v, v1, v2, v3, v4, p)
    v
end
# end A0_ID

# start V0 ID (for double null coordinates used in post-processing)
# classical; real valued fields
# V_ID
# mutates dv
function V0_eq_rhs!(dv::Float64,
                    v::Float64,
                    A0::Float64,
                    ri::Float64,
                    p::Param)
    dv = sqrt(A0)
end

#mutate v
function RK4_V0!(r00::Float64, hr::Float64,
                 A00::Float64, A05::Float64, A1::Float64,
                 #
                 v::Float64,
                 v1::Float64,
                 v2::Float64,
                 v3::Float64,
                 v4::Float64,
                 p::Param)

    v1 = V0_eq_rhs!(v1, v, A00, r00, p)
    v2 = V0_eq_rhs!(v2, v +0.5*hr*v1, A05, r00+0.5hr, p)
    v3 = V0_eq_rhs!(v3, v +0.5*hr*v2, A05, r00+0.5*hr, p)
    v4 = V0_eq_rhs!(v4, v +hr*v3, A1, r00+hr, p)
    v = v + hr*(0.166666667*v1 + 0.333333333*v2 + 0.333333333*v3 + 0.166666667*v4)

    v
end
function RK4_V0(r00::Float64, hr::Float64,
                A00::Float64, A05::Float64, A1::Float64, v::Float64,
                p::Param)

    v1 = 0.0
    v2 = 0.0
    v3 = 0.0
    v4 = 0.0

    v = RK4_V0!(r00, hr, A00, A05, A1, v, v1, v2, v3, v4, p)
    v
end
# end V0_ID

# start classical_ID
"""
the state vector is:
v_classic = [Φ, Π, Ψ, A, B, DB, Utld, K, KB, λ, α, Dα, Θ, Zr, f, g, U, V]^T    
works for both real and complex valued fields
mutates v_classical
"""
function classical_ID(v_classic::Array, sys::System, p::Param)

    r     = sys.r
    Nr    = length(r)
    r_max = r[end]
    hr    = r[4] - r[3] # r[3] = 0.0

    for i in 3:Nr
        # scalar
        v_classic[i,1]  = Φ_ID(r[i], p.amp, p.width, p.rc)
        v_classic[i,2]  = Π_ID(r[i], p.amp, p.width, p.rc)
        v_classic[i,3]  = Ψ_ID(r[i], p.amp, p.width, p.rc)
        # geometry
        v_classic[i,5]  = B_ID(r[i])
        v_classic[i,6]  = DB_ID(r[i])
        v_classic[i,8]  = K_ID(r[i])
        v_classic[i,9]  = KB_ID(r[i])
        v_classic[i,11] = α_ID(r[i])
        v_classic[i,12] = Dα_ID(r[i])
        v_classic[i,13] = Θ_ID(r[i])
        v_classic[i,14] = Zr_ID(r[i])
        v_classic[i,15] = f_ID(r[i])
        v_classic[i,16] = g_ID(r[i])
    end

    # populate ghost points r[1]=-2h, r[2]=-h
    # even | f[3] = f(r=0)
    # f[1] = f[5] = f(r=2h)
    # f[2] = f[4] = f(r=h)
    v_classic[:,1]  = even_ghosts(v_classic[:,1])  # Φ
    v_classic[:,2]  = even_ghosts(v_classic[:,2])  # Π
    v_classic[:,5]  = even_ghosts(v_classic[:,5])  # B
    v_classic[:,8]  = even_ghosts(v_classic[:,8])  # K
    v_classic[:,9]  = even_ghosts(v_classic[:,9])  # KB
    v_classic[:,11] = even_ghosts(v_classic[:,11]) # α
    v_classic[:,13] = even_ghosts(v_classic[:,13]) # Θ
    # for double null
    v_classic[:,15] = even_ghosts(v_classic[:,15]) # f
    v_classic[:,16] = even_ghosts(v_classic[:,16]) # g

    # odd | f[3] = f(r=0)
    # f[1] = -f[5] = -f(r=2h)
    # f[2] = -f[4] = -f(r=h)
    v_classic[:,3]  = odd_ghosts(v_classic[:,3])  # Ψ
    v_classic[:,6]  = odd_ghosts(v_classic[:,6])  # DB
    v_classic[:,12] = odd_ghosts(v_classic[:,12]) # Dα
    v_classic[:,14] = odd_ghosts(v_classic[:,14]) # Zr

    """
    for A we need to integrate in r
    ∂_r A0 = A0*( (1-A0)/r + r*(Ψ0^2)/(2*Mp^2) + r*A0*Λ/(Mp^2) )
    A0(r=0) = 1.0 (ghosts at i=1,2), r=0 is at i=3
    """

    # with RK4_scalar and 3rd order interpolation for intermediate points for Ψ_ID
    # 3rd order accurate (due to interpolation)
    Ψ0_intp = cubic_spline_interpolation(range(-2*hr, r_max, length=Nr), real(v_classic[:,3]))
    v_classic[3,4] = A0_r0() #1.0 + 0.0im
    @inbounds for i in 3:Nr-1
        v_classic[i+1,4] = RK4_A0(r[i], hr,
                                  Ψ0_intp(r[i]),
                                  Ψ0_intp(r[i] + 0.5*hr),
                                  Ψ0_intp(r[i] + hr),
                                  v_classic[i,4],
                                  p)

    end
    # ghosts are populated after radial integration
    # ghosts; A is even
    v_classic[:,4] = even_ghosts(v_classic[:,4])

    """
    for V we need to integrate in r
    ∂_r V0 = Sqrt[A0]
    V0(r=0) = 0.0 (ghosts at i=1,2), r=0 is at i=3
    since dV/dr is even, V is odd
    """

    # with RK4_scalar and 3rd order interpolation for intermediate points for Ψ_ID
    # 3rd order accurate (due to interpolation)
    A0_intp = cubic_spline_interpolation(range(-2*hr, r_max, length=Nr),
                                         real(v_classic[:,4]))
    v_classic[3,18] = V0_r0() # = 0.0
    @inbounds for i in 3:Nr-1
        v_classic[i+1,18] = RK4_V0(r[i], hr,
                                   A0_intp(r[i]),
                                   A0_intp(r[i] +0.5*hr),
                                   A0_intp(r[i] + hr),
                                   v_classic[i,18],
                                   p)

    end
    # ghosts are populated after radial integration
    # ghosts; V is odd
    v_classic[:,18] = odd_ghosts(v_classic[:,18])
    # for U0 we need -V0 since d U0/dr = -sqrt[A0] = -d V0/dr
    v_classic[:,17] = - v_classic[:,18]

    # for λ and U^tilde we need A
    # λ0 = (1- A0/B0)/r
    v_classic[4:end,10] = (1.0 .- v_classic[4:end,4]./v_classic[4:end,5])./r[4:end]
    # λ0|r=0 ; take 1-A/B->0 as r->0; here A=1, B=1, so λ0|r=0 = 0
    v_classic[3,10] = 0.0
    #ghosts; odd
    v_classic[:,10] = odd_ghosts(v_classic[:,10])

    # Utld0 = (∂_r A0 - 4 λ0)/A0
    v_classic[3:end,7]  = (Dr_FD2(v_classic[:,4], sys)[3:end] .- 4.0.*v_classic[3:end,10])./v_classic[3:end,4]
    # ghosts; Utld is odd
    v_classic[:,7] = odd_ghosts(v_classic[:,7])

    v_classic
end
