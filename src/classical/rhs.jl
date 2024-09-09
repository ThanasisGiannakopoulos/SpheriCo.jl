# start quantum rhs

#mutates vf_classic and vf_quantum
function F!(t::Float64,
            vf_classic::Array,
            v_classic::Array, 
            sys::System,
            p::Param)

    r  = sys.r
    Nr = length(r)
    
    """classical"""
    # scalar
    Φ    = v_classic[:,1] # even
    Π    = v_classic[:,2] # even
    Ψ    = v_classic[:,3] # odd
    # geometry
    A    = v_classic[:,4] # even
    B    = v_classic[:,5] # even
    DB   = v_classic[:,6] # odd
    Utld = v_classic[:,7] # odd
    K    = v_classic[:,8] # even
    KB   = v_classic[:,9] # even
    λ    = v_classic[:,10] # odd
    α    = v_classic[:,11] # even
    Dα   = v_classic[:,12] # odd
    # constraint violation
    # Hamiltonian
    Θ    = v_classic[:,13] #  even
    # Momentum (Zr is Z_subscript_r; index down)
    Zr   = v_classic[:,14] #  odd
    # for double null
    f    = v_classic[:,15] #  even
    g    = v_classic[:,16] #  even
    U    = v_classic[:,17] #  odd
    V    = v_classic[:,18] #  odd

    # populate ghost points r[1]=-2h, r[2]=-h
    # even | f[3] = f(r=0)
    # f[1] = f[5] = f(r=2h)
    # f[2] = f[4] = f(r=h)
    Φ  = even_ghosts(Φ)
    Π  = even_ghosts(Π)
    A  = even_ghosts(A)
    B  = even_ghosts(B)
    K  = even_ghosts(K)
    KB = even_ghosts(KB)
    α  = even_ghosts(α)
    Θ  = even_ghosts(Θ)
    f  = even_ghosts(f)
    g  = even_ghosts(g)

    # odd | f[3] = f(r=0)
    # f[1] = -f[5] = -f(r=2h)
    # f[2] = -f[4] = -f(r=h)
    Ψ    = odd_ghosts(Ψ)
    DB   = odd_ghosts(DB)
    Utld = odd_ghosts(Utld)
    λ    = odd_ghosts(λ)
    Dα   = odd_ghosts(Dα)
    Zr   = odd_ghosts(Zr)
    U    = odd_ghosts(U)
    V    = odd_ghosts(V)

    # BC for geometry based on linearized characteristic vars
    KB[end] = (sqrt(A[end])*Θ[end] + Zr[end])/(2.0*sqrt(A[end]))
    DB[end] = sqrt(A[end])*Θ[end] + Zr[end]
    K[end]  = -(1.0/(sqrt(2.0*A[end])*(-2.0 + α[end])))*
    (4.0*sqrt(2.0*A[end])*Θ[end] + 2.0*Dα[end]*sqrt(α[end]) -
     4.0*sqrt(A[end]*α[end])*Θ[end] + 4.0*Zr[end]*sqrt(α[end]) -
     2.0*sqrt(2.0)*Zr[end]*α[end] - Dα[end]*α[end]^(1.5))

    # speed of propagation for sommerfeld diff than 1, for scalars
    vout_rmax = α[end]/sqrt(A[end])

    # FD2 centered r derivatives
    Φ_r    = Dr_FD2(Φ, sys)
    Π_r    = Dr_FD2(Π, sys)
    Ψ_r    = Dr_FD2(Ψ, sys)
    A_r    = Dr_FD2(A, sys)
    B_r    = Dr_FD2(B, sys)
    α_r    = Dr_FD2(α, sys)
    K_r    = Dr_FD2(K, sys)
    KB_r   = Dr_FD2(KB, sys)
    λ_r    = Dr_FD2(λ, sys)
    DB_r   = Dr_FD2(DB, sys)
    Dα_r   = Dr_FD2(Dα, sys)
    Utld_r = Dr_FD2(Utld, sys)
    # constraint violation control vars
    Θ_r    = Dr_FD2(Θ, sys)
    Zr_r   = Dr_FD2(Zr, sys)
    # for double null
    f_r    = Dr_FD2(f, sys)
    g_r    = Dr_FD2(g, sys)

    """ SBP2 operator that approximates to 2nd order accuracy f' +
        L*f/r, for integer L>=1 and odd f """
    Ψ_SBP21    = Dr_SBP2(Ψ, sys, 1)
    Zr_SBP21   = Dr_SBP2(Zr, sys, 1)
    λ_SBP21    = Dr_SBP2(λ, sys, 1)
    DB_SBP21   = Dr_SBP2(DB, sys, 1)
    Dα_SBP21   = Dr_SBP2(Dα, sys, 1)
    Utld_SBP21 = Dr_SBP2(Utld, sys, 1)

    # Kreiss-Oliger dissipation compatible to FD2 (i.e. 3rd order KO diss)
    # scalar
    Φ_diss = -p.sigma* KO_FD2(Φ, sys)
    Π_diss = -p.sigma* KO_FD2(Π, sys)
    Ψ_diss = -p.sigma* KO_FD2(Ψ, sys)
    # geometry
    A_diss    = -p.sigma* KO_FD2(A, sys)
    B_diss    = -p.sigma* KO_FD2(B, sys)
    DB_diss   = -p.sigma* KO_FD2(DB, sys)
    Utld_diss = -p.sigma* KO_FD2(Utld, sys)
    K_diss    = -p.sigma* KO_FD2(K, sys)
    KB_diss   = -p.sigma* KO_FD2(KB, sys)
    λ_diss    = -p.sigma* KO_FD2(λ, sys)
    α_diss    = -p.sigma* KO_FD2(α, sys)
    Dα_diss   = -p.sigma* KO_FD2(Dα, sys)
    # constriant violation
    Θ_diss    = -p.sigma* KO_FD2(Θ, sys)
    Zr_diss   = -p.sigma* KO_FD2(Zr, sys)
    # for double null
    f_diss    = -p.sigma* KO_FD2(f, sys)
    g_diss    = -p.sigma* KO_FD2(g, sys)
    U_diss    = -p.sigma* KO_FD2(U, sys)
    V_diss    = -p.sigma* KO_FD2(V, sys)

    # define 1/A and 1/B metric function; it is faster to multiply than devide
    oA = 1.0 ./ A
    oB = 1.0 ./ B
    # construct DA:
    DA = @. Utld + 2.0*DB + 4.0*B*λ*oA
    # for quantum modes
    DA_r   = Dr_FD2(DA, sys)
    DB_r   = Dr_FD2(DB, sys)
    Dα_r   = Dr_FD2(Dα, sys)
    # define 1/Mp^2 = 8 π G, where G is Netwon's constant and is set to 1 here
    oMp2 = p.overMp2 # 1.0 # 8*π

    Λ = p.CC

    # stress-energy tensor components:
    jA = -(oA.^0.5).*oB.*Π.*Ψ
    SA = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)
    SB = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .- Ψ.^2.0)
    ρ  = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)

    # Hamiltonian constraint
    """ The SBP2 operatores replace:
    1st line: λ/sys.r
    2nd line: DB/sys.r
    3rd line: Utld/sys.r
    4th line: oA*B*4.0*λ/sys.r
    """
    Ham = @. -0.25 * DB * (DB + 2. * Utld) - 2. * B * DB * oA * λ + A * (3. * KB^2 - 2. * KB * K + oMp2 * ρ) + (DB_SBP21) + (Utld_r) - 1. * (Utld_SBP21) - 1. * (A - 4. * B) * oA * ((λ_r) - 1. * (λ_SBP21))

    # Momentum constraint
    # K - 3.0*KB in theory is zero at r=0, but is even, not odd; still let's try SBP21
    Kminus3KB       = K - 3.0*KB
    Kminus3KB_SBP21 = Dr_SBP2(Kminus3KB, sys, 1)
    Kminus3KB_r     = Dr_FD2(Kminus3KB, sys)
    Mom = @. -0.5 * DB * Kminus3KB + 0.5 * jA * oMp2 + (KB_r) + (Kminus3KB_r) - 1. * (Kminus3KB_SBP21)
    # the following fails badly
    #Mom = @. -0.5 * DB * Kminus3KB + 0.5 * jA * oMp2 + (KB_r) +  - (Kminus3KB)/r
    #Mom[3] = -0.5 * DB[3] * Kminus3KB[3] + 0.5 * jA[3] * oMp2 + (KB_r[3]) +  - (Kminus3KB_r[3])

    # EOMS:
    # scalar field in First order in Time and 1st order in Space (FT1S) form
    #∂_t Φ = ; even
    vf_classic[:,1] = @. oA^0.5 * oB * α * Π + Φ_diss
    #∂_t Π = ; even
    vf_classic[:,2] = @. -0.5 * B * oA^0.5 * α * ((DA - 2. * (DB + Dα)) * Ψ + 2. * (Ψ_r) - 4. * (Ψ_SBP21)) + Π_diss
    #∂_t Ψ = ; odd
    vf_classic[:,3] = @. 0.5 * oA^0.5 * oB * α * (-1. * (DA + 2. * DB - 2. * Dα) * Π + 2. * (Π_r)) + Ψ_diss

    # radiative boundary condition at rmax
    # ∂_t f + vout*∂_r f + f/rmax = 0
    vf_classic[end,1] = - vout_rmax*Φ_r[end] + Φ_diss[end] - Φ[end]/sys.r[end]
    vf_classic[end,2] = - vout_rmax*Π_r[end] + Π_diss[end] - Π[end]/sys.r[end]
    vf_classic[end,3] = - vout_rmax*Ψ_r[end] + Ψ_diss[end] - Ψ[end]/sys.r[end]

    # geometry
    #∂_t A = ; even
    vf_classic[:,4] = @. -2. * A * (-2. * KB + K) * α + A_diss
    #∂_t B = ; even
    vf_classic[:,5] = @. -2. * B * KB * α + B_diss
    #∂_t DB = ; odd
    vf_classic[:,6] = @. -2. * α * (Dα * KB + (KB_r)) + DB_diss
    #∂_t Utld = ; odd
    vf_classic[:,7] = @. -2. * oA * α * (A * (6. * DB * KB - 4. * Dα * KB - 2. * DB * K + Dα * K + 2. * jA * oMp2) + 4. * B * (-3. * KB + K) * λ + A * (K_r)) + Utld_diss
    #∂_t K = ; even
    vf_classic[:,8] = @. 0.5 * oA^2. * oB * α * (B * (A * Dα * (-2. * Dα + Utld) + 4. * B * Dα * λ + A^2 * (12. * KB^2 - 8. * KB * K + 2. * K^2 + 6. * p.damping * Θ * p.κ1 * (1. + p.κ2) + oMp2 * (SA + 2. * SB - 2. * Λ + ρ))) + 2. * A * (4. * A * p.damping * (Zr_r) + B * ((Dα_r) - 2. * ((Dα_SBP21) + p.damping * (Zr_r))) - 4. * A * p.damping * (Zr_SBP21))) + K_diss
    #∂_t KB = ; even
    vf_classic[:,9] = @. 0.25 * oA^2. * oB * α * (B * (A * DB * (-2. * Dα + Utld) + 4. * B * DB * λ - 2. * A^2 * (-2. * KB * K + 2. * p.damping * Θ * p.κ1 * (1. + p.κ2) + oMp2 * (-1. * SA + 2.0*Λ + ρ))) + 2. * A * B * (DB_r) - 4. * A * B * (DB_SBP21) + 4. * A * B * (Dα_r) - 4. * A * B * (Dα_SBP21) - 2. * A * B * (Utld_r) + 2. * A * B * (Utld_SBP21) - 8. * A^2 * p.damping * (Zr_r) + 8. * A^2 * p.damping * (Zr_SBP21) + 4. * A * B * (λ_r) - 8. * B^2 * (λ_r) - 4. * (A - 2. * B) * B * (λ_SBP21)) + KB_diss
    #∂_t λ = ; odd
    vf_classic[:,10] = @. A * oB * α * (3. * DB * KB - 1. * DB * K + jA * oMp2 + 2. * (KB_r)) + λ_diss
    #∂_t α = ; even
    vf_classic[:,11] = @. -2. * α * (K - 2. * p.damping * Θ) + α_diss 
    #∂_t Dα = ; odd
    vf_classic[:,12] = @. -2. * (K_r) + 4. * p.damping * (Θ_r) + Dα_diss 
    # damping vars
    # ∂_t Θ = ; even
    vf_classic[:,13] = @. -1. * oA * oB * (B * α * Ham + A * B * Θ * p.κ1 * (2. + p.κ2) + (2. * A - 1. * B) * (Zr_r) - 2. * A * (Zr_SBP21)) + Θ_diss
    # ∂_t Zr = ; odd
    vf_classic[:,14] = @. (-2. * α * Mom - 1. * Zr * p.κ1 + (Θ_r)) + Zr_diss
    # for double null
    # ∂_t f = ; even
    vf_classic[:,15] = @. f*α*(-2*KB + K) - (α*f_r + f*α_r)/(sqrt(A)) + f_diss
    # ∂_t g = ; even
    vf_classic[:,16] = @. -g*α*(-2*KB + K) - (α*g_r + g*α_r)/(sqrt(A)) + g_diss
    # ∂_t U = ; odd
    vf_classic[:,17] = @. f*α + U_diss
    # ∂_t V = ; odd
    vf_classic[:,18] = @. g*α + V_diss

    nothing
end
