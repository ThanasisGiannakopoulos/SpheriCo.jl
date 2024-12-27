# classical and quantum rhs in one function
# MULTIPLE DISPATCH: to treat the non-reg. and reg. quantum versions

# NON-REGULARIZED
#mutates vf_classic and vf_quantum
function F!(t::Float64,
            vf_classic::Array{Float64},
            v_classic::Array{Float64}, 
            vf_quantum::Array{ComplexF64, 4},
            v_quantum::Array{ComplexF64, 4}, 
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
    # set incoming characteristic variables equal to zero
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
    Kminus3KB       = K - 3.0*KB
    Kminus3KB_SBP21 = Dr_SBP2(Kminus3KB, sys, 1)
    Kminus3KB_r     = Dr_FD2(Kminus3KB, sys)
    Mom = @. -0.5 * DB * Kminus3KB + 0.5 * jA * oMp2 + (KB_r) + (Kminus3KB_r) - 1. * (Kminus3KB_SBP21)

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

    """quantum"""
    @threads for k in 1:p.kmax
        # index for uq, ψq, and πq in the state vector (for each l and k)
        uqi =  1
        ψqi =  2
        πqi =  3
        uq = v_quantum[:, 1, Int(k), uqi] # even
        ψq = v_quantum[:, 1, Int(k), ψqi] # odd 
        πq = v_quantum[:, 1, Int(k), πqi] # even
        # populate ghosts
        # even
        uq = even_ghosts(uq)
        πq = even_ghosts(πq)
        # odd
        ψq = odd_ghosts(ψq)

        # for the eoms and BC
        # FD2
        uq_r = Dr_FD2(uq, sys)
        πq_r = Dr_FD2(πq, sys)
        ψq_r = Dr_FD2(ψq, sys)

        # SBP2 with L = p.l + 1
        ψq_SBP2lplus1   = Dr_SBP2(ψq, sys, 1)
        # KO diss
        uq_diss   = -p.sigma* KO_FD2(uq, sys)
        πq_diss   = -p.sigma* KO_FD2(πq, sys)
        ψq_diss   = -p.sigma* KO_FD2(ψq, sys)

        # eoms
        # ∂_t uq = ; even
        vf_quantum[:, 1, Int(k), uqi] = @. oA^0.5 * oB * α * πq + uq_diss
        # ∂_t ψq = ; odd
        vf_quantum[:, 1, Int(k), ψqi] = @. oA^0.5 *
            oB * α * ((-0.5 * DA - 1. * DB + Dα) * πq + (πq_r)) +
            ψq_diss
        # ∂_t πq = ; even
        vf_quantum[:, 1, Int(k), πqi] = @.  0.5 * B *
            oA^0.5 * α * (-1. * DA * ψq + 2. * (DB + Dα) *
            ψq - 2. * (ψq_r) + 4. * (ψq_SBP2lplus1)) + πq_diss
    end # end k loop

    # for l>0
    @threads for l in 1:p.lmax
        for k in 1:p.kmax
            # index for uq, ψq, and πq in the state vector (for each l and k)
            uqi =  1
            ψqi =  2
            πqi =  3
            uq = v_quantum[:, Int(l+1), Int(k), uqi] # even
            ψq = v_quantum[:, Int(l+1), Int(k), ψqi] # odd 
            πq = v_quantum[:, Int(l+1), Int(k), πqi] # even
            # populate ghosts
            # even
            uq = even_ghosts(uq)
            πq = even_ghosts(πq)
            # odd
            ψq = odd_ghosts(ψq)

            # for the eoms and BC
            # FD2
            uq_r = Dr_FD2(uq, sys)
            πq_r = Dr_FD2(πq, sys)
            ψq_r = Dr_FD2(ψq, sys)

            # SBP2 with L = p.l
            λ_SBP2l  = Dr_SBP2(λ, sys, Int(l))
            DA_SBP2l = Dr_SBP2(DA, sys, Int(l))
            DB_SBP2l = Dr_SBP2(DB, sys, Int(l))
            Dα_SBP2l = Dr_SBP2(Dα, sys, Int(l))
            # SBP2 with L = p.l + 1
            ψq_SBP2lplus1   = Dr_SBP2(ψq, sys, Int(l) + 1)
            # KO diss
            uq_diss   = -p.sigma* KO_FD2(uq, sys)
            πq_diss   = -p.sigma* KO_FD2(πq, sys)
            ψq_diss   = -p.sigma* KO_FD2(ψq, sys)

            # eoms
            # ∂_t uq = ; even
            vf_quantum[:, Int(l+1), Int(k), uqi] = @. oA^0.5 * oB * α * πq + uq_diss
            # ∂_t ψq = ; odd
            vf_quantum[:, Int(l+1), Int(k), ψqi] = @. oA^0.5 *
                oB * α * ((-0.5 * DA - 1. * DB + Dα) * πq + (πq_r)) +
                ψq_diss
            # ∂_t πq = ; even
            vf_quantum[:, Int(l+1), Int(k), πqi] = @.  0.5 * B *
                oA^0.5 * α * (-1. * DA * ψq + 2. * (DB + Dα) *
                ψq + uq * ((DA_r) - 1. * (DA_SBP2l) + 2. *
                (DB_SBP2l) + 2. * (Dα_SBP2l) - 2. * ((DB_r) +
                (Dα_r) + (1. + l) * (λ_r)) + 2. * (1. + l) *
                (λ_SBP2l)) - 2. * (ψq_r) + 4. * (ψq_SBP2lplus1)) + πq_diss
        end # end k loop
    end # end l loop

    nothing
end

# REGULARIZED
# mutates vf_classic and vf_quantum
# runs with or without backreaction using "if"
function F!(t::Float64,
            vf_classic::Array{Float64},
            v_classic::Array{Float64}, 
            vf_quantum::Array{ComplexF64, 5},
            v_quantum::Array{ComplexF64, 5}, 
            sys::System,
            p::Param)

    r  = sys.r
    Nr = length(r)
    hbar = p.hbar
    c = 1.0

    # make filter for outer boundary (smooth step function based on tanh)
    filter = zeros(Nr)
    #filter2 = zeros(Nr)
    if p.backreaction==true
        if p.infalling_rmax==false
            # filter the r rescaling because the free BC at rmax make
            # the simulation unstable in the outer part of the domain;
            # this behaviour seems similar also with outgoing BC for
            # quantum modes we filter 2 times; ones the r domain
            # before rescaling and once later in the bilinears the
            # filter is based on a tanh and steepness controls how
            # steep it is at r_cut
            steepness = p.steepness
            r_cut = p.r_cut
            filter[3:end] = 0.5*ones(Nr-2) + 0.5*tanh.(-r[3:end].^steepness .+
                (r_cut^steepness).*ones(Nr-2))
            # steepness can be e.g. 1.5 and it does not work well for r<0, so copy it after
            filter[1] = filter[5]
            filter[2] = filter[4]
            # # filter 2
            # steepness2 = 0.9
            # r_cut2 = r_cut+5
            # filter2[3:end] = 0.5*ones(Nr-2) + 0.5*tanh.(-r[3:end].^steepness2 .+
            #     (r_cut2^steepness2).*ones(Nr-2))
            # filter2[1] = filter2[5]
            # filter2[2] = filter2[4]
            # # multiply filer
            # filter=filter.*filter2
        else
            steepness = p.steepness
            r_cut = p.r_cut - t
            filter[3:end] = 0.5*ones(Nr-2) + 0.5*tanh.(-r[3:end].^steepness .+
                (r_cut^steepness).*ones(Nr-2))
            # steepness can be e.g. 1.5 and it does not work well for r<0, so copy it after
            filter[1] = filter[5]
            filter[2] = filter[4]
        end # if infalling_ramx
    end # if backreaction

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
    oα = 1.0 ./ α
    # construct DA:
    DA = @. Utld + 2.0*DB + 4.0*B*λ*oA
    # for quantum modes
    DA_r   = Dr_FD2(DA, sys)
    DB_r   = Dr_FD2(DB, sys)
    Dα_r   = Dr_FD2(Dα, sys)
    # define 1/Mp^2 = 8 π G, where G is Netwon's constant and is set to 1 here
    oMp2 = p.overMp2 # 1.0 # 8*π

    Λ = p.CC*ones(length(r))

    """ 
    quantum: build the rhs for quantum modes, as well as
    the bilinears for the backreaction (if there is). 
    The loops in l are split in l=0 and l>=1,
    due to the structure of the bilinears.
    The loops in m are also expanded.
    """

    # initiate buffers to multithread bilinears' calculation
    if p.backreaction==true
        bilin = zeros(ComplexF64, (length(r), 5))
        # buffers for the sums in multithreading with
        # format (r, number_of_quantum_bilins, number_of_threads)
        buffers = zeros(ComplexF64, (length(r), 5, Threads.nthreads()))
    end

    # for l=0
    @threads for k in 1:p.kmax
        id = Threads.threadid()

        # uq_mi: even
        uq_m0 = v_quantum[:, 1, Int(k), 1, 1]
        uq_m1 = v_quantum[:, 1, Int(k), 2, 1]
        uq_m2 = v_quantum[:, 1, Int(k), 3, 1]
        uq_m5 = v_quantum[:, 1, Int(k), 4, 1]
        # ψq_mi: odd
        ψq_m0 = v_quantum[:, 1, Int(k), 1, 2]
        ψq_m1 = v_quantum[:, 1, Int(k), 2, 2]
        ψq_m2 = v_quantum[:, 1, Int(k), 3, 2]
        ψq_m5 = v_quantum[:, 1, Int(k), 4, 2]
        # πq_mi: even
        πq_m0 = v_quantum[:, 1, Int(k), 1, 3]
        πq_m1 = v_quantum[:, 1, Int(k), 2, 3]
        πq_m2 = v_quantum[:, 1, Int(k), 3, 3]
        πq_m5 = v_quantum[:, 1, Int(k), 4, 3]

        # populate ghosts
        # even
        uq_m0 = even_ghosts(uq_m0)
        uq_m1 = even_ghosts(uq_m1)
        uq_m2 = even_ghosts(uq_m2)
        uq_m5 = even_ghosts(uq_m5)
        πq_m0 = even_ghosts(πq_m0)
        πq_m1 = even_ghosts(πq_m1)
        πq_m2 = even_ghosts(πq_m2)
        πq_m5 = even_ghosts(πq_m5)
        # odd
        ψq_m0 = odd_ghosts(ψq_m0)
        ψq_m1 = odd_ghosts(ψq_m1)
        ψq_m2 = odd_ghosts(ψq_m2)
        ψq_m5 = odd_ghosts(ψq_m5)

        # for the eoms and BC
        # FD2
        # uq
        uq_m0_r = Dr_FD2(uq_m0, sys)
        uq_m1_r = Dr_FD2(uq_m1, sys)
        uq_m2_r = Dr_FD2(uq_m2, sys)
        uq_m5_r = Dr_FD2(uq_m5, sys)
        # πq
        πq_m0_r = Dr_FD2(πq_m0, sys)
        πq_m1_r = Dr_FD2(πq_m1, sys)
        πq_m2_r = Dr_FD2(πq_m2, sys)
        πq_m5_r = Dr_FD2(πq_m5, sys)
        # ψq
        ψq_m0_r = Dr_FD2(ψq_m0, sys)
        ψq_m1_r = Dr_FD2(ψq_m1, sys)
        ψq_m2_r = Dr_FD2(ψq_m2, sys)
        ψq_m5_r = Dr_FD2(ψq_m5, sys)

        # SBP2 with L = p.l + 1 = 1 here
        ψq_m0_SBP2lplus1 = Dr_SBP2(ψq_m0, sys, 1)
        ψq_m1_SBP2lplus1 = Dr_SBP2(ψq_m1, sys, 1)
        ψq_m2_SBP2lplus1 = Dr_SBP2(ψq_m2, sys, 1)
        ψq_m5_SBP2lplus1 = Dr_SBP2(ψq_m5, sys, 1)
        # KO diss
        # uq
        uq_m0_diss = -p.sigma* KO_FD2(uq_m0, sys)
        uq_m1_diss = -p.sigma* KO_FD2(uq_m1, sys)
        uq_m2_diss = -p.sigma* KO_FD2(uq_m2, sys)
        uq_m5_diss = -p.sigma* KO_FD2(uq_m5, sys)
        # πq
        πq_m0_diss = -p.sigma* KO_FD2(πq_m0, sys)
        πq_m1_diss = -p.sigma* KO_FD2(πq_m1, sys)
        πq_m2_diss = -p.sigma* KO_FD2(πq_m2, sys)
        πq_m5_diss = -p.sigma* KO_FD2(πq_m5, sys)
        # ψq
        ψq_m0_diss = -p.sigma* KO_FD2(ψq_m0, sys)
        ψq_m1_diss = -p.sigma* KO_FD2(ψq_m1, sys)
        ψq_m2_diss = -p.sigma* KO_FD2(ψq_m2, sys)
        ψq_m5_diss = -p.sigma* KO_FD2(ψq_m5, sys)

        # the mass value
        # mi = p.mlist[m]
        # eoms
        # ∂_t uq = ; even
        vf_quantum[:, 1, Int(k), 1, 1] = @. oA^0.5 * oB * α * πq_m0 + uq_m0_diss
        vf_quantum[:, 1, Int(k), 2, 1] = @. oA^0.5 * oB * α * πq_m1 + uq_m1_diss
        vf_quantum[:, 1, Int(k), 3, 1] = @. oA^0.5 * oB * α * πq_m2 + uq_m2_diss
        vf_quantum[:, 1, Int(k), 4, 1] = @. oA^0.5 * oB * α * πq_m5 + uq_m5_diss
        # ∂_t ψq = ; odd
        vf_quantum[:, 1, Int(k), 1, 2] = @. oA^0.5 *
            oB * α * ((-0.5 * DA - 1. * DB + Dα) * πq_m0 + (πq_m0_r)) +
            ψq_m0_diss
        vf_quantum[:, 1, Int(k), 2, 2] = @. oA^0.5 *
            oB * α * ((-0.5 * DA - 1. * DB + Dα) * πq_m1 + (πq_m1_r)) +
            ψq_m1_diss
        vf_quantum[:, 1, Int(k), 3, 2] = @. oA^0.5 *
            oB * α * ((-0.5 * DA - 1. * DB + Dα) * πq_m2 + (πq_m2_r)) +
            ψq_m2_diss
        vf_quantum[:, 1, Int(k), 4, 2] = @. oA^0.5 *
            oB * α * ((-0.5 * DA - 1. * DB + Dα) * πq_m5 + (πq_m5_r)) +
            ψq_m5_diss
        # ∂_t πq = ; even
        vf_quantum[:, 1, Int(k), 1, 3] = @.  0.5 * B *
            oA^0.5 * α * (-1. * DA * ψq_m0 + 2. * (DB + Dα) *
            ψq_m0 - 2. * (ψq_m0_r) + 4. * (ψq_m0_SBP2lplus1)) -
            α*B*A^(0.5)*p.mlist[1]^(2.0)*uq_m0 + πq_m0_diss
        vf_quantum[:, 1, Int(k), 2, 3] = @.  0.5 * B *
            oA^0.5 * α * (-1. * DA * ψq_m1 + 2. * (DB + Dα) *
            ψq_m1 - 2. * (ψq_m1_r) + 4. * (ψq_m1_SBP2lplus1)) -
            α*B*A^(0.5)*p.mlist[2]^(2.0)*uq_m1 + πq_m1_diss
        vf_quantum[:, 1, Int(k), 3, 3] = @.  0.5 * B *
            oA^0.5 * α * (-1. * DA * ψq_m2 + 2. * (DB + Dα) *
            ψq_m2 - 2. * (ψq_m2_r) + 4. * (ψq_m2_SBP2lplus1)) -
            α*B*A^(0.5)*p.mlist[3]^(2.0)*uq_m2 + πq_m2_diss
        vf_quantum[:, 1, Int(k), 4, 3] = @.  0.5 * B *
            oA^0.5 * α * (-1. * DA * ψq_m5 + 2. * (DB + Dα) *
            ψq_m5 - 2. * (ψq_m5_r) + 4. * (ψq_m5_SBP2lplus1)) -
            α*B*A^(0.5)*p.mlist[4]^(2.0)*uq_m5 + πq_m5_diss

        # radiative boundary condition at rmax
        # ∂_t f + vout*∂_r f + f/rmax = 0
        # uq
        vf_quantum[end, 1, Int(k), 1, 1] = - vout_rmax*uq_m0_r[end] + uq_m0_diss[end] - uq_m0[end]/sys.r[end]
        vf_quantum[end, 1, Int(k), 2, 1] = - vout_rmax*uq_m1_r[end] + uq_m1_diss[end] - uq_m1[end]/sys.r[end]
        vf_quantum[end, 1, Int(k), 3, 1] = - vout_rmax*uq_m2_r[end] + uq_m2_diss[end] - uq_m2[end]/sys.r[end]
        vf_quantum[end, 1, Int(k), 4, 1] = - vout_rmax*uq_m5_r[end] + uq_m5_diss[end] - uq_m5[end]/sys.r[end]
        # ψq
        vf_quantum[end, 1, Int(k), 1, 2] = - vout_rmax*ψq_m0_r[end] + ψq_m0_diss[end] - ψq_m0[end]/sys.r[end]
        vf_quantum[end, 1, Int(k), 2, 2] = - vout_rmax*ψq_m1_r[end] + ψq_m1_diss[end] - ψq_m1[end]/sys.r[end]
        vf_quantum[end, 1, Int(k), 3, 2] = - vout_rmax*ψq_m2_r[end] + ψq_m2_diss[end] - ψq_m2[end]/sys.r[end]
        vf_quantum[end, 1, Int(k), 4, 2] = - vout_rmax*ψq_m5_r[end] + ψq_m5_diss[end] - ψq_m5[end]/sys.r[end]
        # πq
        vf_quantum[end, 1, Int(k), 1, 3] = - vout_rmax*πq_m0_r[end] + πq_m0_diss[end] - πq_m0[end]/sys.r[end]
        vf_quantum[end, 1, Int(k), 2, 3] = - vout_rmax*πq_m1_r[end] + πq_m1_diss[end] - πq_m1[end]/sys.r[end]
        vf_quantum[end, 1, Int(k), 3, 3] = - vout_rmax*πq_m2_r[end] + πq_m2_diss[end] - πq_m2[end]/sys.r[end]
        vf_quantum[end, 1, Int(k), 4, 3] = - vout_rmax*πq_m5_r[end] + πq_m5_diss[end] - πq_m5[end]/sys.r[end]

        if p.backreaction==true
            # utld_kl = r^l*uq; for l=0, utld_kl = uq
            utld_kl_m0 = uq_m0
            utld_kl_m1 = uq_m1
            utld_kl_m2 = uq_m2
            utld_kl_m5 = uq_m5
            # ∂r_utld_kl = l*r^(l-1)*uq + r^l*ψq = ψq for l=0
            dr_utld_kl_m0 = ψq_m0
            dr_utld_kl_m1 = ψq_m1
            dr_utld_kl_m2 = ψq_m2
            dr_utld_kl_m5 = ψq_m5
            # ∂t_utld_kl = r^l*α*πq/(B*A^0.5)
            dt_utld_kl_m0 = @. α* πq_m0* oB* oA^0.5
            dt_utld_kl_m1 = @. α* πq_m1* oB* oA^0.5
            dt_utld_kl_m2 = @. α* πq_m2* oB* oA^0.5
            dt_utld_kl_m5 = @. α* πq_m5* oB* oA^0.5

            # bilin[:,1]
            # (4*π)/(hbar*c^2)* [μ^2*<Φ^2>]_quantum = Sum_{klm} dk*μ^2[m]*weight[m]*|utld_klm|^2;
            # where μ=m is the quantum scalar mass, p.mlist = [0*mPV, 1*mPV, sqrt(3)*mPV, 2*mPV]
            # and weight[m] = [1,-2,2,-1]
            buffers[:,1,id] .+= @. p.dk*(p.mlist[1]^2 *utld_kl_m0*conj(utld_kl_m0) -
                p.mlist[2]^2 *2.0*utld_kl_m1*conj(utld_kl_m1) +
                p.mlist[3]^2 *2.0*utld_kl_m2*conj(utld_kl_m2) -
                p.mlist[4]^2 *utld_kl_m5*conj(utld_kl_m5))
            # bilin[:,2]
            # (4*π*α^2)/(hbar*c^2*A*B^2)* <Π^2>_quantum =
            #  Sum_{klm} dk*(2*l+1)*weight[m]*|∂_t utld_klm|^2
            # l=0 here and weight[m] = [1,-2,2,-1]
            buffers[:,2,id] .+= @. p.dk*(dt_utld_kl_m0*conj(dt_utld_kl_m0) -
                2.0*dt_utld_kl_m1*conj(dt_utld_kl_m1) +
                2.0*dt_utld_kl_m2*conj(dt_utld_kl_m2) -
                dt_utld_kl_m5*conj(dt_utld_kl_m5))
            # bilin[:,3]
            # (4*π)/(hbar*c^2)* <Ψ^2>_quantum =
            #  Sum_{klm} dk*(2l+1)*weight[m]*|∂_r utld_klm|^2
            # l=0 here and weight[m] = [1,-2,2,-1]
            buffers[:,3,id] .+= @. p.dk*(dr_utld_kl_m0*conj(dr_utld_kl_m0) -
                2.0*dr_utld_kl_m1*conj(dr_utld_kl_m1) +
                2.0*dr_utld_kl_m2*conj(dr_utld_kl_m2) -
                dr_utld_kl_m5*conj(dr_utld_kl_m5))
            # bilin[:,4]
            # (4*π*α)/(hbar*c^2*A^(1/2)*B)* <Π*Ψ>_quantum =
            #  Sum_{klm} dk*0.5*(2*l+1)*weight[m]*[(∂_r utld_klm)*conj(∂_t utld_klm) +
            #    conj(∂_r utld_klm)(∂_t utld_klm)]
            # l=0 here and weight[m] = [1,-2,2,-1]
            buffers[:,4,id] .+= @. p.dk*(0.5*dr_utld_kl_m0*conj(dt_utld_kl_m0) +
                0.5*dt_utld_kl_m0*conj(dr_utld_kl_m0) -
                dr_utld_kl_m1*conj(dt_utld_kl_m1) -
                dt_utld_kl_m1*conj(dr_utld_kl_m1) +
                dr_utld_kl_m2*conj(dt_utld_kl_m2) +
                dt_utld_kl_m2*conj(dr_utld_kl_m2) -
                0.5*dr_utld_kl_m5*conj(dt_utld_kl_m5) -
                0.5*dt_utld_kl_m5*conj(dr_utld_kl_m5))
            # bilin[:,5]
            # (4*π)/(hbar*c^2)* [(1/r^2)*<∂_θ Φ^2>]_quantum =
            #  Sum_{klm} 0.5*weight[m]*l*(l+1)*(2l+1)*|r^(l-1)*u_klm|^2;
            # l=0 here so this contribution is 0
            buffers[:,5,id] = buffers[:,5,id] .+ zeros(ComplexF64, length(bilin[:,5]))
        end # if backreaction
    end # end loops in k

    # sum for l>0
    @threads for l in 1:p.lmax
        for k in 1:p.kmax
            id = Threads.threadid()

            uq_m0 = v_quantum[:, Int(l+1), Int(k), 1, 1]
            uq_m1 = v_quantum[:, Int(l+1), Int(k), 2, 1]
            uq_m2 = v_quantum[:, Int(l+1), Int(k), 3, 1]
            uq_m5 = v_quantum[:, Int(l+1), Int(k), 4, 1]

            ψq_m0 = v_quantum[:, Int(l+1), Int(k), 1, 2]
            ψq_m1 = v_quantum[:, Int(l+1), Int(k), 2, 2]
            ψq_m2 = v_quantum[:, Int(l+1), Int(k), 3, 2]
            ψq_m5 = v_quantum[:, Int(l+1), Int(k), 4, 2]

            πq_m0 = v_quantum[:, Int(l+1), Int(k), 1, 3]
            πq_m1 = v_quantum[:, Int(l+1), Int(k), 2, 3]
            πq_m2 = v_quantum[:, Int(l+1), Int(k), 3, 3]
            πq_m5 = v_quantum[:, Int(l+1), Int(k), 4, 3]

            # populate ghosts
            # even
            uq_m0 = even_ghosts(uq_m0)
            uq_m1 = even_ghosts(uq_m1)
            uq_m2 = even_ghosts(uq_m2)
            uq_m5 = even_ghosts(uq_m5)
            πq_m0 = even_ghosts(πq_m0)
            πq_m1 = even_ghosts(πq_m1)
            πq_m2 = even_ghosts(πq_m2)
            πq_m5 = even_ghosts(πq_m5)
            # odd
            ψq_m0 = odd_ghosts(ψq_m0)
            ψq_m1 = odd_ghosts(ψq_m1)
            ψq_m2 = odd_ghosts(ψq_m2)
            ψq_m5 = odd_ghosts(ψq_m5)

            # for the eoms and BC
            # FD2
            uq_m0_r = Dr_FD2(uq_m0, sys)
            uq_m1_r = Dr_FD2(uq_m1, sys)
            uq_m2_r = Dr_FD2(uq_m2, sys)
            uq_m5_r = Dr_FD2(uq_m5, sys)
            πq_m0_r = Dr_FD2(πq_m0, sys)
            πq_m1_r = Dr_FD2(πq_m1, sys)
            πq_m2_r = Dr_FD2(πq_m2, sys)
            πq_m5_r = Dr_FD2(πq_m5, sys)
            ψq_m0_r = Dr_FD2(ψq_m0, sys)
            ψq_m1_r = Dr_FD2(ψq_m1, sys)
            ψq_m2_r = Dr_FD2(ψq_m2, sys)
            ψq_m5_r = Dr_FD2(ψq_m5, sys)

            # SBP2 with L = p.l
            λ_SBP2l  = Dr_SBP2(λ, sys, Int(l))
            DA_SBP2l = Dr_SBP2(DA, sys, Int(l))
            DB_SBP2l = Dr_SBP2(DB, sys, Int(l))
            Dα_SBP2l = Dr_SBP2(Dα, sys, Int(l))
            # SBP2 with L = p.l + 1 = 1 here
            ψq_m0_SBP2lplus1 = Dr_SBP2(ψq_m0, sys, Int(l+1))
            ψq_m1_SBP2lplus1 = Dr_SBP2(ψq_m1, sys, Int(l+1))
            ψq_m2_SBP2lplus1 = Dr_SBP2(ψq_m2, sys, Int(l+1))
            ψq_m5_SBP2lplus1 = Dr_SBP2(ψq_m5, sys, Int(l+1))
            # KO diss
            uq_m0_diss = -p.sigma* KO_FD2(uq_m0, sys)
            uq_m1_diss = -p.sigma* KO_FD2(uq_m1, sys)
            uq_m2_diss = -p.sigma* KO_FD2(uq_m2, sys)
            uq_m5_diss = -p.sigma* KO_FD2(uq_m5, sys)
            πq_m0_diss = -p.sigma* KO_FD2(πq_m0, sys)
            πq_m1_diss = -p.sigma* KO_FD2(πq_m1, sys)
            πq_m2_diss = -p.sigma* KO_FD2(πq_m2, sys)
            πq_m5_diss = -p.sigma* KO_FD2(πq_m5, sys)
            ψq_m0_diss = -p.sigma* KO_FD2(ψq_m0, sys)
            ψq_m1_diss = -p.sigma* KO_FD2(ψq_m1, sys)
            ψq_m2_diss = -p.sigma* KO_FD2(ψq_m2, sys)
            ψq_m5_diss = -p.sigma* KO_FD2(ψq_m5, sys)

            # the mass value
            # mi = p.mlist[m]
            # eoms
            # ∂_t uq = ; even
            vf_quantum[:, Int(l+1), Int(k), 1, 1] = @. oA^0.5 * oB * α * πq_m0 + uq_m0_diss
            vf_quantum[:, Int(l+1), Int(k), 2, 1] = @. oA^0.5 * oB * α * πq_m1 + uq_m1_diss
            vf_quantum[:, Int(l+1), Int(k), 3, 1] = @. oA^0.5 * oB * α * πq_m2 + uq_m2_diss
            vf_quantum[:, Int(l+1), Int(k), 4, 1] = @. oA^0.5 * oB * α * πq_m5 + uq_m5_diss
            # ∂_t ψq = ; odd
            vf_quantum[:, Int(l+1), Int(k), 1, 2] = @. oA^0.5 *
                oB * α * ((-0.5 * DA - 1. * DB + Dα) * πq_m0 + (πq_m0_r)) +
                ψq_m0_diss
            vf_quantum[:, Int(l+1), Int(k), 2, 2] = @. oA^0.5 *
                oB * α * ((-0.5 * DA - 1. * DB + Dα) * πq_m1 + (πq_m1_r)) +
                ψq_m1_diss
            vf_quantum[:, Int(l+1), Int(k), 3, 2] = @. oA^0.5 *
                oB * α * ((-0.5 * DA - 1. * DB + Dα) * πq_m2 + (πq_m2_r)) +
                ψq_m2_diss
            vf_quantum[:, Int(l+1), Int(k), 4, 2] = @. oA^0.5 *
                oB * α * ((-0.5 * DA - 1. * DB + Dα) * πq_m5 + (πq_m5_r)) +
                ψq_m5_diss
            # ∂_t πq = ; even
            vf_quantum[:, Int(l+1), Int(k), 1, 3] = @.  0.5 * B *
                oA^0.5 * α * (-1. * DA * ψq_m0 + 2. * (DB + Dα) *
                ψq_m0 + uq_m0 * ((DA_r) - 1. * (DA_SBP2l) + 2. *
                (DB_SBP2l) + 2. * (Dα_SBP2l) - 2. * ((DB_r) +
                (Dα_r) + (1. + l) * (λ_r)) + 2. * (1. + l) *
                (λ_SBP2l)) - 2. * (ψq_m0_r) + 4. * (ψq_m0_SBP2lplus1)) -
                α*B*A^(0.5)*p.mlist[1]^(2.0)*uq_m0 + πq_m0_diss
            vf_quantum[:, Int(l+1), Int(k), 2, 3] = @.  0.5 * B *
                oA^0.5 * α * (-1. * DA * ψq_m1 + 2. * (DB + Dα) *
                ψq_m1 + uq_m1 * ((DA_r) - 1. * (DA_SBP2l) + 2. *
                (DB_SBP2l) + 2. * (Dα_SBP2l) - 2. * ((DB_r) +
                (Dα_r) + (1. + l) * (λ_r)) + 2. * (1. + l) *
                (λ_SBP2l)) - 2. * (ψq_m1_r) + 4. * (ψq_m1_SBP2lplus1)) -
                α*B*A^(0.5)*p.mlist[2]^(2.0)*uq_m1 + πq_m1_diss
            vf_quantum[:, Int(l+1), Int(k), 3, 3] = @.  0.5 * B *
                oA^0.5 * α * (-1. * DA * ψq_m2 + 2. * (DB + Dα) *
                ψq_m2 + uq_m2 * ((DA_r) - 1. * (DA_SBP2l) + 2. *
                (DB_SBP2l) + 2. * (Dα_SBP2l) - 2. * ((DB_r) +
                (Dα_r) + (1. + l) * (λ_r)) + 2. * (1. + l) *
                (λ_SBP2l)) - 2. * (ψq_m2_r) + 4. * (ψq_m2_SBP2lplus1)) -
                α*B*A^(0.5)*p.mlist[3]^(2.0)*uq_m2 + πq_m2_diss
            vf_quantum[:, Int(l+1), Int(k), 4, 3] = @.  0.5 * B *
                oA^0.5 * α * (-1. * DA * ψq_m5 + 2. * (DB + Dα) *
                ψq_m5 + uq_m5 * ((DA_r) - 1. * (DA_SBP2l) + 2. *
                (DB_SBP2l) + 2. * (Dα_SBP2l) - 2. * ((DB_r) +
                (Dα_r) + (1. + l) * (λ_r)) + 2. * (1. + l) *
                (λ_SBP2l)) - 2. * (ψq_m5_r) + 4. * (ψq_m5_SBP2lplus1)) -
                α*B*A^(0.5)*p.mlist[4]^(2.0)*uq_m5 + πq_m5_diss

            # radiative boundary condition at rmax
            # ∂_t f + vout*∂_r f + f/rmax = 0
            # uq
            vf_quantum[end, Int(l+1), Int(k), 1, 1] = - vout_rmax*uq_m0_r[end] + uq_m0_diss[end] - uq_m0[end]/sys.r[end]
            vf_quantum[end, Int(l+1), Int(k), 2, 1] = - vout_rmax*uq_m1_r[end] + uq_m1_diss[end] - uq_m1[end]/sys.r[end]
            vf_quantum[end, Int(l+1), Int(k), 3, 1] = - vout_rmax*uq_m2_r[end] + uq_m2_diss[end] - uq_m2[end]/sys.r[end]
            vf_quantum[end, Int(l+1), Int(k), 4, 1] = - vout_rmax*uq_m5_r[end] + uq_m5_diss[end] - uq_m5[end]/sys.r[end]
            # ψq
            vf_quantum[end, Int(l+1), Int(k), 1, 2] = - vout_rmax*ψq_m0_r[end] + ψq_m0_diss[end] - ψq_m0[end]/sys.r[end]
            vf_quantum[end, Int(l+1), Int(k), 2, 2] = - vout_rmax*ψq_m1_r[end] + ψq_m1_diss[end] - ψq_m1[end]/sys.r[end]
            vf_quantum[end, Int(l+1), Int(k), 3, 2] = - vout_rmax*ψq_m2_r[end] + ψq_m2_diss[end] - ψq_m2[end]/sys.r[end]
            vf_quantum[end, Int(l+1), Int(k), 4, 2] = - vout_rmax*ψq_m5_r[end] + ψq_m5_diss[end] - ψq_m5[end]/sys.r[end]
            # πq
            vf_quantum[end, Int(l+1), Int(k), 1, 3] = - vout_rmax*πq_m0_r[end] + πq_m0_diss[end] - πq_m0[end]/sys.r[end]
            vf_quantum[end, Int(l+1), Int(k), 2, 3] = - vout_rmax*πq_m1_r[end] + πq_m1_diss[end] - πq_m1[end]/sys.r[end]
            vf_quantum[end, Int(l+1), Int(k), 3, 3] = - vout_rmax*πq_m2_r[end] + πq_m2_diss[end] - πq_m2[end]/sys.r[end]
            vf_quantum[end, Int(l+1), Int(k), 4, 3] = - vout_rmax*πq_m5_r[end] + πq_m5_diss[end] - πq_m5[end]/sys.r[end]

            if p.backreaction==true

                # apply the filter to r as well
                rfilt = r.*filter
                # rlminus1_u_kl_mi := r^(l-1)*u_klm
                rlminus1_u_kl_m0 = @. rfilt^(l-1.0)*uq_m0
                rlminus1_u_kl_m1 = @. rfilt^(l-1.0)*uq_m1
                rlminus1_u_kl_m2 = @. rfilt^(l-1.0)*uq_m2
                rlminus1_u_kl_m5 = @. rfilt^(l-1.0)*uq_m5

                # ∂t_utld_kl = r^l*α*πq/(B*A^0.5)
                dt_utld_kl_m0 = @. rfilt^l* α* πq_m0/(B* A^0.5)
                dt_utld_kl_m1 = @. rfilt^l* α* πq_m1/(B* A^0.5)
                dt_utld_kl_m2 = @. rfilt^l* α* πq_m2/(B* A^0.5)
                dt_utld_kl_m5 = @. rfilt^l* α* πq_m5/(B* A^0.5)

                # ∂r_utld_kl = l*r^(l-1)*uq + r^l*ψq
                dr_utld_kl_m0 = @. (l*rfilt^(l-1)*uq_m0 + ψq_m0*rfilt^l)
                dr_utld_kl_m1 = @. (l*rfilt^(l-1)*uq_m1 + ψq_m1*rfilt^l)
                dr_utld_kl_m2 = @. (l*rfilt^(l-1)*uq_m2 + ψq_m2*rfilt^l)
                dr_utld_kl_m5 = @. (l*rfilt^(l-1)*uq_m5 + ψq_m5*rfilt^l)

                # utld_kl = r^l*uq
                utld_kl_m0 = @. rfilt^l*uq_m0
                utld_kl_m1 = @. rfilt^l*uq_m1
                utld_kl_m2 = @. rfilt^l*uq_m2
                utld_kl_m5 = @. rfilt^l*uq_m5

                # bilin[:,1]
                # (4*π)/(hbar*c^2)*[μ^2*<Φ^2>]_quantum=Sum_{klm} dk*μ^2[m]*weight[m]*|utld_klm|^2;
                # μ=m is the quantum scalar mass, p.mlist = [0*mPV, 1*mPV, sqrt(3)*mPV, 2*mPV]
                # weight[m] = [1,-2,2,-1]
                buffers[:,1,id] .+= @. p.dk*(2.0*l+1.0)*(
                    p.mlist[1]^2 *utld_kl_m0*conj(utld_kl_m0) -
                    p.mlist[2]^2 *2.0*utld_kl_m1*conj(utld_kl_m1) +
                    p.mlist[3]^2 *2.0*utld_kl_m2*conj(utld_kl_m2) -
                    p.mlist[4]^2 *utld_kl_m5*conj(utld_kl_m5))
                # bilin[:,2]
                # (4*π*α^2)/(hbar*c^2*A*B^2)* <Π^2>_quantum =
                #  Sum_{klm} dk*(2*l+1)*weight[m]*|∂_t utld_klm|^2;
                # weight[m] = [1,-2,2,-1]
                buffers[:,2,id] .+= @. p.dk*(2.0*l+1.0)*(
                    dt_utld_kl_m0*conj(dt_utld_kl_m0) -
                    2.0*dt_utld_kl_m1*conj(dt_utld_kl_m1) +
                    2.0*dt_utld_kl_m2*conj(dt_utld_kl_m2) -
                    dt_utld_kl_m5*conj(dt_utld_kl_m5))
                # bilin[:,3]
                # (4*π)/(hbar*c^2)* <Ψ^2>_quantum = Sum_{klm} dk*(2l+1)*weight[m]*|∂_r utld_klm|^2;
                # weight[m] = [1,-2,2,-1]
                buffers[:,3,id] .+= @. p.dk*(2.0*l+1.0)*(
                    dr_utld_kl_m0*conj(dr_utld_kl_m0) -
                    2.0*dr_utld_kl_m1*conj(dr_utld_kl_m1) +
                    2.0*dr_utld_kl_m2*conj(dr_utld_kl_m2) -
                    dr_utld_kl_m5*conj(dr_utld_kl_m5))
                # bilin[:,4]
                # (4*π*α)/(hbar*c^2*A^(1/2)*B)* <Π*Ψ>_quantum =
                #  Sum_{klm} dk*0.5*(2*l+1)*weight[m]*[(∂_r utld_klm)*conj(∂_t utld_klm) +
                #    conj(∂_r utld_klm)(∂_t utld_klm)];
                # weight[m] = [1,-2,2,-1]
                buffers[:,4,id] .+= @. p.dk*(2.0*l+1.0)*(
                    0.5*dr_utld_kl_m0*conj(dt_utld_kl_m0) +
                    0.5*dt_utld_kl_m0*conj(dr_utld_kl_m0) -
                    dr_utld_kl_m1*conj(dt_utld_kl_m1) -
                    dt_utld_kl_m1*conj(dr_utld_kl_m1) +
                    dr_utld_kl_m2*conj(dt_utld_kl_m2) +
                    dt_utld_kl_m2*conj(dr_utld_kl_m2) -
                    0.5*dr_utld_kl_m5*conj(dt_utld_kl_m5) -
                    0.5*dt_utld_kl_m5*conj(dr_utld_kl_m5))
                # bilin[:,5]
                # (4*π)/(hbar*c^2)* [(1/r^2)*<∂_θ Φ^2>]_quantum =
                #  Sum_{klm} dk*0.5*weight[m]*l*(l+1)*(2l+1)*|r^(l-1)*u_klm|^2;
                # weight[m] = [1,-2,2,-1]
                buffers[:,5,id] .+= @. p.dk*l*(l+1.0)*(2.0*l+1.0)*(
                    0.5*rlminus1_u_kl_m0*conj(rlminus1_u_kl_m0) -
                        rlminus1_u_kl_m1*conj(rlminus1_u_kl_m1) +
                        rlminus1_u_kl_m2*conj(rlminus1_u_kl_m2) -
                        0.5*rlminus1_u_kl_m5*conj(rlminus1_u_kl_m5))
            end # if backreaction

        end # end loops in k
    end # end loops in l

    # classic stress-energy tensor components:
    jA = -(oA.^0.5).*oB.*Π.*Ψ
    SA = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)
    SB = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .- Ψ.^2.0)
    ρ  = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)

    """
    bln[:,1] = (4*π)/(hbar*c^2)* [μ^2*<Φ^2>]_quantum = Sum_{klm} μ^2[m]*weight[m]*|utld_klm|^2;
    bln[:,2] = (4*π*α^2)/(hbar*c^2*A*B^2)* <Π^2>_quantum = Sum_{klm} (2*l+1)*weight[m]*
                                                                     |∂_t utld_klm|^2 
    bln[:,3] = (4*π)/(hbar*c^2)* <Ψ^2>_quantum = Sum_{klm} (2l+1)*weight[m]*|∂_r utld_klm|^2
    bln[:,4] = (4*π*α)/(hbar*c^2*A^(1/2)*B)* <Π*Ψ>_quantum = Sum_{klm} 0.5*(2*l+1)*weight[m]*
                            [(∂_r utld_klm)*conj(∂_t utld_klm) + conj(∂_r utld_klm)(∂_t utld_klm)]
    bln[:,5] =  (4*π)/(hbar*c^2)* [(1/r^2)*<∂_θ Φ^2>]_quantum = Sum_{klm} 0.5*weight[m]*l*(l+1)*
                                                                (2l+1)*|r^(l-1)*u_klm|^2
    """
    # if there is backreaction, add its contribution
    if p.backreaction==true
        bilin[:,1] .= sum(buffers[:,1,:], dims=2)
        bilin[:,2] .= sum(buffers[:,2,:], dims=2)
        bilin[:,3] .= sum(buffers[:,3,:], dims=2)
        bilin[:,4] .= sum(buffers[:,4,:], dims=2)
        bilin[:,5] .= sum(buffers[:,5,:], dims=2)

        # populate ghosts
        bilin[:,1]  = even_ghosts(bilin[:,1])
        bilin[:,2]  = even_ghosts(bilin[:,2])
        bilin[:,3]  = even_ghosts(bilin[:,3])
        bilin[:,4]  = odd_ghosts(bilin[:,4])
        bilin[:,5]  = even_ghosts(bilin[:,5])

        ρ_quantum = @. (hbar*c^2/(4.0*π))*(
            (0.5* oα^2)*bilin[:,2] + (0.5*oA)*bilin[:,3] + oB*bilin[:,5] + 0.5*bilin[:,1] )
        jA_quantum = @. -(hbar*c^2/(4.0*π))*oα*bilin[:,4]
        SA_quantum = @. (hbar*c^2/(4.0*π))*( (0.5*oα^2)*bilin[:,2] +
            (0.5*oA)*bilin[:,3] - oB*bilin[:,5] - 0.5*bilin[:,1] )
        SB_quantum = @. (hbar*c^2/(4.0*π))*( (0.5*oα^2)*bilin[:,2] -
            (0.5*oA)*bilin[:,3] - 0.5*bilin[:,1] )

        # check if imaginary part is above machine precision
        if maximum(abs.(imag.(ρ_quantum).*filter)) > 10^(-15)
            println("max|Im(ρ_quantum)| = ", maximum(abs.(imag.(ρ_quantum))))
        end
        if maximum(abs.(imag.(jA_quantum).*filter)) > 10^(-15)
            println("max|Im(jA_quantum)| = ", maximum(abs.(imag.(jA_quantum))))
        end
        if maximum(abs.(imag.(SA_quantum).*filter)) > 10^(-15)
            println("max|Im(SA_quantum)| = ", maximum(abs.(imag.(SA_quantum))))
        end
        if maximum(abs.(imag.(SB_quantum).*filter)) > 10^(-15)
            println("max|Im(SB_quantum)| = ", maximum(abs.(imag.(SB_quantum))))
        end

        """
        we filter out the backreaction ~in the part of the r-domain
        that is causally connected to rmax (otherwise we are ruining the asymptotics)
        """

        # the filter is step function (based on a tanh) that is 1 from r=0
        # until ~ r=rmax-tmax (causally connecter)
        # and the turns zero smoothly
        ρ  = ρ  .+ real.(ρ_quantum).*filter
        jA = jA .+ real.(jA_quantum).*filter
        SA = SA .+ real.(SA_quantum).*filter
        SB = SB .+ real.(SB_quantum).*filter
    end # if backreaction

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

    # otherwise unstable
    if p.backreaction==true
        Ham = 0.0*Ham
        Mom = 0.0*Mom
    end

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
    #∂_t K = ; even | apply filter in the matter fields-> in backreaction case mismatch of bilins and CC lead to bad BC at rmax for classical eoms
    vf_classic[:,8] = @. 0.5 * oA^2. * oB * α * (B * (A * Dα * (-2. * Dα + Utld) + 4. * B * Dα * λ + A^2 * (12. * KB^2 - 8. * KB * K + 2. * K^2 + 6. * p.damping * Θ * p.κ1 * (1. + p.κ2) + oMp2 * (SA + 2. * SB - 2. * Λ * filter + ρ))) + 2. * A * (4. * A * p.damping * (Zr_r) + B * ((Dα_r) - 2. * ((Dα_SBP21) + p.damping * (Zr_r))) - 4. * A * p.damping * (Zr_SBP21))) + K_diss
    #println("source K[3] = ", (oMp2 * (SA .+ 2. * SB .- 2. * Λ .+ ρ))[3])
    #∂_t KB = ; even | apply filter in the matter fields-> in backreaction case mismatch of bilins and CC lead to bad BC at rmax for classical eoms
    vf_classic[:,9] = @. 0.25 * oA^2. * oB * α * (B * (A * DB * (-2. * Dα + Utld) + 4. * B * DB * λ - 2. * A^2 * (-2. * KB * K + 2. * p.damping * Θ * p.κ1 * (1. + p.κ2) + oMp2 * (-1. * SA + 2.0 * Λ * filter + ρ))) + 2. * A * B * (DB_r) - 4. * A * B * (DB_SBP21) + 4. * A * B * (Dα_r) - 4. * A * B * (Dα_SBP21) - 2. * A * B * (Utld_r) + 2. * A * B * (Utld_SBP21) - 8. * A^2 * p.damping * (Zr_r) + 8. * A^2 * p.damping * (Zr_SBP21) + 4. * A * B * (λ_r) - 8. * B^2 * (λ_r) - 4. * (A - 2. * B) * B * (λ_SBP21)) + KB_diss
    #println("source KB[3] = ", (oMp2 * (-1. * SA .+ 2.0*Λ + ρ))[3] )
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
