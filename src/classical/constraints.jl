# start constraints

# classical case; real valued

# v = [Φ, Π, Ψ, A, B, DB, Utld, K, KB, λ, α, Dα, Θ, Zr]^T
"""
TODO: check the function below; it is not really used somewhere in the evolution
only in postprocessing (def below without System and Param)
"""
function constraints(v::Array, sys::System, p::Param)

    # scalar
    Φ    = v[:,1] # even
    Π    = v[:,2] # even
    Ψ    = v[:,3] # odd
    # geometry
    A    = v[:,4] # even
    B    = v[:,5] # even
    DB   = v[:,6] # odd
    Utld = v[:,7] # odd
    K    = v[:,8] # even
    KB   = v[:,9] # even
    λ    = v[:,10] # odd
    α    = v[:,11] # even
    Dα   = v[:,12] # odd

    # populate ghost points r[1]=-2h, r[2]=-h
    # even
    Φ  = even_ghosts(Φ)
    Π  = even_ghosts(Π)
    A  = even_ghosts(A)
    B  = even_ghosts(B)
    K  = even_ghosts(K)
    KB = even_ghosts(KB)
    α  = even_ghosts(α)

    # odd
    Ψ    = odd_ghosts(Ψ)
    DB   = odd_ghosts(DB)
    Utld = odd_ghosts(Utld)
    λ    = odd_ghosts(λ)
    Dα   = odd_ghosts(Dα)

    # FD2 centered r derivatives
    λ_r    = Dr_FD2(λ, sys)
    DB_r   = Dr_FD2(DB, sys)
    Utld_r = Dr_FD2(Utld, sys)
    K_r    = Dr_FD2(K, sys)
    KB_r   = Dr_FD2(KB, sys)

    # define 1/A and 1/B metric function; it is faster to multiply than devide
    oA = 1.0 ./ A
    oB = 1.0 ./ B
    # 1/Mp2 convention
    oMp2 = p.overMp2
    # stress-energy tensor components:
    jA = -(oA.^0.5).*oB.*Π.*Ψ
    SA = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)
    SB = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .- Ψ.^2.0)
    ρ  = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)

    # Hamiltonian constraint 
    Ham = @. DB_r + (λ + DB - Utld - 4.0*λ*B/A)/sys.r -
    DB*(0.25*DB + 0.5*Utld + 2.0*λ*B/A) -
    A*KB*(2.0*K - 3.0*KB) + A*ρ*oMp2
    # for r=0, l 'Hospital for 1/r term
    Ham[3] = DB_r[3] + (λ_r[3] + DB_r[3] - Utld_r[3] - 4.0*λ_r[3]*B[3]/A[3]) -
    DB[3]*(0.25*DB[3] + 0.5*Utld[3] + 2.0*λ[3]*B[3]/A[3]) -
    A[3]*KB[3]*(2.0*K[3] - 3.0*KB[3]) + A[3]*ρ[3]*oMp2

    # Momentum constraint
    Mom = @. KB_r - (1.0/sys.r + DB*0.5)*(K - 3.0*KB) + 0.5*jA*oMp2
    # for r=0, use that KA = KB = 0?
    Mom[3] = KB_r[3] + 0.5*jA[3]*oMp2

    Ham, Mom
end

# without the need of System and Param; useful for postprocessing
function constraints(v::Array, dr::Float64, oMp2::Float64)

    # scalar
    Φ    = v[:,1] # even
    Π    = v[:,2] # even
    Ψ    = v[:,3] # odd
    # geometry
    A    = v[:,4] # even
    B    = v[:,5] # even
    DB   = v[:,6] # odd
    Utld = v[:,7] # odd
    K    = v[:,8] # even
    KB   = v[:,9] # even
    λ    = v[:,10] # odd
    α    = v[:,11] # even
    Dα   = v[:,12] # odd

    # populate ghost points r[1]=-2h, r[2]=-h
    # even
    Φ  = even_ghosts(Φ)
    Π  = even_ghosts(Π)
    A  = even_ghosts(A)
    B  = even_ghosts(B)
    K  = even_ghosts(K)
    KB = even_ghosts(KB)
    α  = even_ghosts(α)

    # odd
    Ψ    = odd_ghosts(Ψ)
    DB   = odd_ghosts(DB)
    Utld = odd_ghosts(Utld)
    λ    = odd_ghosts(λ)
    Dα   = odd_ghosts(Dα)

    # FD2 centered r derivatives
    λ_r    = Dr_FD2(λ, dr)
    DB_r   = Dr_FD2(DB, dr)
    Utld_r = Dr_FD2(Utld, dr)
    K_r    = Dr_FD2(K, dr)
    KB_r   = Dr_FD2(KB, dr)

    λ_SBP21    = Dr_SBP2(λ, dr, 1)
    DB_SBP21   = Dr_SBP2(DB, dr, 1)
    Utld_SBP21 = Dr_SBP2(Utld, dr, 1)

    # define 1/A and 1/B metric function; it is faster to multiply than devide
    oA = 1.0 ./ A
    oB = 1.0 ./ B
    # 1/Mp2 convention
    #oMp2 = p.overMp2
    # stress-energy tensor components:
    jA = -(oA.^0.5).*oB.*Π.*Ψ
    SA = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)
    SB = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .- Ψ.^2.0)
    ρ  = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)

    Ham = @. -0.25 * DB * (DB + 2. * Utld) - 2. * B * DB * oA * λ + A * (3. * KB^2 - 2. * KB * K + oMp2 * ρ) + (DB_SBP21) + (Utld_r) - 1. * (Utld_SBP21) - 1. * (A - 4. * B) * oA * ((λ_r) - 1. * (λ_SBP21))

    # Momentum constraint
    # K - 3.0*KB in theory is zero at r=0, but is even, not odd; still let's try SBP21
    Kminus3KB       = K - 3.0*KB
    Kminus3KB_SBP21 = Dr_SBP2(Kminus3KB, dr, 1)
    Kminus3KB_r     = Dr_FD2(Kminus3KB, dr)
    Mom = @. -0.5 * DB * Kminus3KB + 0.5 * jA * oMp2 + (KB_r) + (Kminus3KB_r) - 1. * (Kminus3KB_SBP21)

    Ham, Mom
end

# end constraints
