function find_AH(sys::System, B::Vector, A::Vector, KB::Vector)

    r = sys.r
    B_r = Dr_FD2(B, sys)
    Θ = (2.0./r .+ B_r./B)./(A.^0.5) - 2.0.*KB

    test = sign(Θ[end])
    j = length(r)
    while sign(Θ[j])==test && j>=2 # is the 1st ghost point on the left
    test = sign(Θ[j])
    j = j-1
    end

    # we take the midpoint where the sign of the expansion changes. If this is <0, then there is no AH
    rAH = (r[j] + r[j+1])/2
    rAH
end
# multiple dispatch
# for tools
function find_AH(r::Vector, B::Vector, A::Vector, KB::Vector)

    hr = r[4] - r[3]
    B_r = Dr_FD2(B, hr)
    Θ = (2.0./r .+ B_r./B)./(A.^0.5) - 2.0.*KB

    test = sign(Θ[end])
    j = length(r)
    while sign(Θ[j])==test && j>=2 # is the 1st ghost point on the left
    test = sign(Θ[j])
    j = j-1
    end

    # we take the midpoint where the sign of the expansion changes. If this is <0, then there is no AH
    rAH = (r[j] + r[j+1])/2
    rAH
end

# calculate the total initial energy that is Integral_rmin^rmax ρ*r^2*dr
# with ρ  = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)
function energy(sys::System, v_classic::Array)
    En = sum((0.5.*((v_classic[3:end,2].^2.0)./v_classic[3:end,4].^2 .+ v_classic[3:end,1].^2.0)./v_classic[3:end,4]).*sys.r[3:end].^2)*(sys.r[4]-sys.r[3])

    En
end

# # for tools
# function constraints_sbp21(Φ::Vector, Π::Vector, Ψ::Vector,
#                            A::Vector, B::Vector, DB::Vector,
#                            Utld::Vector, K::Vector, KB::Vector,
#                            λ::Vector, α::Vector, Dα::Vector, r::Vector,
#                            oMp2::Float64)

#     dr = r[2] - r[1]

#     # populate ghost points r[1]=-2h, r[2]=-h
#     # even | f[3] = f(r=0)
#     # f[1] = f[5] = f(r=2h)
#     # f[2] = f[4] = f(r=h)
#     Φ[1]  = Φ[5]
#     Φ[2]  = Φ[4]
#     Π[1]  = Π[5]
#     Π[2]  = Π[4]
#     A[1]  = A[5]
#     A[2]  = A[4]
#     B[1]  = B[5]
#     B[2]  = B[4]
#     K[1]  = K[5]
#     K[2]  = K[4]
#     KB[1] = KB[5]
#     KB[2] = KB[4]
#     α[1]  = α[5]
#     α[2]  = α[4]

#     # odd | f[3] = f(r=0)
#     # f[1] = -f[5] = -f(r=2h)
#     # f[2] = -f[4] = -f(r=h)
#     Ψ[1]    = -Ψ[5]
#     Ψ[2]    = -Ψ[4]
#     DB[1]   = -DB[5]
#     DB[2]   = -DB[4]
#     Utld[1] = -Utld[5]
#     Utld[2] = -Utld[4]
#     λ[1]    = -λ[5]
#     λ[2]    = -λ[4]
#     Dα[1]   = -Dα[5]
#     Dα[2]   = -Dα[4]

#     # FD2 centered r derivatives
#     λ_r    = Dr_FD2(λ, dr)
#     DB_r   = Dr_FD2(DB, dr)
#     Utld_r = Dr_FD2(Utld, dr)
#     K_r    = Dr_FD2(K, dr)
#     KB_r   = Dr_FD2(KB, dr)

#     Ψ_SBP21    = Dr_SBP2(Ψ, dr, 1)
#     λ_SBP21    = Dr_SBP2(λ, dr, 1)
#     DB_SBP21   = Dr_SBP2(DB, dr, 1)
#     Dα_SBP21   = Dr_SBP2(Dα, dr, 1)
#     Utld_SBP21 = Dr_SBP2(Utld, dr, 1)

#     # define 1/A and 1/B metric function; it is faster to multiply than devide
#     oA = 1.0 ./ A
#     oB = 1.0 ./ B
#     # construct DA:
#     DA = @. Utld + 2.0*DB + 4.0*B*λ*oA
#     # for quantum modes
#     DA_r   = Dr_FD2(DA, dr)
#     DB_r   = Dr_FD2(DB, dr)
#     Dα_r   = Dr_FD2(Dα, dr)
#     # define 1/Mp^2 = 8 π G, where G is Netwon's constant and is set to 1 here
#     #oMp2 = 1.0 # 8*π
#     #Λ = 0.0

#     # stress-energy tensor components:
#     jA = -(oA.^0.5).*oB.*Π.*Ψ
#     SA = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)
#     SB = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .- Ψ.^2.0)
#     ρ  = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)

#     # Hamiltonian constraint
#     Ham = @. -0.25 * DB * (DB + 2. * Utld) - 2. * B * DB * oA * λ +
#     A * (3. * KB^2 - 2. * KB * K + oMp2 * ρ) + (DB_SBP21) + (Utld_r) - 
#     1. * (Utld_SBP21) - 1. * (A - 4. * B) * oA * ((λ_r) - 1. * (λ_SBP21))

#     # Momentum constraint
#     Kminus3KB       = K - 3.0*KB
#     Kminus3KB_SBP21 = Dr_SBP2(Kminus3KB, dr, 1)
#     Kminus3KB_r     = Dr_FD2(Kminus3KB, dr)
#     Mom = @. -0.5 * DB * Kminus3KB + 0.5 * jA * oMp2 + (KB_r) + (Kminus3KB_r) - 1. * (Kminus3KB_SBP21)

#     Ham, Mom
# end



# the function lists all h5 files in a directory with prefix e.g. "data_"
function list_h5_files(foldername::String; prefix::String)
    path     = abspath(foldername)
    allfiles = readdir(path)

    Ns = length(prefix)

    its_names = Tuple[]
    # append only the files whose names start with the given prefix
    for file in allfiles
        try
            if (file[1:Ns] == prefix && file[end-2:end] == ".h5")
                # extract iteration
                it_str = file[Ns+1:end-3]
                fullname = joinpath(path, file)
                # add to list of tuples with iteration and name
                push!(its_names, (parse(Int64, it_str), fullname))
            end
        catch ex
            if isa(ex, BoundsError)
                # probably triggered by string comparison; do nothing
            else
                throw(ex)
            end
        end
    end

    # sort according to iteration
    sort!(its_names)
    # and extract the list of filenames and iterations
    filenames = [name for (it, name) in its_names]
    its       = [it for (it, name) in its_names]

    (its, filenames)
end

function Ricci_scalar(v_classic::Array, r::Vector, oMp2::Float64)

    dr = r[2] - r[1]

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

    # odd | f[3] = f(r=0)
    # f[1] = -f[5] = -f(r=2h)
    # f[2] = -f[4] = -f(r=h)
    Ψ    = odd_ghosts(Ψ)
    DB   = odd_ghosts(DB)
    Utld = odd_ghosts(Utld)
    λ    = odd_ghosts(λ)
    Dα   = odd_ghosts(Dα)

    # define 1/A and 1/B metric function; it is faster to multiply than devide
    oA = 1.0 ./ A
    oB = 1.0 ./ B
    # construct DA:
    DA = @. Utld + 2.0*DB + 4.0*B*λ*oA

    # FD2 centered r derivatives
    λ_r    = Dr_FD2(λ, dr)
    DB_r   = Dr_FD2(DB, dr)
    Utld_r = Dr_FD2(Utld, dr)
    K_r    = Dr_FD2(K, dr)
    KB_r   = Dr_FD2(KB, dr)

    Ψ_SBP21    = Dr_SBP2(Ψ, dr, 1)
    λ_SBP21    = Dr_SBP2(λ, dr, 1)
    DB_SBP21   = Dr_SBP2(DB, dr, 1)
    DA_SBP21   = Dr_SBP2(DA, dr, 1)
    Dα_SBP21   = Dr_SBP2(Dα, dr, 1)
    Utld_SBP21 = Dr_SBP2(Utld, dr, 1)

    # for quantum modes
    DA_r   = Dr_FD2(DA, dr)
    DB_r   = Dr_FD2(DB, dr)
    Dα_r   = Dr_FD2(Dα, dr)

    Λ = 0.0

    # stress-energy tensor components:
    jA = -(oA.^0.5).*oB.*Π.*Ψ
    SA = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)
    SB = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .- Ψ.^2.0)
    ρ  = 0.5.*oA.*((oB.^2.0).*(Π.^2.0) .+ Ψ.^2.0)

    ric = zeros(length(A))

    # Ricci scalar
    ric[4:end] = @. (6. * KB^2 - 4. * KB * K - 1. * DA * DB * oA + 
    1.5 * DB^2 * oA + (2. * oA)/r^2 - (2. * oB)/r^2 - 25.132741228718345 * SA - 
    50.26548245743669 * SB + 75.39822368615503 * ρ + 2. * oA * (DA_r) -
    2. * oA * (DA_SBP21) - 4. * oA * (DB_r) + 6. * oA * (DB_SBP21))[4:end]

    ric[3] = 6. * KB[3]^2 - 4. * KB[3] * K[3] - 1. * DA[3] * DB[3] * oA[3] + 1.5 * DB[3]^2 * oA[3] -
    25.132741228718345 * SA[3] - 50.26548245743669 * SB[3] + 75.39822368615503 * ρ[3] +
    2. * oA[3] * (DA_r[3]) + 2. * oB[3] * (DA_r[3]) - 2. * oA[3] * (DA_SBP21[3]) -
    2. * oB[3] * (DA_SBP21[3]) - 6. * oA[3] * (DB_r[3]) + 8. * oA[3] * (DB_SBP21[3])

    ric
end
