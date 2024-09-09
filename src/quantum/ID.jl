# start quantum_ID
# MULTIPLE DISPATCH

# quantum fields initial data
function uq_ID(r::Float64, k::Float64, m::Float64, l::Float64)
    ω = sqrt(k^2 + m^2)
    if r==0.0
        (k^(l+1)/sqrt(ω))/(gamma(1.5 + l)*2.0^(1+l))
    else
        (k/sqrt(π*ω))*sphericalbesselj(Int(l), k*r)/r^l
    end
end
function ψq_ID(r::Float64, k::Float64, m::Float64, l::Float64)
    ω = sqrt(k^2 + m^2)
    if r==0.0
        0.0 + 0.0im
    else
        (k/sqrt(π*ω))*0.5*(
            k*r*sphericalbesselj(Int(l)-1,k*r) -
            (2*l+1)*sphericalbesselj(Int(l),k*r) -
            k*r*sphericalbesselj(Int(l)+1,k*r) 
        )/(r^(l+1))
    end
end
function πq_ID(r::Float64, k::Float64, m::Float64, l::Float64)
    ω = sqrt(k^2 + m^2)
    if r==0.0
        -((sqrt(ω)*k^(l+1))/(gamma(1.5 + l)*2.0^(1+l)))im
    else
        -(k*sqrt(ω/π)*sphericalbesselj(Int(l), k*r)/r^l)im
    end
end

# non-regularized version (only mass=0):
# one field for each k in [1,kmax], l in [0,lmax], 1 + 2 reduction vars 
# v_quantum = zeros(ComplexF64, (length(rr), Int(p.lmax+1), Int(p.kmax), 3) )
function quantum_ID(v_quantum::Array{ComplexF64, 4}, sys::System, p::Param)

    # number or points in r-grid
    Nr = length(sys.r)
    # mass of quantum fields
    m = 0.0

    """
    The ID for the quantum modes are organized in the v matrix of size
    v_quantum = size(Nr, (lmax+1), kmax, 3),
    where the lamx+1 is because l includes l=0, wheareas k starts for
    k=1, in the following way:

    v_quantum[:,1, 1, 1:3]  are 3 quantum fields for l=0, k=1
    v_quantum[:,1, 2, 1:3]  are 3 quantum fields for l=0, k=2
    """
    l_checked = zeros(Threads.nthreads()) # to check the progress of generating ID
    @threads for ll in 0:p.lmax
        id = Threads.threadid()
        l_checked[id] += 1.0
        print("\rgenerating quantum ID: $(round(sum(l_checked)/(p.lmax+1.0)*100.0, digits=2)) %")
        for kk in 1:p.kmax
            # index for uq, ψq, and πq
            uqi = 1
            ψqi = 2
            πqi = 3
            for i in 3:Nr
                # find the k for this l, for which the ψ ID is zero at rmax
                k=kk*p.dk
                # uq; even
                v_quantum[i, Int(ll+1), Int(kk), uqi] = uq_ID(sys.r[i], k, m, Float64(ll) )
                # ψq; odd
                v_quantum[i, Int(ll+1), Int(kk), ψqi] = ψq_ID(sys.r[i], k, m, Float64(ll) )
                # πq; even
                v_quantum[i, Int(ll+1), Int(kk), πqi] = πq_ID(sys.r[i], k, m, Float64(ll) )
            end
            # uq ghosts; even
            v_quantum[:, Int(ll+1), Int(kk), uqi] = even_ghosts(v_quantum[:, Int(ll+1), Int(kk), uqi])
            # ψq ghosts; odd
            v_quantum[:, Int(ll+1), Int(kk), ψqi] = odd_ghosts(v_quantum[:, Int(ll+1), Int(kk), ψqi])
            # πq ghosts; even
            v_quantum[:, Int(ll+1), Int(kk), πqi] = even_ghosts(v_quantum[:, Int(ll+1), Int(kk), πqi])
        end
    end
    println()

    v_quantum
end

# regularized version: one for each k in [1,kmax], l in [0,lmax], and m in p.mlist
# v_quantum = zeros(ComplexF64, (length(rr), Int(p.lmax+1), Int(p.kmax), length(p.mlist), 3) )
function quantum_ID(v_quantum::Array{ComplexF64, 5}, sys::System, p::Param)

    # number or points in r-grid
    Nr = length(sys.r)

    """ 
    The ID for the quantum modes are organized in the v matrix of size
    v_quantum = size(Nr, (lmax+1), kmax, length(p.mlist), 3),
    where the lamx+1 is because l includes l=0, wheareas k starts for
    k=1, in the following way:

    v_quantum[:,1, 1, 1, 1:3]  are 3 quantum fields for l=0, k=1, mi=m0
    v_quantum[:,1, 1, 2, 1:3]  are 3 quantum fields for l=0, k=1, mi=m1
        """
    l_checked = zeros(Threads.nthreads()) # to check the progress of generating ID
    @threads for ll in 0:p.lmax
        id = Threads.threadid()
        l_checked[id] += 1.0
        print("\rgenerating quantum ID: $(round(sum(l_checked)/(p.lmax+1.0)*100.0, digits=2)) %")
        for kk in 1:p.kmax
            for mm in 1:length(p.mlist)
                # index for uq, ψq, and πq
                uqi = 1
                ψqi = 2
                πqi = 3
                mi = p.mlist[mm]
                for i in 3:Nr
                    k=kk*p.dk
                    # uq; even
                    v_quantum[i, Int(ll+1), Int(kk), mm, uqi] = uq_ID(sys.r[i], k, mi, Float64(ll) )
                    # ψq; odd
                    v_quantum[i, Int(ll+1), Int(kk), mm, ψqi] = ψq_ID(sys.r[i], k, mi, Float64(ll) )
                    # πq; even
                    v_quantum[i, Int(ll+1), Int(kk), mm, πqi] = πq_ID(sys.r[i], k, mi, Float64(ll) )
                end
                # uq ghosts; even
                v_quantum[:, Int(ll+1), Int(kk), mm, uqi] = even_ghosts(v_quantum[:, Int(ll+1), Int(kk), mm, uqi])
                # ψq ghosts; odd
                v_quantum[:, Int(ll+1), Int(kk), mm, ψqi] = odd_ghosts(v_quantum[:, Int(ll+1), Int(kk), mm, ψqi])
                # πq ghosts; even
                v_quantum[:, Int(ll+1), Int(kk), mm, πqi] = even_ghosts(v_quantum[:, Int(ll+1), Int(kk), mm, πqi])
            end
        end
    end
    println()

    return v_quantum
end
