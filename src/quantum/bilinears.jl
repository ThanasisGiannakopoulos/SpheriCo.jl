# start bilinears; only for quantum case
# MULTIPLE DISPATCH

# for the stress tensor (classical+quantum) aka bilinears
# mutates bilin

# no regularization
function bilinears(t::Float64,
                   v_classic::Array{Float64, 2}, v_quantum::Array{ComplexF64, 4},
                   p::Param, r::Vector, bilin::Array{ComplexF64})

    # load classical scalar field vars, for bilin_classic
    Φ = v_classic[:, 1]
    Π = v_classic[:, 2]
    Ψ = v_classic[:, 3]
    # load needed geometric quantities
    α = v_classic[:, 11]
    A = v_classic[:, 4]
    B = v_classic[:, 5]

    # buffers for the sums in multithreading with
    # format (r, number_of_quantum_bilins, number_of_threads)
    buffers = zeros(ComplexF64, (length(bilin[:,1]), 5, Threads.nthreads()))

    # sum for ll=0
    @threads for kk in 1:p.kmax
        id = Threads.threadid()

        # utld_kl = r^l*uq; for l=0, utld_kl = uq
        utld_kl = v_quantum[:, 1, Int(kk), 1]
        # ∂r_utld_kl = l*r^(l-1)*uq + r^l*ψq = ψq for l=0
        dr_utld_kl = v_quantum[:, 1, Int(kk), 2]  
        πq = v_quantum[:, 1, Int(kk), 3]
        # ∂t_utld_kl = r^l*α*πq/(B*A^0.5)
        dt_utld_kl = @. α* πq/(B* A^0.5)

        # bilin[:,1]
        # (4*π)/(hbar*c^2)* [μ^2*<Φ^2>]_quantum = Sum_{kl} μ^2*(2l+1)*|utld_kl|^2;
        # where μ=m is the quantum scalar mass here μ=0,
        # so it has a zero contribution
        buffers[:,1,id] = buffers[:,1,id] + zeros(ComplexF64, length(bilin[:,1]))
        # bilin[:,2]
        # (4*π*α^2)/(hbar*c^2*A*B^2)* <Π^2>_quantum = Sum_{kl} dk*(2*l+1)*|∂_t utld_kl|^2
        # l=0 here
        buffers[:,2,id] .+= @. dt_utld_kl*conj(dt_utld_kl)*p.dk
        # bilin[:,3]
        # (4*π)/(hbar*c^2)* <Ψ^2>_quantum = Sum_{kl} dk*(2l+1)*|∂_r utld_kl|^2
        # l=0 here
        buffers[:,3,id] .+= @. dr_utld_kl*conj(dr_utld_kl)*p.dk
        # bilin[:,4]
        # (4*π*α)/(hbar*c^2*A^(1/2)*B)* <Π*Ψ>_quantum =
        #  Sum_{kl} dk*0.5*(2*l+1)[(∂_r utld_kl)*conj(∂_t utld_kl) +
        #    conj(∂_r utld_kl)(∂_t utld_kl)]
        # l=0 here
        buffers[:,4,id] .+= @. p.dk*(0.5*dr_utld_kl*conj(dt_utld_kl) +
            0.5*dt_utld_kl*conj(dr_utld_kl))
        # bilin[:,5]
        # (4*π)/(hbar*c^2)* [(1/r^2)*<∂_θ Φ^2>]_quantum =
        #  Sum_{kl} dk*0.5*l*(l+1)*(2l+1)*|r^(l-1)*u_kl|^2;
        # l=0 here so this contribution is 0
        buffers[:,5,id] = buffers[:,5,id] .+ zeros(ComplexF64, length(bilin[:,5]))
    end
    # end loops in k

    # sum for l>0
    @threads for l in 1:p.lmax
        for k in 1:p.kmax
            id = Threads.threadid()

            uq = v_quantum[:, Int(l+1), Int(k), 1]
            ψq = v_quantum[:, Int(l+1), Int(k), 2]  
            πq = v_quantum[:, Int(l+1), Int(k), 3]
            # ∂t_utld_kl = r^l*α*πq/(B*A^0.5)
            dt_utld_kl = @. r^l* α* πq/(B* A^0.5)
            # ∂r_utld_kl = l*r^(l-1)*uq + r^l*ψq
            dr_utld_kl = @. l*r^(l-1)*uq + ψq*r^l             
            # utld_kl = r^l*uq
            utld_kl = @. r^l*uq
            # rlminus1_u_kl := r^(l-1)*u_kl
            rlminus1_u_kl = r.^(l-1.0).*uq

            # bilin[:,1]
            # (4*π)/(hbar*c^2)* [μ^2*<Φ^2>]_quantum = Sum_{kl} dk*μ^2*(2l+1)*|utld_kl|^2;
            # where μ=0 here, so this contribution is zero
            buffers[:,1,id] = buffers[:,1,id] + zeros(ComplexF64, length(bilin[:,1]))
            # bilin[:,2]
            # (4*π*α^2)/(hbar*c^2*A*B^2)* <Π^2>_quantum = Sum_{kl} dk*(2*l+1)*|∂_t utld_kl|^2;
            buffers[:,2,id] .+= @. p.dk*(2.0*l+1)*(dt_utld_kl*conj(dt_utld_kl))
            # bilin[:,3]
            # (4*π)/(hbar*c^2)* <Ψ^2>_quantum = Sum_{kl} dk*(2l+1)*|∂_r utld_kl|^2;
            buffers[:,3,id] .+= @. p.dk*(2.0*l+1.0)*(dr_utld_kl*conj(dr_utld_kl))
            # bilin[:,4]
            # (4*π*α)/(hbar*c^2*A^(1/2)*B)* <Π*Ψ>_quantum = Sum_{kl} dk*0.5*(2*l+1)*[(∂_r utld_kl)*conj(∂_t utld_kl) + conj(∂_r utld_kl)(∂_t utld_kl)];
            buffers[:,4,id] .+= @. p.dk* (2.0*l+1.0)*(0.5*dr_utld_kl*conj(dt_utld_kl) +
                0.5*dt_utld_kl*conj(dr_utld_kl))
            # bilin[:,5]
            # (4*π)/(hbar*c^2)* [(1/r^2)*<∂_θ Φ^2>]_quantum = Sum_{kl} dk*0.5*l*(l+1)*(2l+1)*|r^(l-1)*u_kl|^2;
            buffers[:,5,id] .+= @. p.dk*l*(l+1.0)*(2.0*l+1.0)*(
                0.5*rlminus1_u_kl*conj(rlminus1_u_kl))
        end # end loops in k
    end # end loops in l

    bilin[:,1] .= sum(buffers[:,1,:], dims=2)
    bilin[:,2] .= sum(buffers[:,2,:], dims=2)
    bilin[:,3] .= sum(buffers[:,3,:], dims=2)
    bilin[:,4] .= sum(buffers[:,4,:], dims=2)
    bilin[:,5] .= sum(buffers[:,5,:], dims=2)

    nothing
end

# with regularization
function bilinears(t::Float64,
                   v_classic::Array{Float64, 2}, v_quantum::Array{ComplexF64, 5},
                   p::Param, r::Vector, bilin::Array{ComplexF64, 2})

    # load classical scalar field vars, for bilin_classic
    Φ = v_classic[:, 1]
    Π = v_classic[:, 2]
    Ψ = v_classic[:, 3]
    # load needed geometric quantities
    α = v_classic[:, 11]
    A = v_classic[:, 4]
    B = v_classic[:, 5]

    # buffers for the sums in multithreading with
    # format (r, number_of_quantum_bilins, number_of_threads)
    buffers = zeros(ComplexF64, (length(bilin[:,1]), 5, Threads.nthreads()))

    # sum for ll=0
    @threads for kk in 1:p.kmax
        id = Threads.threadid()

        # utld_kl = r^l*uq; for l=0, utld_kl = uq
        utld_kl_m0 = v_quantum[:, 1, Int(kk), 1, 1]
        utld_kl_m1 = v_quantum[:, 1, Int(kk), 2, 1]
        utld_kl_m2 = v_quantum[:, 1, Int(kk), 3, 1]
        utld_kl_m5 = v_quantum[:, 1, Int(kk), 4, 1]

        # ∂r_utld_kl = l*r^(l-1)*uq + r^l*ψq = ψq for l=0
        dr_utld_kl_m0 = v_quantum[:, 1, Int(kk), 1, 2]
        dr_utld_kl_m1 = v_quantum[:, 1, Int(kk), 2, 2]
        dr_utld_kl_m2 = v_quantum[:, 1, Int(kk), 3, 2]
        dr_utld_kl_m5 = v_quantum[:, 1, Int(kk), 4, 2]

        πq_m0 = v_quantum[:, 1, Int(kk), 1, 3]
        πq_m1 = v_quantum[:, 1, Int(kk), 2, 3]
        πq_m2 = v_quantum[:, 1, Int(kk), 3, 3]
        πq_m5 = v_quantum[:, 1, Int(kk), 4, 3]

        # ∂t_utld_kl = r^l*α*πq/(B*A^0.5)
        dt_utld_kl_m0 = @. α* πq_m0/(B* A^0.5)
        dt_utld_kl_m1 = @. α* πq_m1/(B* A^0.5)
        dt_utld_kl_m2 = @. α* πq_m2/(B* A^0.5)
        dt_utld_kl_m5 = @. α* πq_m5/(B* A^0.5)

        # bilin[:,1]
        # (4*π)/(hbar*c^2)* [μ^2*<Φ^2>]_quantum = Sum_{klm} dk*μ^2[m]*weight[m]*|utld_klm|^2;
        # where μ=m is the quantum scalar mass here p.mlist = [0*mPV, 1*mPV, sqrt(3)*mPV, 2*mPV]
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
    end # end loops in k

    # sum for ll>0
    @threads for ll in 1:p.lmax
        for kk in 1:p.kmax
            id = Threads.threadid()

            uq_m0 = v_quantum[:, Int(ll+1), Int(kk), 1, 1]
            uq_m1 = v_quantum[:, Int(ll+1), Int(kk), 2, 1]
            uq_m2 = v_quantum[:, Int(ll+1), Int(kk), 3, 1]
            uq_m5 = v_quantum[:, Int(ll+1), Int(kk), 4, 1]

            ψq_m0 = v_quantum[:, Int(ll+1), Int(kk), 1, 2]
            ψq_m1 = v_quantum[:, Int(ll+1), Int(kk), 2, 2]
            ψq_m2 = v_quantum[:, Int(ll+1), Int(kk), 3, 2]
            ψq_m5 = v_quantum[:, Int(ll+1), Int(kk), 4, 2]

            πq_m0 = v_quantum[:, Int(ll+1), Int(kk), 1, 3]
            πq_m1 = v_quantum[:, Int(ll+1), Int(kk), 2, 3]
            πq_m2 = v_quantum[:, Int(ll+1), Int(kk), 3, 3]
            πq_m5 = v_quantum[:, Int(ll+1), Int(kk), 4, 3]

            # rlminus1_u_kl_mi := r^(l-1)*u_klm
            rlminus1_u_kl_m0 = @. r^(ll-1.0)*uq_m0
            rlminus1_u_kl_m1 = @. r^(ll-1.0)*uq_m1
            rlminus1_u_kl_m2 = @. r^(ll-1.0)*uq_m2
            rlminus1_u_kl_m5 = @. r^(ll-1.0)*uq_m5

            # ∂t_utld_kl = r^l*α*πq/(B*A^0.5)
            dt_utld_kl_m0 = @. r^ll* α* πq_m0/(B* A^0.5)
            dt_utld_kl_m1 = @. r^ll* α* πq_m1/(B* A^0.5)
            dt_utld_kl_m2 = @. r^ll* α* πq_m2/(B* A^0.5)
            dt_utld_kl_m5 = @. r^ll* α* πq_m5/(B* A^0.5)

            # ∂r_utld_kl = l*r^(l-1)*uq + r^l*ψq
            dr_utld_kl_m0 = @. (ll*r^(ll-1)*uq_m0 + ψq_m0*r^ll)
            dr_utld_kl_m1 = @. (ll*r^(ll-1)*uq_m1 + ψq_m1*r^ll)
            dr_utld_kl_m2 = @. (ll*r^(ll-1)*uq_m2 + ψq_m2*r^ll)
            dr_utld_kl_m5 = @. (ll*r^(ll-1)*uq_m5 + ψq_m5*r^ll)

            # utld_kl = r^l*uq
            utld_kl_m0 = @. r^ll*uq_m0
            utld_kl_m1 = @. r^ll*uq_m1
            utld_kl_m2 = @. r^ll*uq_m2
            utld_kl_m5 = @. r^ll*uq_m5

            # bilin[:,1]
            # (4*π)/(hbar*c^2)* [μ^2*<Φ^2>]_quantum = Sum_{klm} dk*μ^2[m]*weight[m]*|utld_klm|^2;
            # where μ=m is the quantum scalar mass, p.mlist = [0*mPV, 1*mPV, sqrt(3)*mPV, 2*mPV]
            # weight[m] = [1,-2,2,-1]
            buffers[:,1,id] .+= @. p.dk*(2.0*ll+1.0)*(p.mlist[1]^2 *utld_kl_m0*conj(utld_kl_m0) -
                p.mlist[2]^2 *2.0*utld_kl_m1*conj(utld_kl_m1) +
                p.mlist[3]^2 *2.0*utld_kl_m2*conj(utld_kl_m2) -
                p.mlist[4]^2 *utld_kl_m5*conj(utld_kl_m5))
            # bilin[:,2]
            # (4*π*α^2)/(hbar*c^2*A*B^2)* <Π^2>_quantum =
            #  Sum_{klm} dk*(2*l+1)*weight[m]*|∂_t utld_klm|^2;
            # weight[m] = [1,-2,2,-1]
            buffers[:,2,id] .+= @. p.dk*(2.0*ll+1.0)*(dt_utld_kl_m0*conj(dt_utld_kl_m0) -
                2.0*dt_utld_kl_m1*conj(dt_utld_kl_m1) +
                2.0*dt_utld_kl_m2*conj(dt_utld_kl_m2) -
                dt_utld_kl_m5*conj(dt_utld_kl_m5))
            # bilin[:,3]
            # (4*π)/(hbar*c^2)* <Ψ^2>_quantum = Sum_{klm} dk*(2l+1)*weight[m]*|∂_r utld_klm|^2;
            # weight[m] = [1,-2,2,-1]
            buffers[:,3,id] .+= @. p.dk*(2.0*ll+1.0)*(dr_utld_kl_m0*conj(dr_utld_kl_m0) -
                2.0*dr_utld_kl_m1*conj(dr_utld_kl_m1) +
                2.0*dr_utld_kl_m2*conj(dr_utld_kl_m2) -
                dr_utld_kl_m5*conj(dr_utld_kl_m5))
            # bilin[:,4]
            # (4*π*α)/(hbar*c^2*A^(1/2)*B)* <Π*Ψ>_quantum =
            #  Sum_{klm} dk*0.5*(2*l+1)*weight[m]*[(∂_r utld_klm)*conj(∂_t utld_klm) +
            #    conj(∂_r utld_klm)(∂_t utld_klm)];
            # weight[m] = [1,-2,2,-1]
            buffers[:,4,id] .+= @. p.dk*(2.0*ll+1.0)*(0.5*dr_utld_kl_m0*conj(dt_utld_kl_m0) +
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
            buffers[:,5,id] .+= @. p.dk*ll*(ll+1.0)*(2.0*ll+1.0)*(
                0.5*rlminus1_u_kl_m0*conj(rlminus1_u_kl_m0) -
                rlminus1_u_kl_m1*conj(rlminus1_u_kl_m1) +
                rlminus1_u_kl_m2*conj(rlminus1_u_kl_m2) -
                0.5*rlminus1_u_kl_m5*conj(rlminus1_u_kl_m5))
        end # end loops in k
    end # end loops in l

    bilin[:,1] .= sum(buffers[:,1,:], dims=2)
    bilin[:,2] .= sum(buffers[:,2,:], dims=2)
    bilin[:,3] .= sum(buffers[:,3,:], dims=2)
    bilin[:,4] .= sum(buffers[:,4,:], dims=2)
    bilin[:,5] .= sum(buffers[:,5,:], dims=2)

    nothing
end

# end bilinears
