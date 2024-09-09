# start correlators
# MULTIPLE DISPATCH

function correlators(v_classic::Array{Float64, 2}, v_quantum::Array{ComplexF64, 4},
                     p::Param, r::Vector, correlators::Array{ComplexF64})

    # load needed geometric quantities
    α = v_classic[:, 11]
    A = v_classic[:, 4]
    B = v_classic[:, 5]
    # length of r-grid
    Nr = length(A)

    # buffers for the sums in multithreading with
    # format (r, r, number_of_correlators, number_of_threads)
    buffers = zeros(ComplexF64, (length(A), length(A), 5, Threads.nthreads()))

    # sum for ll=0
    @threads for kk in 1:p.kmax
        id = Threads.threadid()

        # utld_kl = r^l*uq; for l=0, utld_kl = uq
        uq =  v_quantum[:, 1, Int(kk), 1]
        utld_kl = uq
        # ∂r_utld_kl = l*r^(l-1)*uq + r^l*ψq; for l=0 it is just ψq
        ψq = v_quantum[:, 1, Int(kk), 2] 
        dr_utld_kl = ψq
        # ∂t_utld_kl = r^l*α*πq/(B*A^0.5)
        πq = v_quantum[:, 1, Int(kk), 3]
        dt_utld_kl = @. α* πq/(B* A^0.5)

        for i in 1:Nr
            for j in 1:Nr
                #correlators[i,j,1]
                buffers[i,j,1,id] += p.dk*dt_utld_kl[i]*conj(dt_utld_kl[j])

                #correlators[i,j,2]
                buffers[i,j,2,id] += p.dk*dr_utld_kl[i]*conj(dr_utld_kl[j])

                #correlators[i,j,3]
                buffers[i,j,3,id] += p.dk*utld_kl[i]*conj(utld_kl[j])

                #correlators[i,j,4]
                buffers[i,j,4,id] += p.dk*πq[i]*conj(πq[j])

                #correlators[i,j,5]
                buffers[i,j,5,id] += p.dk*ψq[i]*conj(ψq[j])
            end # end loops in Nr
        end # end loops in Nr
    end # end loops in k

    # sum for l>0
    @threads for l in 1:p.lmax-1
        id = Threads.threadid()
        for k in 1:p.kmax
            uq = v_quantum[:, Int(l)+1, Int(k), 1]
            ψq = v_quantum[:, Int(l)+1, Int(k), 2]
            πq = v_quantum[:, Int(l)+1, Int(k), 3]

            # ∂t_utld_kl = r^l*α*πq/(B*A^0.5)
            dt_utld_kl = @. r^l* α* πq/(B* A^0.5)
            # ∂r_utld_kl = l*r^(l-1)*uq + r^l*ψq
            dr_utld_kl = @. l*r^(l-1)*uq + ψq*r^l             
            utld_kl = @. r^l*uq
            πtld_kl = @. r^l*πq
            ψtld_kl = @. r^l*ψq

            for i in 1:Nr
                for j in 1:Nr
                    #correlators[i,j,1]
                    buffers[i,j,1,id] += p.dk*dt_utld_kl[i]*conj(dt_utld_kl[j])*(2.0*l+1)

                    #correlators[i,j,2]
                    buffers[i,j,2,id] += p.dk*dr_utld_kl[i]*conj(dr_utld_kl[j])*(2.0*l+1)

                    #correlators[i,j,3]
                    buffers[i,j,3,id] += p.dk*utld_kl[i]*conj(utld_kl[j])*(2.0*l+1)

                    #correlators[i,j,4]
                    buffers[i,j,4,id] +=  p.dk*πtld_kl[i]*conj(πtld_kl[j])

                    #correlators[i,j,5]
                    buffers[i,j,5,id] +=  p.dk*ψtld_kl[i]*conj(ψtld_kl[j])
                end # end loops in Nr
            end # end loops in Nr
        end # end loops in k
    end # end loops in l

    correlators[:,:,1] .= sum(buffers[:,:,1,:], dims=3)
    correlators[:,:,2] .= sum(buffers[:,:,2,:], dims=3)
    correlators[:,:,3] .= sum(buffers[:,:,3,:], dims=3)
    correlators[:,:,4] .= sum(buffers[:,:,4,:], dims=3)
    correlators[:,:,5] .= sum(buffers[:,:,5,:], dims=3)

    nothing
end

#mutates correlators
function correlators(v_classic::Array{Float64, 2}, v_quantum::Array{ComplexF64, 5},
                     p::Param, r::Vector, correlators::Array{ComplexF64})

    # load needed geometric quantities
    α = v_classic[:, 11]
    A = v_classic[:, 4]
    B = v_classic[:, 5]
    # length of r-grid
    Nr = length(A)

    # buffers for the sums in multithreading with
    # format (r, r, number_of_correlators, number_of_threads)
    buffers = zeros(ComplexF64, (length(A), length(A), 5, Threads.nthreads()))

    # sum for l=0
    @threads for k in 1:p.kmax
        id = Threads.threadid()

        # utld_kl = r^l*uq; for l=0, utld_kl = uq
        utld_kl_m0 = v_quantum[:, 1, Int(k), 1, 1]
        utld_kl_m1 = v_quantum[:, 1, Int(k), 2, 1]
        utld_kl_m2 = v_quantum[:, 1, Int(k), 3, 1]
        utld_kl_m5 = v_quantum[:, 1, Int(k), 4, 1]

        # ∂r_utld_kl = l*r^(l-1)*uq + r^l*ψq; for l=0 it is just ψq
        ψq_m0 = v_quantum[:, 1, Int(k), 1, 2]
        ψq_m1 = v_quantum[:, 1, Int(k), 2, 2]
        ψq_m2 = v_quantum[:, 1, Int(k), 3, 2]
        ψq_m5 = v_quantum[:, 1, Int(k), 4, 2]
        dr_utld_kl_m0 = ψq_m0
        dr_utld_kl_m1 = ψq_m1
        dr_utld_kl_m2 = ψq_m2
        dr_utld_kl_m5 = ψq_m5

        # ∂t_utld_kl = r^l*α*πq/(B*A^0.5)
        πq_m0 = v_quantum[:, 1, Int(k), 1, 3]
        πq_m1 = v_quantum[:, 1, Int(k), 2, 3]
        πq_m2 = v_quantum[:, 1, Int(k), 3, 3]
        πq_m5 = v_quantum[:, 1, Int(k), 4, 3]
        dt_utld_kl_m0 = α.* πq_m0./(B.* A.^0.5)
        dt_utld_kl_m1 = α.* πq_m1./(B.* A.^0.5)
        dt_utld_kl_m2 = α.* πq_m2./(B.* A.^0.5)
        dt_utld_kl_m5 = α.* πq_m5./(B.* A.^0.5)

        for i in 1:Nr
            for j in 1:Nr
                #correlators[i,j,1]
                buffers[i,j,1,id] += p.dk*(
                    dt_utld_kl_m0[i]*conj(dt_utld_kl_m0[j]) -
                    2.0*dt_utld_kl_m1[i]*conj(dt_utld_kl_m1[j]) +
                    2.0*dt_utld_kl_m2[i]*conj(dt_utld_kl_m2[j]) -
                    dt_utld_kl_m5[i]*conj(dt_utld_kl_m5[j]))

                #correlators[i,j,2]
                buffers[i,j,2,id] += p.dk*(
                    dr_utld_kl_m0[i]*conj(dr_utld_kl_m0[j]) -
                    2.0*dr_utld_kl_m1[i]*conj(dr_utld_kl_m1[j]) +
                    2.0*dr_utld_kl_m2[i]*conj(dr_utld_kl_m2[j]) -
                    dr_utld_kl_m5[i]*conj(dr_utld_kl_m5[j]))

                #correlators[i,j,3]
                buffers[i,j,3,id] += p.dk*(
                    utld_kl_m0[i]*conj(utld_kl_m0[j]) -
                    2.0*utld_kl_m1[i]*conj(utld_kl_m1[j]) +
                    2.0*utld_kl_m2[i]*conj(utld_kl_m2[j]) -
                    utld_kl_m5[i]*conj(utld_kl_m5[j]))

                #correlators[i,j,4]
                buffers[i,j,4,id] += p.dk*(πq_m0[i]*conj(πq_m0[j]) -
                    2.0*πq_m1[i]*conj(πq_m1[j]) +
                    2.0*πq_m2[i]*conj(πq_m2[j]) -
                    πq_m5[i]*conj(πq_m5[j]))

                #correlators[i,j,5]
                buffers[i,j,5,id] += p.dk*(ψq_m0[i]*conj(ψq_m0[j]) -
                    2.0*ψq_m1[i]*conj(ψq_m1[j]) +
                    2.0*ψq_m2[i]*conj(ψq_m2[j]) -
                    ψq_m5[i]*conj(ψq_m5[j]))
            end # end loops in Nr
        end # end loops in Nr
    end # end loop in k

    # sum for l>0
    @threads for l in 1:p.lmax
        id = Threads.threadid()
        for k in 1:p.kmax
            uq_m0 = v_quantum[:, Int(l)+1, Int(k), 1, 1]
            uq_m1 = v_quantum[:, Int(l)+1, Int(k), 2, 1]
            uq_m2 = v_quantum[:, Int(l)+1, Int(k), 3, 1]
            uq_m5 = v_quantum[:, Int(l)+1, Int(k), 4, 1]

            ψq_m0 = v_quantum[:, Int(l)+1, Int(k), 1, 2]
            ψq_m1 = v_quantum[:, Int(l)+1, Int(k), 2, 2]
            ψq_m2 = v_quantum[:, Int(l)+1, Int(k), 3, 2]
            ψq_m5 = v_quantum[:, Int(l)+1, Int(k), 4, 2]

            πq_m0 = v_quantum[:, Int(l)+1, Int(k), 1, 3]
            πq_m1 = v_quantum[:, Int(l)+1, Int(k), 2, 3]
            πq_m2 = v_quantum[:, Int(l)+1, Int(k), 3, 3]
            πq_m5 = v_quantum[:, Int(l)+1, Int(k), 4, 3]

            # ∂t_utld_kl = r^l*α*πq/(B*A^0.5)
            dt_utld_kl_m0 = @. r^l* α* πq_m0/(B* A^0.5)
            dt_utld_kl_m1 = @. r^l* α* πq_m1/(B* A^0.5)
            dt_utld_kl_m2 = @. r^l* α* πq_m2/(B* A^0.5)
            dt_utld_kl_m5 = @. r^l* α* πq_m5/(B* A^0.5)

            # ∂r_utld_kl = l*r^(l-1)*uq + r^l*ψq
            dr_utld_kl_m0 = @. l*r^(l-1)*uq_m0 + ψq_m0*r^l
            dr_utld_kl_m1 = @. l*r^(l-1)*uq_m1 + ψq_m1*r^l
            dr_utld_kl_m2 = @. l*r^(l-1)*uq_m2 + ψq_m2*r^l
            dr_utld_kl_m5 = @. l*r^(l-1)*uq_m5 + ψq_m5*r^l

            utld_kl_m0 = @. r^l*uq_m0
            utld_kl_m1 = @. r^l*uq_m1
            utld_kl_m2 = @. r^l*uq_m2
            utld_kl_m5 = @. r^l*uq_m5

            for i in 1:Nr
                for j in 1:Nr
                    #correlators[i,j,1]
                    buffers[i,j,1,id] += p.dk*(2.0*l+1.0)*(
                        dt_utld_kl_m0[i]*conj(dt_utld_kl_m0[j]) -
                        2.0*dt_utld_kl_m1[i]*conj(dt_utld_kl_m1[j]) +
                        2.0*dt_utld_kl_m2[i]*conj(dt_utld_kl_m2[j]) -
                        dt_utld_kl_m5[i]*conj(dt_utld_kl_m5[j]) )

                    #correlators[i,j,2]
                    buffers[i,j,2,id] += p.dk*(2.0*l+1.0)*(
                        dr_utld_kl_m0[i]*conj(dr_utld_kl_m0[j]) -
                        2.0*dr_utld_kl_m1[i]*conj(dr_utld_kl_m1[j]) +
                        2.0*dr_utld_kl_m2[i]*conj(dr_utld_kl_m2[j]) -
                        dr_utld_kl_m5[i]*conj(dr_utld_kl_m5[j]) )

                    #correlators[i,j,3]
                    buffers[i,j,3,id] += p.dk*(2.0*l+1.0)*(
                        utld_kl_m0[i]*conj(utld_kl_m0[j]) -
                        2.0*utld_kl_m1[i]*conj(utld_kl_m1[j]) +
                        2.0*utld_kl_m2[i]*conj(utld_kl_m2[j]) -
                        utld_kl_m5[i]*conj(utld_kl_m5[j]) )

                    #correlators[i,j,4]
                    buffers[i,j,4,id] += p.dk*(
                        πq_m0[i]*conj(πq_m0[j]) -
                        2.0*πq_m1[i]*conj(πq_m1[j]) +
                        2.0*πq_m2[i]*conj(πq_m2[j]) -
                        πq_m5[i]*conj(πq_m5[j]))

                    #correlators[i,j,5]
                    buffers[i,j,5,id] += p.dk*(
                        ψq_m0[i]*conj(ψq_m0[j]) -
                        2.0*ψq_m1[i]*conj(ψq_m1[j]) +
                        2.0*ψq_m2[i]*conj(ψq_m2[j]) -
                        ψq_m5[i]*conj(ψq_m5[j]))
                end # end loops in Nr
            end # end loops in Nr
        end # end loops in k
    end # end loops in l

    correlators[:,:,1] .= sum(buffers[:,:,1,:], dims=3)
    correlators[:,:,2] .= sum(buffers[:,:,2,:], dims=3)
    correlators[:,:,3] .= sum(buffers[:,:,3,:], dims=3)
    correlators[:,:,4] .= sum(buffers[:,:,4,:], dims=3)
    correlators[:,:,5] .= sum(buffers[:,:,5,:], dims=3)

    nothing
end

#end correlators
