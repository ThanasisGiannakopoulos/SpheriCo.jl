
# function that performs the time evolution
function run_quantum(grid::Grid, p::Param)

    t1 = now() # for human readable total runtime
    tstart = time() # for checkpoint runtime
    last_checkpoint_time = 0.0

    println("Start run_quantum with number of threads = ", Threads.nthreads())
    println("CC = ", p.CC)
    println("mlist = ", p.mlist)

    # pass the parameters of the system
    sys = System(grid)
    rr = sys.r
    Nr = grid.Nr

    # create the folders where data are saved
    data_dir = p.out_dir
    """
    A termination folder with the name "stop" can be created during the 
    simulation if one wants to stop the simulations (a checkpoint is saved) 
    """
    # remove termination folder if it exists from precious run
    if isdir(data_dir*"/stop")
        rm(data_dir*"/stop")
    end

    # timestep. careful with CFL condition
    dt   = p.cfl* sys.hr

    # initiate the variables
    # the state vector is:
    # v_classic = [Φ, Π, Ψ, A, B, DB, Utld, K, KB, λ, α, Dα, Θ, Zr, f, g, U, V]^T
    # allow for complex value to check if the backreaction is done correctly
    # should be real-valued (to machine precision?)
    v_classic = zeros(Float64, ( length(rr), 18) )
    # quantum state vector:
    if p.PV_reg==false
        # non-regularized version: one for each k in [1,kmax], l in [0,lmax], 1 + 2 reduction vars 
        v_quantum = zeros(ComplexF64, (length(rr), Int(p.lmax+1), Int(p.kmax), 3))
    else
        # regularized version: one for each k in [1,kmax], l in [0,lmax], and m in p.mlist
        v_quantum = zeros(ComplexF64, (length(rr), Int(p.lmax+1), Int(p.kmax), length(p.mlist), 3))
    end

    # search for checkpoints
    (its, all_filenames) = list_h5_files(data_dir, prefix="checkpoint_")
    if length(its) == 0 # no checkpoint
        # classic initial data
        v_classic = classical_ID(v_classic, sys, p)
        # quantum initial data
        v_quantum = quantum_ID(v_quantum, sys, p)
        # write the coordinates.
        h5write(data_dir*"/r.h5", "r", sys.r)
        # first iteration
        it = 0
        t  = 0.0
    else # load from checkpoint
        println("Loading from checkpoint...")
        it        = its[end]
        it_str    = lpad(it, 4, "0")
        t         = h5readattr(data_dir*"checkpoint_$(it_str).h5", "/")["time"]
        v_classic = h5read(data_dir*"checkpoint_$(it_str).h5","v_classic")
        v_quantum = h5read(data_dir*"checkpoint_$(it_str).h5","v_quantum")
        # load the correct r-grid
        new_rmax = h5read(data_dir*"checkpoint_$(it_str).h5","r")[end]
        grid.r_max = new_rmax
        # pass the parameters of the system
        sys = System(grid)
        rr = sys.r
    end

    println("|uqmax| = ", maximum(abs.(real.(v_quantum[:, end, end, end, 1]))))

    # iteration speed
    speed = 0.00
    println("-------------------------------------------------------------")
    if p.infalling_rmax == false
        println("speed (steps/hour) |   Iteration      Time |     |Φ|     |   |uqmax|  |     |Θ|      |     |Zr|   |            α             |")
        println("                   |                       |     max     |     max    |     max      |     max    |    min           max     |")
    else
        println("speed (steps/hour) |   Iteration      Time |     |Φ|     |   |uqmax|  |     |Θ|      |     |Zr|   |            α             |    r_AH      |     rmax    |")
        println("                   |                       |     max     |     max   |     max      |     max    |    min           max     |              |             |")
    end
    println("-------------------------------------------------------------")
    # find AH radius
    r_AH = find_AH(sys, v_classic[:,5], v_classic[:,4], v_classic[:,9])
    # print info and save data
    out(p, speed, it, t, v_classic, v_quantum, rr, r_AH, data_dir)
    # check if AH is found and exit code if this option is chosen; saves a checkpoint as well
    AH_break(r_AH, p, data_dir, it, t, v_classic, v_quantum, rr)

    #for AB3
    F0_classic = similar(v_classic)
    F1_classic = similar(v_classic)
    F2_classic = similar(v_classic)
    F0_quantum = similar(v_quantum)
    F1_quantum = similar(v_quantum)
    F2_quantum = similar(v_quantum)

    # for infalling rmax
    if p.infalling_rmax == true
        #for projection
        F0_classic_on_rr = similar(v_classic)
        F1_classic_on_rr = similar(v_classic)
        F2_classic_on_rr = similar(v_classic)
        F0_quantum_on_rr = similar(v_quantum)
        F1_quantum_on_rr = similar(v_quantum)
        F2_quantum_on_rr = similar(v_quantum)
        r0 = similar(rr)
        r1 = similar(rr)
        r2 = similar(rr)
    end

    # from it=0
    # calculate and save the rhs with ID
    F!(t, F0_classic, v_classic, F0_quantum, v_quantum, sys, p)

    # for infalling rmax
    if p.infalling_rmax == true
        r0 .= rr
    end

    #1st step
    it += 1
    t  += dt
    # mutates v
    ti = now()
    RK4(t, dt, v_classic, v_quantum, sys, p)
    tf = now()
    speed = round(3.6e+6/Dates.value(tf-ti),digits=2)
    # find AH radius
    r_AH = find_AH(sys, v_classic[:,5], v_classic[:,4], v_classic[:,9])
    # print info and save data
    out(p, speed, it, t, v_classic, v_quantum, rr, r_AH, data_dir)
    # check if AH is found and exit code if this option is chosen; saves a checkpoint as well
    AH_break(r_AH, p, data_dir, it, t, v_classic, v_quantum, rr)

    #for it=1, calculate and save rhs
    F!(t, F1_classic, v_classic, F1_quantum, v_quantum, sys, p)

    # for infalling rmax
    if p.infalling_rmax == true
        r1 .= rr
    end

    # 2nd step
    it += 1
    t  += dt
    # mutates v
    ti = now()
    RK4(t, dt, v_classic, v_quantum, sys, p)
    tf = now()
    speed = round(3.6e+6/Dates.value(tf-ti),digits=2)
    # find AH radius
    r_AH = find_AH(sys, v_classic[:,5], v_classic[:,4], v_classic[:,9])
    # print info and save data
    out(p, speed, it, t, v_classic, v_quantum, rr, r_AH, data_dir)
    # check if AH is found and exit code if this option is chosen; saves a checkpoint as well
    AH_break(r_AH, p, data_dir, it, t, v_classic, v_quantum, rr)

    #for it=2, calculate and save rhs
    F!(t, F2_classic, v_classic, F2_quantum, v_quantum, sys, p)

    # for infalling rmax
    if p.infalling_rmax == true
        r2 .= rr
    end

    # start time evolution with AB3
    while t <= p.t_max
        it += 1
        t  += dt

        ti = now()

        if p.infalling_rmax == false

            AB3!(v_classic, F2_classic, F1_classic, F0_classic,
                 v_quantum, F2_quantum, F1_quantum, F0_quantum,
                 dt)

            """
            roll-over the F0, F1, and F2, that is mutate them
            for the next AB3
            """
            F0_classic .= F1_classic
            F1_classic .= F2_classic
            F0_quantum .= F1_quantum
            F1_quantum .= F2_quantum
            F!(t, F2_classic, v_classic, F2_quantum, v_quantum, sys, p)

        else # for infalling rmax

            # update rmax in params such that it falls at the speed of light
            # drops the r grid on the timestep you land after AB3, by t
            new_rmax = rr[end] - dt
            grid.r_max = new_rmax
            # pass the parameters of the system
            sys = System(grid)
            rr = sys.r

            # projections on new smaller grid using:
            # project_on_smaller_grid!(r_old::Array, r_new::Array, f::Array)
            # classical
            @threads for i in 1:length(v_classic[1,:])
                v_classic[:,i] .= project_on_smaller_grid(r2, rr, v_classic[:,i])
                F2_classic_on_rr[:,i] .= project_on_smaller_grid(r2, rr, F2_classic[:,i])
                F1_classic_on_rr[:,i] .= project_on_smaller_grid(r1, rr, F1_classic[:,i])
                F0_classic_on_rr[:,i] .= project_on_smaller_grid(r0, rr, F0_classic[:,i])
            end
            # quantum
            if p.PV_reg==false
                @threads for l in 0:p.lmax
                    for k in 1:p.kmax
                        for i in 1:3 # for the 3 quantum fields uq_klm, ψq_klm, πq_klm
                            v_quantum[:,Int(l+1), Int(k), i] .= project_on_smaller_grid(
                                r2, rr, v_quantum[:,Int(l+1), Int(k), i])
                            F2_quantum_on_rr[:,Int(l+1), Int(k), i] .=
                                project_on_smaller_grid(
                                    r2, rr, F2_quantum[:,Int(l+1), Int(k), i])
                            F1_quantum_on_rr[:,Int(l+1), Int(k), i] .=
                                project_on_smaller_grid(
                                    r1, rr, F1_quantum[:,Int(l+1), Int(k), i])
                            F0_quantum_on_rr[:,Int(l+1), Int(k), i] .=
                                project_on_smaller_grid(
                                    r0, rr, F0_quantum[:,Int(l+1), Int(k), i])
                        end #  end i loop
                    end # end k loop
                end # end l loop
            else
                @threads for l in 0:p.lmax
                    for k in 1:p.kmax
                        for m in 1:length(p.mlist)
                            for i in 1:3 # for the 3 quantum fields uq_klm, ψq_klm, πq_klm
                                v_quantum[:,Int(l+1), Int(k), m, i] .= project_on_smaller_grid(
                                    r2, rr, v_quantum[:,Int(l+1), Int(k), m, i])
                                F2_quantum_on_rr[:,Int(l+1), Int(k), m, i] .=
                                    project_on_smaller_grid(
                                        r2, rr, F2_quantum[:,Int(l+1), Int(k), m, i])
                                F1_quantum_on_rr[:,Int(l+1), Int(k), m, i] .=
                                    project_on_smaller_grid(
                                        r1, rr, F1_quantum[:,Int(l+1), Int(k), m, i])
                                F0_quantum_on_rr[:,Int(l+1), Int(k), m, i] .=
                                    project_on_smaller_grid(
                                        r0, rr, F0_quantum[:,Int(l+1), Int(k), m, i])
                            end #  end i loop
                        end # end m loop
                    end # end k loop
                end # end l loop
            end # end if

            AB3!(v_classic, F2_classic_on_rr, F1_classic_on_rr, F0_classic_on_rr,
                 v_quantum, F2_quantum_on_rr, F1_quantum_on_rr, F0_quantum_on_rr,
                 dt)
            """
            roll-over the F0, F1, and F2, that is mutate them
            for the next AB3
            """
            F0_classic .= F1_classic
            F1_classic .= F2_classic
            F0_quantum .= F1_quantum
            F1_quantum .= F2_quantum
            F!(t, F2_classic, v_classic, F2_quantum, v_quantum, sys, p)

            # adapt dt to new r-grid
            dt   = p.cfl* sys.hr

        end # end infalling rmax check

        tf = now()
        speed = round(3.6e+6/Dates.value(tf-ti),digits=2)
        # find AH radius
        r_AH = find_AH(sys, v_classic[:,5], v_classic[:,4], v_classic[:,9])
        # print info and save data
        out(p, speed, it, t, v_classic, v_quantum, rr, r_AH, data_dir)
        # check if AH is found and exit code if this option is chosen; saves a checkpoint as well
        AH_break(r_AH, p, data_dir, it, t, v_classic, v_quantum, rr)

        # for checkpoint
        runtime = time() - tstart # in seconds; Float64
        if p.save_checkpoint == true && runtime >= last_checkpoint_time + p.checkpoint_every*3600
            last_checkpoint_time = runtime
            print("\rwriting checkpoint...")
            ti = now()
            write_checkpoint(it, t, data_dir, v_classic, v_quantum, rr)
            tf = now()
            print("\rwrote checkpoint in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
            println()
        end

        # check for termination folder
        if isdir(data_dir*"/stop")
            println("Termination instruction found. I'm wrapping up...")
            print("\rwriting checkpoint...")
            ti = now()
            write_checkpoint(it, t, data_dir, v_classic, v_quantum, rr)
            tf = now()
            print("\rwrote checkpoint in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
            exit()
        end

    end # end while

    # write final checkpoint
    write_checkpoint(it, t, data_dir, v_classic, v_quantum, rr)

    # print AH radius
    println("r_AH = $(r_AH)")
    println()

    println("-------------------------------------------------------------")
    println("Done.")

    t2 = now()
    println("total run time = ", Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(t2)-Dates.DateTime(t1))))

end
# end evolution function
