
# function that performs the time evolution
function run_classical(grid::Grid, p::Param)

    t1 = now() # for human readable total runtime
    tstart = time() # for checkpoint runtime
    last_checkpoint_time = 0.0

    #println("BLAS.get_num_threads() = ", BLAS.get_num_threads())
    println("Start run_classical with number of threads = ", Threads.nthreads())

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
    v_classic = zeros(Float64, ( length(rr), 18) )

    # search for checkpoints
    (its, all_filenames) = list_h5_files(data_dir, prefix="checkpoint_")
    if length(its) == 0 # no checkpoint
        # classic initial data
        v_classic = classical_ID(v_classic, sys, p)
        # if random is on
        if p.rand == true
            println("random true")
            v_classic .= v_classic .+ p.A_rand*randn( (length(rr), 18) )
        end
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
        v_classic = h5read(data_dir*"checkpoint_$(it_str).h5","v")
        # load the correct r-grid
        new_rmax = h5read(data_dir*"checkpoint_$(it_str).h5","r")[end]
        grid.r_max = new_rmax
        # pass the parameters of the system
        sys = System(grid)
        rr = sys.r
    end

    # iteration speed
    speed = 0.00
    println("-------------------------------------------------------------")
    if p.infalling_rmax == false
        println("speed (steps/hour) |   Iteration      Time |      Φ      |      Φ     |     |Θ|      |     |Zr|   |            α             |")
        println("                   |                       |     max     |     min    |     max      |     max    |    min           max     |")
    else
        println("speed (steps/hour) |   Iteration      Time |      Φ      |      Φ     |     |Θ|      |     |Zr|   |            α             |    r_AH      |     rmax    |")
        println("                   |                       |     max     |     min    |     max      |     max    |    min           max     |              |             |")
    end
    println("-------------------------------------------------------------")
    # find AH radius
    r_AH = find_AH(sys, v_classic[:,5], v_classic[:,4], v_classic[:,9])
    # print info and save data
    out(p, speed, it, t, v_classic, rr, r_AH, data_dir)
    # check if AH is found and exit code if this option is chosen; saves a checkpoint as well
    AH_break(r_AH, p, data_dir, it, t, v_classic, rr)

    #for AB3
    F0_classic = similar(v_classic)
    F1_classic = similar(v_classic)
    F2_classic = similar(v_classic)

    # for infalling rmax
    if p.infalling_rmax == true
        #for projection
        F0_classic_on_rr = similar(v_classic)
        F1_classic_on_rr = similar(v_classic)
        F2_classic_on_rr = similar(v_classic)
        r0 = similar(rr)
        r1 = similar(rr)
        r2 = similar(rr)
    end

    # from it=0
    # calculate and save the rhs with ID
    F!(t, F0_classic, v_classic, sys, p)

    # for infalling rmax
    if p.infalling_rmax == true
        r0 .= rr
    end

    #1st step
    it += 1
    t  += dt
    # mutates v
    ti = now()
    RK4(t, dt, v_classic, sys, p)
    tf = now()
    speed = round(3.6e+6/Dates.value(tf-ti),digits=2)
    # find AH radius
    r_AH = find_AH(sys, v_classic[:,5], v_classic[:,4], v_classic[:,9])
    # print info and save data
    out(p, speed, it, t, v_classic, rr, r_AH, data_dir)
    # check if AH is found and exit code if this option is chosen; saves a checkpoint as well
    AH_break(r_AH, p, data_dir, it, t, v_classic, rr)
    
    #for it=1, calculate and save rhs
    F!(t, F1_classic, v_classic, sys, p)

    # for infalling rmax
    if p.infalling_rmax == true
        r1 .= rr
    end

    # 2nd step
    it += 1
    t  += dt
    # mutates v
    ti = now()
    RK4(t, dt, v_classic, sys, p)
    tf = now()
    speed = round(3.6e+6/Dates.value(tf-ti),digits=2)
    # find AH radius
    r_AH = find_AH(sys, v_classic[:,5], v_classic[:,4], v_classic[:,9])
    # print info and save data
    out(p, speed, it, t, v_classic, rr, r_AH, data_dir)
    # check if AH is found and exit code if this option is chosen; saves a checkpoint as well
    AH_break(r_AH, p, data_dir, it, t, v_classic, rr)

    #for it=2, calculate and save rhs
    F!(t, F2_classic, v_classic, sys, p)

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

            AB3!(v_classic, F2_classic, F1_classic, F0_classic, dt)

            """
            roll-over the F0, F1, and F2, that is mutate them
            for the next AB3
            """
            F0_classic .= F1_classic
            F1_classic .= F2_classic
            F!(t, F2_classic, v_classic, sys, p)

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
            @threads for i in 1:length(v_classic[1,:])
                v_classic[:,i] .= project_on_smaller_grid(r2, rr, v_classic[:,i])
                F2_classic_on_rr[:,i] .= project_on_smaller_grid(r2, rr, F2_classic[:,i])
                F1_classic_on_rr[:,i] .= project_on_smaller_grid(r1, rr, F1_classic[:,i])
                F0_classic_on_rr[:,i] .= project_on_smaller_grid(r0, rr, F0_classic[:,i])
            end

            AB3!(v_classic, F2_classic_on_rr, F1_classic_on_rr, F0_classic_on_rr, dt)

            """
            roll-over the F0, F1, and F2, that is mutate them
            for the next AB3
            """
            F0_classic .= F1_classic
            F1_classic .= F2_classic
            F!(t, F2_classic, v_classic, sys, p)
            r0 .= r1
            r1 .= r2
            r2 .= rr

            # adapt dt to new r-grid
            dt   = p.cfl* sys.hr

        end # end infalling rmax check

        tf = now()
        speed = round(3.6e+6/Dates.value(tf-ti),digits=2)
        # find AH radius
        r_AH = find_AH(sys, v_classic[:,5], v_classic[:,4], v_classic[:,9])
        # print info and save data
        out(p, speed, it, t, v_classic, rr, r_AH, data_dir)
        # check if AH is found and exit code if this option is chosen; saves a checkpoint as well
        AH_break(r_AH, p, data_dir, it, t, v_classic, rr)

        # for checkpoint
        runtime = time() - tstart # in seconds; Float64
        if p.save_checkpoint == true && runtime >= last_checkpoint_time + p.checkpoint_every*3600
            last_checkpoint_time = runtime
            print("\rwriting checkpoint...")
            ti = now()
            write_checkpoint(it, t, data_dir, v_classic, rr)
            tf = now()
            print("\rwrote checkpoint in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
            println()
        end

        # check for termination folder
        if isdir(data_dir*"/stop")
            println("Termination instruction found. I'm wrapping up...")
            print("\rwriting checkpoint...")
            ti = now()
            write_checkpoint(it, t, data_dir, v_classic, rr)
            tf = now()
            print("\rwrote checkpoint in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
            exit()
        end

    end
    #end while

    # write final checkpoint
    write_checkpoint(it, t, data_dir, v_classic, rr)

    # print AH radius
    println("r_AH = $(r_AH)")
    # if there is an AH, save its position in txt; needed for bisection if AH_break in false
    if r_AH > 0
        outfile = joinpath(data_dir, "r_AH.txt")
        writedlm(outfile, r_AH)
    end
    println()

    println("-------------------------------------------------------------")
    println("Done.")

    t2 = now()
    println("total run time = ", Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(t2)-Dates.DateTime(t1))))

end
#end evolution function
