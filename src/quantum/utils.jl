# used to save the various types of data
function save_all(p::Param, it::Int64, t::Float64,
                  v_classic::Array,
                  v_quantum::Array{ComplexF64},
                  r::Array, data_dir::String)

    if p.save_data == true && it % p.data_every == 0
        print("\rwriting data...")
        # don't save r grid if rmax is fixed
        if p.infalling_rmax == false
            ti = now()
            write_data(it, t, data_dir, v_classic)
            tf = now()
        else
            ti = now()
            write_data(it, t, data_dir, v_classic, r)
            tf = now()
        end
        print("\rwrote data in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
        println()
    end

    if p.save_quantum == true && it % p.quantum_every == 0
        print("\rwriting quantum...")
        # don't save r grid if rmax is fixed
        if p.infalling_rmax == false
            ti = now()
            write_quantum(it, t, data_dir, v_quantum)
            tf = now()
        else
            ti = now()
            write_quantum(it, t, data_dir, v_quantum, r)
            tf = now()
        end
        print("\rwrote quantum in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
        println()
    end

    if p.save_data_r0 == true && it % p.data_r0_every == 0
        print("\rwriting data_r0...")
        ti = now()
        write_data_r0(it, t, data_dir, v_classic, r)
        tf = now()
        print("\rwrote data_r0 in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
        println()
    end

    if p.save_quantum_r0 == true && it % p.quantum_r0_every == 0
        print("\rwriting quantum_r0...")
        ti = now()
        write_quantum_r0(it, t, data_dir, v_quantum, r)
        tf = now()
        print("\rwrote quantum_r0 in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
        println()
    end
 
    if p.save_bilinears == true && it % p.bilinears_every == 0
        bln = zeros(ComplexF64, (length(r), 5))
        # calculate the bilinears
        print("\rcalculating bilinears...")
        ti = now()
        bilinears(t, v_classic, v_quantum, p, r, bln)
        tf = now()
        print("\rcalculated bilinears in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
        println()
        print("\rwriting bilinears...")
        ti = now()
        # write bilinears
        write_bilinears(it, t, data_dir, bln)
        tf = now()
        print("\rwrote bilinears in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
        println()
    end

    if p.save_correlators == true && it % p.correlators_every == 0
        crlt = zeros(ComplexF64, (length(r), length(r), 5))
        # calculate correlators
        print("\rcalculating correlators...")
        ti = now()
        correlators(v_classic, v_quantum, p, r, crlt)
        tf = now()
        print("\rcalculated correlators in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
        println()
        print("\rwriting correlators...")
        ti = now()
        # write bilinears and correlators
        write_correlators(it, t, data_dir, crlt)
        tf = now()
        print("\rwrote correlators in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
        println()
    end

end

# function that prints output to inform the user while the code is running
# MULTIPLE DISPATCH

# take a real-valued grid function f and return
# the Maximum, of the Absolute value (at each r_i), abbreviated as "MA"
function MA(f::Vector{Float64})
    maximum(abs.(f))
end

# non-regularized version
function out(p::Param, speed::Float64, it::Int64, t::Float64,
             v_classic::Array{Float64, 2},
             v_quantum::Array{ComplexF64, 4}, # non-reg. no ghosts
             rr::Array, r_AH::Float64, data_dir::String)

    # uqmax is the real part of the highest k,l mode: maximum(abs.(real((v_quantum[:, end, end, 1])))
    if p.infalling_rmax == false
        @printf "         %9.2f | %9d   %9.3f |  %9.4g  |  %9.4g |   %9.4g  | %9.4g  |  %9.4g     %9.4g |\n" speed it t MA(v_classic[:,1]) MA(real.(v_quantum[:, end, end, 1])) MA(v_classic[:,13]) MA(v_classic[:,14]) minimum(v_classic[:,11]) maximum(v_classic[:,11])
    else
        @printf "         %9.2f | %9d   %9.3f |  %9.4g  |  %9.4g |   %9.4g  | %9.4g  |  %9.4g     %9.4g |   %9.8g  |   %9.3f  |\n" speed it t MA(v_classic[:,1]) MA(real.(v_quantum[:, end, end, 1])) MA(v_classic[:,13]) MA(v_classic[:,14]) minimum(v_classic[:,11]) maximum(v_classic[:,11]) r_AH rr[end]
    end

    save_all(p, it, t, v_classic, v_quantum, rr, data_dir)

end

# regularized version
function out(p::Param, speed::Float64, it::Int64, t::Float64,
             v_classic::Array{Float64, 2},
             v_quantum::Array{ComplexF64, 5}, # regularized; extra dim. for ghosts
             rr::Array, r_AH::Float64, data_dir::String)

    # uqmax is the real part of the highest k,l, and most massive: maximum(abs.(real((v_quantum[:, end, end, end, 1])))
    if p.infalling_rmax == false
        @printf "         %9.2f | %9d   %9.3f |  %9.4g  |  %9.4g |   %9.4g  | %9.4g  |  %9.4g     %9.4g |\n" speed it t MA(v_classic[:,1]) MA(real.(v_quantum[:, end, end, end, 1])) MA(v_classic[:,13]) MA(v_classic[:,14]) minimum(v_classic[:,11]) maximum(v_classic[:,11])
    else
        @printf "         %9.2f | %9d   %9.3f |  %9.4g  |  %9.4g |   %9.4g  | %9.4g  |  %9.4g     %9.4g |   %9.8g  |   %9.3f  |\n" speed it t MA(v_classic[:,1]) MA(real.(v_quantum[:, end, end, end, 1])) MA(v_classic[:,13]) MA(v_classic[:,14]) minimum(v_classic[:,11]) maximum(v_classic[:,11]) r_AH rr[end]
    end

    save_all(p, it, t, v_classic, v_quantum, rr, data_dir)

end

# function that checks if Ah is found, saves its radius and exit the code if p.AH == true
function AH_break(r_AH::Float64, p::Param, data_dir::String,
                  it::Int64, t::Float64,
                  v_classic::Array, v_quantum::Array,
                  rr::Array)
    if r_AH > 0 && p.AH == true
        println("AH found at r = ", r_AH)
        outfile = joinpath(data_dir, "r_AH.txt")
        writedlm(outfile, r_AH)
        write_checkpoint(it, t, data_dir, v_classic, v_quantum, rr)
        println("Exiting code")
        exit()
    end
end
