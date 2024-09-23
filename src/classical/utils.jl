# function that prints output to inform the user while the code is running
function out(p::Param, speed::Float64, it::Int64, t::Float64, v_classic::Array,
             rr::Array, r_AH::Float64, data_dir::String)

    if p.infalling_rmax == false
        @printf "         %9.2f | %9d   %9.3f |  %9.4g  |  %9.4g |   %9.4g  | %9.4g  |  %9.4g     %9.4g |\n" speed it t maximum(v_classic[:,1]) minimum(v_classic[:,1]) maximum(abs.(v_classic[:,13])) maximum(abs.(v_classic[:,14])) minimum(v_classic[:,11]) maximum(v_classic[:,11])
    else
        @printf "         %9.2f | %9d   %9.3f |  %9.4g  |  %9.4g |   %9.4g  | %9.4g  |  %9.4g     %9.4g |   %9.8g  |   %9.3f  |\n" speed it t maximum(v_classic[:,1]) minimum(v_classic[:,1]) maximum(abs.(v_classic[:,13])) maximum(abs.(v_classic[:,14])) minimum(v_classic[:,11]) maximum(v_classic[:,11]) r_AH rr[end]
    end

    if p.save_data == true && it % p.data_every == 0
        print("\rwriting data...")
        ti = now()
        write_data(it, t, data_dir, v_classic, rr)
        tf = now()
        print("\rwrote data in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
        println()
    end

    if p.save_data_r0 == true && it % p.data_r0_every == 0
        print("\rwriting data_r0...")
        ti = now()
        write_data_r0(it, t, data_dir, v_classic, rr)
        tf = now()
        print("\rwrote data_r0 in $(Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(tf)-Dates.DateTime(ti))))")
        println()
    end

end

# function that checks if an AH is found
# saves its radius and exit the code if p.AH == true
function AH_break(r_AH::Float64, p::Param, data_dir::String,
                  it::Int64, t::Float64, v_classic::Array, rr::Array)
    if r_AH > 0 && p.AH == true
        println("AH found at r = ", r_AH)
        outfile = joinpath(data_dir, "r_AH.txt")
        writedlm(outfile, r_AH)
        write_checkpoint(it, t, data_dir, v_classic, rr)
        println("Exiting code")
        exit()
    end
end
