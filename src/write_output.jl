# MULTIPLE DISPATCH
# start write_output

# saves the state vector v, and the timestep t
function write_data(it::Int, t, data_dir::String,
                    v::Array)
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "data_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "v", v)
        attrs(file)["time"] = t
    end
    nothing
end
# saves r grid as well; useful when rmax decreases
function write_data(it::Int, t, data_dir::String,
                    v::Array, r::Array)
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "data_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "v", v)
        write(file, "r", r)
        attrs(file)["time"] = t
    end
    nothing
end

# saves the quantum state vector
# works for both non-regularized and regularized version
function write_quantum(it::Int, t, data_dir::String,
                       v_quantum::Array)
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "quantum_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "v", v_quantum)
        attrs(file)["time"] = t
    end
    nothing
end
# with r grid
function write_quantum(it::Int, t, data_dir::String,
                       v_quantum::Array, r::Array)
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "quantum_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "v", v_quantum)
        write(file, "r", r)
        attrs(file)["time"] = t
    end
    nothing
end

# save only the 5 points around r=0 for all elements of the state vector v,
# the respective r-grid, that is r = (-2h, -h, 0, h, 2h)
# s.t. we can take derivatives in postprocessing, as well as timestep t
function write_data_r0(it::Int, t, data_dir::String,
                       v::Array, r::Array)
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "r0data_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "v" , v[1:5,:])
        write(file, "r" , r[1:5])
        attrs(file)["time"] = t
    end
    nothing
end
# quantum; MULTIPLE DISPATCH
# non-regularized
function write_quantum_r0(it::Int, t, data_dir::String,
                          v_quantum::Array{ComplexF64, 4}, r::Array)
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "r0quantum_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "v" , v_quantum[1:5,:,:,:])
        write(file, "r" , r[1:5])
        attrs(file)["time"] = t
    end
    nothing
end
# regularized
function write_quantum_r0(it::Int, t, data_dir::String,
                          v_quantum::Array{ComplexF64, 5}, r::Array)
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "r0quantum_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "v", v_quantum[1:5,:,:,:,:])
        write(file, "r", r[1:5])
        attrs(file)["time"] = t
    end
    nothing
end

# write bilinears (1-point function; an 1D array for each t)
function write_bilinears(it::Int, t, data_dir::String,
                         bilinears::Array)
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "bilinears_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "bilinears" ,  bilinears)
        attrs(file)["time"] = t
    end
    nothing
end

# write correlators (2-point function; a 2D array for each t)
function write_correlators(it::Int, t, data_dir::String,
                           correlators::Array)
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "correlators_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "correlators" ,  correlators)
        attrs(file)["time"] = t
    end
    nothing
end

# MULTIPLE DISPATCH

# write checkpoint of the full state vector v
# for classical v = v_classical
function write_checkpoint(it::Int, t::Float64, data_dir::String,
                          v::Array{Float64}, r::Array{Float64})
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "checkpoint_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "v" ,  v)
        write(file, "r" ,  r)
        attrs(file)["time"] = t
    end
    nothing
end

# for quantum it works for both regularized and non-regularized versions
# saves both classical anc quantum state vectors (separately)
function write_checkpoint(it::Int, t::Float64, data_dir::String,
                          v_classic::Array,
                          v_quantum::Array{ComplexF64},
                          r::Array{Float64})
    it_str  = lpad(it, 4, "0")
    outfile = joinpath(data_dir, "checkpoint_$(it_str).h5")
    h5open(outfile, "w") do file
        write(file, "v_classic" ,  v_classic)
        write(file, "v_quantum" ,  v_quantum)
        write(file, "r" ,  r)
        attrs(file)["time"] = t
    end
    nothing
end

# end write_output

