using HDF5
using LaTeXStrings
using Plots#; pythonplot()
using DelimitedFiles
using SpheriCo

#using ProgressBars
import Base.Threads.@threads

println("running with number of threads = ", Threads.nthreads())

# give the directory where the data from all the runs are saved
dir = "../examples/quantum_runs/"
par = "rmax30.0_cfl0.0625_sigma0.0_amp0.0_width2.0_rc5.0_infalling_rmax_false_kmax10.0_lmax60.0_reg_true_backreact_false_mPV1.0"
your_dir = dir*par

# set manually here; make sure it's correct
mPV = 1.0;
dk = π/30.0;
mlist = [0.0, 1.0*mPV, sqrt(3.0)*mPV, 2.0*mPV];

# make the output dir
out_dir = "./convergence_plots/standing_waves/"*par*"/norms/"
if ispath(out_dir)==false
    mkpath(out_dir)
end

# analytical solutions (standing waves)
function uq_analytic_t(r, t, k, m, l)
    ω = sqrt(k^2 + m^2)
    SpheriCo.quantum.uq_ID(r, k, m, l)*exp(-im*ω*t)
end
function ψq_analytic_t(r, t, k, m, l)
    ω = sqrt(k^2 + m^2)
    SpheriCo.quantum.ψq_ID(r, k, m, l)*exp(-im*ω*t)
end
function πq_analytic_t(r, t, k, m, l)
    ω = sqrt(k^2 + m^2)
    SpheriCo.quantum.πq_ID(r, k, m, l)*exp(-im*ω*t)
end

# set manually to that of the coarse simulation
D = 2
Nr = 128*2^D + 3 # the overal coarse graining

# load the r grid
r1 = h5read(your_dir*"/data_$(Nr)/r.h5","r")
dr1 = r1[2] - r1[1]

r2 = h5read(your_dir*"/data_$((Nr-3)*2 +3)/r.h5","r")
dr2 = r2[2] - r2[1]

r3 = h5read(your_dir*"/data_$((Nr-3)*4 +3)/r.h5","r")
dr3 = r3[2] - r3[1]

@assert r1[3:end] ≈ r2[3:2:end]
@assert r1[3:end] ≈ r3[3:4:end]
@assert dr1 ≈ 2.0*dr2
@assert dr1 ≈ 4.0*dr3

# list all available iterations (and corresponding files)
(its_quantum1, all_filenames_quantum1) = list_h5_files(your_dir*"/data_$(Nr)", prefix="quantum_");
(its_quantum2, all_filenames_quantum2) = list_h5_files(your_dir*"/data_$((Nr-3)*2 +3)", prefix="quantum_");
(its_quantum3, all_filenames_quantum3) = list_h5_files(your_dir*"/data_$((Nr-3)*4 +3)", prefix="quantum_");

@assert length(its_quantum1) ≈ length(its_quantum2)
@assert length(its_quantum1) ≈ length(its_quantum3)

# get the saved (common) timesteps in a list
println("getting t_list")
t_list = zeros(length(its_quantum1))
for i in 1:length(its_quantum1)
    it_quantum1 = its_quantum1[i]
    it_str_quantum1  = lpad(it_quantum1, 4, "0")  
    it_quantum2 = its_quantum2[i]
    it_str_quantum2  = lpad(it_quantum2, 4, "0")
    it_quantum3 = its_quantum3[i]
    it_str_quantum3  = lpad(it_quantum3, 4, "0")

    t1 =  h5readattr(your_dir*"/data_$(Nr)/quantum_$(it_str_quantum1).h5", "/")["time"]
    t2 =  h5readattr(your_dir*"/data_$((Nr-3)*2 +3)/quantum_$(it_str_quantum2).h5", "/")["time"]
    t3 =  h5readattr(your_dir*"/data_$((Nr-3)*4 +3)/quantum_$(it_str_quantum3).h5", "/")["time"]

    @assert t1 ≈ t2
    @assert t1 ≈ t3

    t_list[i] = t1
end

# save t_list
outfile = joinpath(out_dir, "t_list.h5")
h5open(outfile, "w") do file
    write(file, "t_list", t_list)
end

# find ri for the first causally disconnected r point
# (towards the outer boundary) from the boundary
function get_ri_causal(rcausal)
    ri = length(r1)
    while r1[ri] >= rcausal
        ri = ri - 1
    end
    return ri
end

ri = get_ri_causal(r1[end] - t_list[end])
println("rcausal = ", r1[end] - t_list[end])
println("ri = ", ri)
@assert r1[3:ri] ≈ r2[3:2:ri*2-3]
@assert r1[3:ri] ≈ r3[3:4:ri*4-9]
    
# build the norms for specific l and m
# norm2 has size (number_of_different_resols, number_of_saved_timesteps)
function norm2_lm(norm2, l, m, ri_max)
    ri1 = ri_max
    ri2 = ri_max*2-3
    ri3 = ri_max*4-9
    
    # buffers for the sums in multithreading with
    # format (number_of_different_resols, number_of_saved_timesteps, number_of_threads)
    buffers = zeros(Float64, (3, length(its_quantum1), Threads.nthreads()))

    # to check the progress
    i_checked = 0.0
    for i in 1:length(its_quantum1)
        i_checked += 1.0
        print("\rnorm2_lm: $(round(i_checked/(length(its_quantum1))*100.0, digits=2)) % ")

        it_quantum1 = its_quantum1[i]
        it_str_quantum1 = lpad(it_quantum1, 4, "0")
        it_quantum2 = its_quantum2[i]
        it_str_quantum2 = lpad(it_quantum2, 4, "0")
        it_quantum3 = its_quantum3[i]
        it_str_quantum3 = lpad(it_quantum3, 4, "0")

        v1 = h5read(your_dir*"/data_$(Nr)/quantum_$(it_str_quantum1).h5","v")
        v2 = h5read(your_dir*"/data_$((Nr-3)*2+3)/quantum_$(it_str_quantum2).h5","v")
        v3 = h5read(your_dir*"/data_$((Nr-3)*4+3)/quantum_$(it_str_quantum3).h5","v")

        t1 = h5readattr(your_dir*"/data_$(Nr)/quantum_$(it_str_quantum1).h5","/")["time"]
        t2 = h5readattr(your_dir*"/data_$((Nr-3)*2+3)/quantum_$(it_str_quantum2).h5","/")["time"]
        t3 = h5readattr(your_dir*"/data_$((Nr-3)*4+3)/quantum_$(it_str_quantum3).h5","/")["time"]

        @assert t1 ≈ t2
        @assert t1 ≈ t3
        
        @threads for k in 1:length(v1[1,1,:,1,1])
            id = Threads.threadid()
            # uq diff against analytic
            uq1_diff = v1[3:ri1, Int(l)+1, Int(k), m, 1] .-
                uq_analytic_t.(r1[3:ri1], t1, dk*k, mlist[m], Float64(l))
            uq2_diff = v2[3:ri2, Int(l)+1, Int(k), m, 1] .-
                uq_analytic_t.(r2[3:ri2], t2, dk*k, mlist[m], Float64(l))
            uq3_diff = v3[3:ri3, Int(l)+1, Int(k), m, 1] .-
                uq_analytic_t.(r3[3:ri3], t3, dk*k, mlist[m], Float64(l))
            # ψq diff against analytic
            ψq1_diff = v1[3:ri1, Int(l)+1, Int(k), m, 2] .-
                ψq_analytic_t.(r1[3:ri1], t1, dk*k, mlist[m], Float64(l))
            ψq2_diff = v2[3:ri2, Int(l)+1, Int(k), m, 2] .-
                ψq_analytic_t.(r2[3:ri2], t2, dk*k, mlist[m], Float64(l))
            ψq3_diff = v3[3:ri3, Int(l)+1, Int(k), m, 2] .-
                ψq_analytic_t.(r3[3:ri3], t3, dk*k, mlist[m], Float64(l))
            # πq diff against analytic
            πq1_diff = v1[3:ri1, Int(l)+1, Int(k), m, 3] .-
                πq_analytic_t.(r1[3:ri1], t1, dk*k, mlist[m], Float64(l))
            πq2_diff = v2[3:ri2, Int(l)+1, Int(k), m, 3] .-
                πq_analytic_t.(r2[3:ri2], t2, dk*k, mlist[m], Float64(l))
            πq3_diff = v3[3:ri3, Int(l)+1, Int(k), m, 3] .-
                πq_analytic_t.(r3[3:ri3], t3, dk*k, mlist[m], Float64(l))
            # norm2 = Sum(uq_diff^2 + ψq_diff^2 + πq_diff^2)*r^l*dr
            buffers[1, i, id] += sum( dr1* r1[3:ri1].^l .*(uq1_diff.*conj.(uq1_diff) .+ ψq1_diff.*conj.(ψq1_diff) .+ πq1_diff.*conj.(πq1_diff)) )
            buffers[2, i, id] += sum( dr2* r2[3:ri2].^l .*(uq2_diff.*conj.(uq2_diff) .+ ψq2_diff.*conj.(ψq2_diff) .+ πq2_diff.*conj.(πq2_diff)) )
            buffers[3, i, id] += sum( dr3* r3[3:ri3].^l .*(uq3_diff.*conj.(uq3_diff) .+ ψq3_diff.*conj.(ψq3_diff) .+ πq3_diff.*conj.(πq3_diff)) )

        end # k loops
    end # i loop

    norm2 .= sum(buffers, dims=3)
    nothing
end

# build the norms for all l and m
# norm2 has size (number_of_different_resols, number_of_saved_timesteps)
function norm2(norm2, ri_max)
    ri1 = ri_max
    ri2 = ri_max*2-3
    ri3 = ri_max*4-9
    
    # buffers for the sums in multithreading with
    # format (number_of_different_resols, number_of_saved_timesteps, number_of_threads)
    buffers = zeros(Float64, (3, length(its_quantum1), Threads.nthreads()))

    # to check the progress
    i_checked = 0.0
    for i in 1:length(its_quantum1)
        i_checked += 1.0
        print("\rnorm2: $(round(i_checked/(length(its_quantum1))*100.0, digits=2)) % ")

        it_quantum1 = its_quantum1[i]
        it_str_quantum1 = lpad(it_quantum1, 4, "0")
        it_quantum2 = its_quantum2[i]
        it_str_quantum2 = lpad(it_quantum2, 4, "0")
        it_quantum3 = its_quantum3[i]
        it_str_quantum3 = lpad(it_quantum3, 4, "0")

        v1 = h5read(your_dir*"/data_$(Nr)/quantum_$(it_str_quantum1).h5","v")
        v2 = h5read(your_dir*"/data_$((Nr-3)*2+3)/quantum_$(it_str_quantum2).h5","v")
        v3 = h5read(your_dir*"/data_$((Nr-3)*4+3)/quantum_$(it_str_quantum3).h5","v")

        t1 = h5readattr(your_dir*"/data_$(Nr)/quantum_$(it_str_quantum1).h5","/")["time"]
        t2 = h5readattr(your_dir*"/data_$((Nr-3)*2+3)/quantum_$(it_str_quantum2).h5","/")["time"]
        t3 = h5readattr(your_dir*"/data_$((Nr-3)*4+3)/quantum_$(it_str_quantum3).h5","/")["time"]

        @assert t1 ≈ t2
        @assert t1 ≈ t3
        
        @threads for l in 0:length(v1[1,:,1,1,1])-1
            id = Threads.threadid()
            for k in 1:length(v1[1,1,:,1,1])
                for m in 1:length(mlist)
                    # uq diff against analytic
                    uq1_diff = v1[3:ri1, Int(l)+1, Int(k), m, 1] .-
                        uq_analytic_t.(r1[3:ri1], t1, dk*k, mlist[m], Float64(l))
                    uq2_diff = v2[3:ri2, Int(l)+1, Int(k), m, 1] .-
                        uq_analytic_t.(r2[3:ri2], t2, dk*k, mlist[m], Float64(l))
                    uq3_diff = v3[3:ri3, Int(l)+1, Int(k), m, 1] .-
                        uq_analytic_t.(r3[3:ri3], t3, dk*k, mlist[m], Float64(l))
                    # ψq diff against analytic
                    ψq1_diff = v1[3:ri1, Int(l)+1, Int(k), m, 2] .-
                        ψq_analytic_t.(r1[3:ri1], t1, dk*k, mlist[m], Float64(l))
                    ψq2_diff = v2[3:ri2, Int(l)+1, Int(k), m, 2] .-
                        ψq_analytic_t.(r2[3:ri2], t2, dk*k, mlist[m], Float64(l))
                    ψq3_diff = v3[3:ri3, Int(l)+1, Int(k), m, 2] .-
                        ψq_analytic_t.(r3[3:ri3], t3, dk*k, mlist[m], Float64(l))
                    # πq diff against analytic
                    πq1_diff = v1[3:ri1, Int(l)+1, Int(k), m, 3] .-
                        πq_analytic_t.(r1[3:ri1], t1, dk*k, mlist[m], Float64(l))
                    πq2_diff = v2[3:ri2, Int(l)+1, Int(k), m, 3] .-
                        πq_analytic_t.(r2[3:ri2], t2, dk*k, mlist[m], Float64(l))
                    πq3_diff = v3[3:ri3, Int(l)+1, Int(k), m, 3] .-
                        πq_analytic_t.(r3[3:ri3], t3, dk*k, mlist[m], Float64(l))
                    # norm2 = Sum(uq_diff^2 + ψq_diff^2 + πq_diff^2)*r^l*dr
                    buffers[1, i, id] += sum( dr1* r1[3:ri1].^l .*(uq1_diff.*conj.(uq1_diff) .+ 
                        ψq1_diff.*conj.(ψq1_diff) .+ πq1_diff.*conj.(πq1_diff)) )
                    buffers[2, i, id] += sum( dr2* r2[3:ri2].^l .*(uq2_diff.*conj.(uq2_diff) .+ 
                        ψq2_diff.*conj.(ψq2_diff) .+ πq2_diff.*conj.(πq2_diff)) )
                    buffers[3, i, id] += sum( dr3* r3[3:ri3].^l .*(uq3_diff.*conj.(uq3_diff) .+ 
                        ψq3_diff.*conj.(ψq3_diff) .+ πq3_diff.*conj.(πq3_diff)) )
                end # m loops
            end # k loops
        end # l loops
    end # i loops

    norm2 .= sum(buffers, dims=3)
    nothing
end

# calculate and save norms

# l=60, for higher mass ghost (fastest oscillations)
norm2_604 = zeros(3, length(t_list))
@time norm2_lm(norm2_604, 60, 4, ri)
println()
# save norm
outfile = joinpath(out_dir, "norm2_604.h5")
h5open(outfile, "w") do file
    write(file, "norm2", norm2_604)
    attrs(file)["r_causal"] = r1[ri]
end

# all modes included
norm2_all = zeros(3, length(t_list))
@time norm2(norm2_all, ri)
# save norm
outfile = joinpath(out_dir, "norm2_all.h5")
h5open(outfile, "w") do file
    write(file, "norm2", norm2_all)
    attrs(file)["r_causal"] = r1[ri]
end
