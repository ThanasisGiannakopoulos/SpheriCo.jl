using HDF5
using LaTeXStrings
using Plots#; pythonplot()
using DelimitedFiles
using SpheriCo
using ProgressBars

# give the directory where the data from all the runs are saved
dir = "../examples/classical_runs/"
par = "a0.5_b5.0_c2.0_rmax30_tmax15_cfl0.125_sigma0.02_infalling_rmax_false_rand_false_overMp2_25.132741228718345_damp0"
your_dir = dir*par

out_dir = "./convergence_plots/"*par
if ispath(out_dir)==false
    mkpath(out_dir)
end

# which are the coarse resolutions?
# you need 3 resols to have a self-convergence rate
coarse = [3] #[2, 3, 4, 5]

# initiate l2 matrices for calculations later
D1 = coarse[1]
Nr1 = 128*2^D1 + 3 # the overal course graining
# list all available iterations (and corresponding files)
(its_c1, all_filenames_c1) = list_h5_files(your_dir*"/data_$(Nr1)", prefix="data_")

l2_cm  = zeros( (length(coarse), length(its_c1) ))
l2_mf  = zeros( (length(coarse), length(its_c1) ))
t_list = zeros( (length(coarse), length(its_c1) ))
t_AH_list = zeros(length(coarse))

function v_cmf(D::Int, i::Int,
               its_c::Array, its_m::Array, its_f::Array)
    it_c = its_c[i]
    it_m = its_m[i]
    it_f = its_f[i] 
    it_str_c  = lpad(it_c, 4, "0")
    it_str_m  = lpad(it_m, 4, "0")
    it_str_f  = lpad(it_f, 4, "0")

    Nrc = 128*2^D + 3
    Nrm = 128*2^(D+1) + 3
    Nrf = 128*2^(D+2) + 3

    vc = h5read(your_dir*"/data_$(Nrc)/data_$(it_str_c).h5","v")
    vm = h5read(your_dir*"/data_$(Nrm)/data_$(it_str_m).h5","v")
    vf = h5read(your_dir*"/data_$(Nrf)/data_$(it_str_f).h5","v")

    tc = h5readattr(your_dir*"/data_$(Nrc)/data_$(it_str_c).h5", "/")["time"]
    tm = h5readattr(your_dir*"/data_$(Nrm)/data_$(it_str_m).h5", "/")["time"]
    tf = h5readattr(your_dir*"/data_$(Nrf)/data_$(it_str_f).h5", "/")["time"]

    # make sure we're comparing the same timestep
    @assert tc ≈ tm
    @assert tc ≈ tf
    
    vcm = vc[3:1:end,:] - vm[3:2:end,:]
    vmf = vm[3:2:end,:] - vf[3:4:end,:]
    return vcm, vmf
end

# get index of r_c in coarse resol, for which  r[index] is the first causally connected point
function find_ri_max_causal(Nr::Int, its_c::Array)
    it_c      = its_c[end]
    it_str_c  = lpad(it_c, 4, "0")
    tc_end    =  h5readattr(your_dir*"/data_$(Nr)/data_$(it_str_c).h5", "/")["time"]
    ii = 1
    r_c = h5read(your_dir*"/data_$(Nr)/r.h5","r")
    while r_c[ii] != r_c[end] - tc_end
        ii = ii + 1
    end
    ri_max_causal = ii-1
    drc = r_c[4] - r_c[3]
    r_c[ri_max_causal] == r_c[end] - tc_end - drc
    return ri_max_causal
end
# get time that AH first appears
function find_t_AH(Nr::Int, its_c::Array)
    j = 1
    it_c = its_c[j]
    it_str_c  = lpad(it_c, 4, "0")

    vc    = h5read(your_dir*"/data_$(Nr)/data_$(it_str_c).h5","v")
    tc    = h5readattr(your_dir*"/data_$(Nr)/data_$(it_str_c).h5", "/")["time"]
    r_c = h5read(your_dir*"/data_$(Nr)/r.h5","r")
    Ac    = vc[:,4]
    Bc    = vc[:,5]
    KBc   = vc[:,9]
    AH    = find_AH(r_c, Bc, Ac, KBc)

    while AH < 0
        it_c = its_c[j]
        it_str_c  = lpad(it_c, 4, "0")
        vc = h5read(your_dir*"/data_$(Nr)/data_$(it_str_c).h5","v")
        tc = h5readattr(your_dir*"/data_$(Nr)/data_$(it_str_c).h5", "/")["time"]
        Ac    = vc[:,4]
        Bc    = vc[:,5]
        KBc   = vc[:,9]
        AH = find_AH(r_c, Bc, Ac, KBc)
        j = j + 1  
    end
    t_AH = tc
    return t_AH
end

println("Calculate norms")
for Di in ProgressBar(1:length(coarse))
    D = coarse[Di]
    Nrc = 128*2^D + 3 # coarse
    Nrm = 128*2^(D+1) + 3 # medium
    Nrf = 128*2^(D+2) + 3 # fine
    
    # load the r grid
    r_c = h5read(your_dir*"/data_$(Nrc)/r.h5","r")
    drc = r_c[2] - r_c[1]
    r_m = h5read(your_dir*"/data_$(Nrm)/r.h5","r")
    r_f = h5read(your_dir*"/data_$(Nrf)/r.h5","r")

    # list all available iterations (and corresponding files)
    (its_c, all_filenames_c) = list_h5_files(your_dir*"/data_$(Nrc)", prefix="data_")
    (its_m, all_filenames_m) = list_h5_files(your_dir*"/data_$(Nrm)", prefix="data_")
    (its_f, all_filenames_f) = list_h5_files(your_dir*"/data_$(Nrf)", prefix="data_")

    # lists used in checking if the projection on the coarse grid is done correctly
    #rm_c = zeros(length(rc));
    @assert r_c[3:end] ≈ r_m[3:2:end]
    @assert r_c[3:end] ≈ r_f[3:4:end]

    # get index of r_c in coarse resol, for which  r[index] is the first causally connected point
    ri_max_causal = find_ri_max_causal(Nrc, its_c)
    t_AH_list[Di] = find_t_AH(Nrc, its_c)
    
    # it starts from r_c[3], so shift the ri_max_causal by 3
    for ti ∈ 1:length(its_c)
        v_cm, v_mf = v_cmf(D, ti, its_c, its_m, its_f)
        it_c       = its_c[ti]
        it_str_c   = lpad(it_c, 4, "0")
        l2_cm[Di, ti]   = sqrt(sum(v_cm[1:ri_max_causal-3].^2)*drc)
        l2_mf[Di, ti]   = sqrt(sum(v_mf[1:ri_max_causal-3].^2)*drc)
        t_list[Di, ti]  =  h5readattr(your_dir*"/data_$(Nrc)/data_$(it_str_c).h5", "/")["time"]
    end

end

function save_plots(t_list::Array, l2_cm::Array,
                    l2_mf::Array, t_AH_list::Array)
    println("making plot")
    Dc = coarse[1]
    Dm = Dc+1
    Df = Dm+1
    plot(t_list[1, 2:end], log2.(l2_cm[1, 2:end]./l2_mf[1, 2:end]),
         label="D=$Dc,$Dm,$Df",
         linewidth = 2,
         ylim =(0.0,4.0),  frame = true, legend =:bottomright,
         xlabel = "time",
         ylabel=L"\log_2 \frac{ ||\mathbf{v}_c-\mathbf{v}_m||}{||\mathbf{v}_m-\mathbf{v}_f||}")
    plt = plot!([t_AH_list[1]], seriestype="vline", label="t_AH", linewidth=1, color = "black") 
    if length(coarse)>1
        for i in 2:length(coarse)
            Dc = coarse[i]
            Dm = Dc+1
            Df = Dm+1
            plt = plot!(t_list[i, 2:end], log2.(l2_cm[i, 2:end]./l2_mf[i, 2:end]),
                        label="D=$Dc,$Dm,$Df",
                        linewidth = 2,
                        ylim =(0.0,4.0), frame = true, legend =:bottomright,
                        xlabel = "time",
                        ylabel=L"\log_2 \frac{ ||\mathbf{v}_c-\mathbf{v}_m||}{ ||\mathbf{v}_m-\mathbf{v}_f|| }")
        end
    end
    println("saving plot")
    savefig(plt, out_dir*"/L2_norms.pdf")
end

save_plots(t_list, l2_cm, l2_mf, t_AH_list)
