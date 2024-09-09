using HDF5
using LaTeXStrings
using Plots; pythonplot()
using DelimitedFiles
using SpheriCo

using ProgressBars

# give the directory where the data from all the runs are saved
dir = "../examples/classical_runs/"
par = "a0.5_b5.0_c2.0_rmax30_tmax15_cfl0.125_sigma0.0_infalling_rmax_false_rand_false_overMp2_25.132741228718345_damp0"
your_dir = dir*par
# the convention
oMp2 = 8*π #1.0

# set D to your coarse resolution
D = 3
out_dir = "./convergence_plots/"*par*"/resol_D$(D)$(D+1)$(D+2)"
if ispath(out_dir)==false
    mkpath(out_dir)
end

Nr = 128*2^D + 3 # the overal course graining
n = 0

# load the r grid
r_c = h5read(your_dir*"/data_$((Nr-3)*2^n+3)/r.h5","r")
r = r_c
drc = r_c[2] - r_c[1]
r_m = h5read(your_dir*"/data_$((Nr-3)*2^(n+1)+3)/r.h5","r")
r_f = h5read(your_dir*"/data_$((Nr-3)*2^(n+2)+3)/r.h5","r")

# list all available iterations (and corresponding files)
(its_c, all_filenames_c) = list_h5_files(your_dir*"/data_$((Nr-3)*2^n+3)", prefix="data_")
(its_m, all_filenames_m) = list_h5_files(your_dir*"/data_$((Nr-3)*2^(n+1)+3)", prefix="data_")
(its_f, all_filenames_f) = list_h5_files(your_dir*"/data_$((Nr-3)*2^(n+2)+3)", prefix="data_")

# lists used in checking if the projection on the coarse grid is done correctly
#rm_c = zeros(length(rc));
@assert r_c[3:end] ≈ r_m[3:2:end]
@assert r_c[3:end] ≈ r_f[3:4:end]

# v_classic = [Φ, Π, Ψ, A, B, DB, Utld, K, KB, λ, α, Dα, Θ, Zr, f, g, U, V]^T
v_classic_labels = ["Φ", "Π", "Ψ", "A", "B", "DB", "Utld", "K", "KB", "λ", "α", "Dα", "Θ", "Zr", "f", "g", "U", "V"]
# tune fi to choose the function you want from the classic state vector
fi =1

# get index of r_c in coarse resol, for which  r[index] is the first causally connected point
function find_ri_max_causal()
    it_c      = its_c[end]
    it_str_c  = lpad(it_c, 4, "0")
    tc_end    =  h5readattr(your_dir*"/data_$((Nr-3)*2^n + 3)/data_$(it_str_c).h5", "/")["time"]
    ii = 1
    while r_c[ii] != r_c[end] - tc_end
        ii = ii + 1
    end
    ri_max_causal = ii-1
    drc = r_c[4] - r_c[3]
    r_c[ri_max_causal] == r_c[end] - tc_end# - drc
    return ri_max_causal
end
ri_max_causal = find_ri_max_causal()

function f_evol(ti::Int, fi::Int, ri_min, ri_max)
    it_f = its_f[ti] 
    it_str_f  = lpad(it_f, 4, "0")
    # fine resol
    n = 2

    v = h5read(your_dir*"/data_$((Nr-3)*2^n + 3)/data_$(it_str_f).h5","v")
    t = h5readattr(your_dir*"/data_$((Nr-3)*2^n + 3)/data_$(it_str_f).h5", "/")["time"]
    r = h5read(your_dir*"/data_$((Nr-3)*2^n + 3)/r.h5","r")

    t     = h5readattr(your_dir*"/data_$((Nr-3)*2^n + 3)/data_$(it_str_f).h5", "/")["time"]
    f     = v[:,fi]
    fstrg = v_classic_labels[fi]

    A  = v[:,4]
    B  = v[:,5]
    KB = v[:,9]

    AH = find_AH(r, B, A, KB)

    plot(r[ri_min:ri_max], real(f[ri_min:ri_max]), title = "t=$(t)",
         label=fstrg*"(r)", linewidth = 3,
         ylim = (1.1*minimum(f[ri_min:ri_max]), 1.1*maximum(f[ri_min:ri_max])),
         #color = "royalblue",
         frame = true, wsize = (800,600))
    plot!([AH], seriestype="vline", label="AH", linewidth=1.5, color = "black",
         frame = true, wsize = (800,600))

end

function f_conv(i::Int, fi::Int, ri_min, ri_max)
    it_c = its_c[i]
    it_m = its_m[i]
    it_f = its_f[i] 
    it_str_c  = lpad(it_c, 4, "0")
    it_str_m  = lpad(it_m, 4, "0")
    it_str_f  = lpad(it_f, 4, "0")

    tc    =  h5readattr(your_dir*"/data_$((Nr-3)*2^n + 3)/data_$(it_str_c).h5", "/")["time"]
    vc    =  h5read(your_dir*"/data_$((Nr-3)*2^n + 3)/data_$(it_str_c).h5","v")
    vm    =  h5read(your_dir*"/data_$((Nr-3)*2^(n+1) + 3)/data_$(it_str_m).h5","v")
    vf    =  h5read(your_dir*"/data_$((Nr-3)*2^(n+2) + 3)/data_$(it_str_f).h5","v")
    fstrg =  v_classic_labels[fi]
    fc    =  vc[:,fi]
    fm    =  vm[:,fi]
    ff    =  vf[:,fi]
    Af    =  vf[:,4]
    Bf    =  vf[:,5]
    KBf   =  vf[:,9]

    AH = find_AH(r_f, Bf, Af, KBf)

    fcm = fc[3:end] - fm[3:2:end]
    fmf = fm[3:2:end] - ff[3:4:end]

    plot(r[3:ri_max], real(fcm[1:ri_max-2]), title = "t=$(tc)",
         label=fstrg*L"{}_{cm}(r)", linewidth = 2.5,
         ylim = (1.1*minimum(fcm[1:ri_max-2]), 1.1*maximum(fcm[1:ri_max-2])),
         #color = "royalblue",
         frame = true)
    plot!(r[3:ri_max], 4*real(fmf[1:ri_max-2]), title = "t=$(tc)",
          label=L"4\,"*fstrg*L"{}_{mf}(r)", linewidth = 3, 
          ylim =(1.1*minimum(4*real(fmf[1:ri_max-2])), 1.1*maximum(4*real(fmf[1:ri_max-2]))),
          style = :dash,
          #color = "coral",
          frame = true)
    plot!([AH], seriestype="vline", label="AH", linewidth=1, color = "black", frame = true)

end

function plot_cmf(i::Int, ri_min, ri_max)
    it_c = its_c[i]
    it_m = its_m[i]
    it_f = its_f[i] 
    it_str_c  = lpad(it_c, 4, "0")
    it_str_m  = lpad(it_m, 4, "0")
    it_str_f  = lpad(it_f, 4, "0")

    vc = h5read(your_dir*"/data_$((Nr-3)*2^n + 3)/data_$(it_str_c).h5","v")
    vm = h5read(your_dir*"/data_$((Nr-3)*2^(n+1) + 3)/data_$(it_str_m).h5","v")
    vf = h5read(your_dir*"/data_$((Nr-3)*2^(n+2) + 3)/data_$(it_str_f).h5","v")

    tc = h5readattr(your_dir*"/data_$((Nr-3)*2^n + 3)/data_$(it_str_c).h5", "/")["time"]
    tm = h5readattr(your_dir*"/data_$((Nr-3)*2^(n+1) + 3)/data_$(it_str_m).h5", "/")["time"]
    tf = h5readattr(your_dir*"/data_$((Nr-3)*2^(n+2) + 3)/data_$(it_str_f).h5", "/")["time"]

    Φc    = vc[:,1]
    Πc    = vc[:,2]
    Ψc    = vc[:,3]
    Ac    = vc[:,4]
    Bc    = vc[:,5]
    DBc   = vc[:,6]
    Utldc = vc[:,7]
    Kc    = vc[:,8]
    KBc   = vc[:,9]
    λc    = vc[:,10]
    αc    = vc[:,11]
    Dαc   = vc[:,12]

    Φm    = vm[:,1]
    Πm    = vm[:,2]
    Ψm    = vm[:,3]
    Am    = vm[:,4]
    Bm    = vm[:,5]
    DBm   = vm[:,6]
    Utldm = vm[:,7]
    Km    = vm[:,8]
    KBm   = vm[:,9]
    λm    = vm[:,10]
    αm    = vm[:,11]
    Dαm   = vm[:,12]

    Φf    = vf[:,1]
    Πf    = vf[:,2]
    Ψf    = vf[:,3]
    Af    = vf[:,4]
    Bf    = vf[:,5]
    DBf   = vf[:,6]
    Utldf = vf[:,7]
    Kf    = vf[:,8]
    KBf   = vf[:,9]
    λf    = vf[:,10]
    αf    = vf[:,11]
    Dαf   = vf[:,12]

    AH = find_AH(r_f, Bf, Af, KBf)

    Hc, Pc  = constraints_sbp21(Φc, Πc, Ψc, Ac, Bc, DBc, Utldc, Kc, KBc, λc, αc, Dαc, r_c, oMp2);
    Hm, Pm  = constraints_sbp21(Φm, Πm, Ψm, Am, Bm, DBm, Utldm, Km, KBm, λm, αm, Dαm, r_m, oMp2);
    Hf, Pf  = constraints_sbp21(Φf, Πf, Ψf, Af, Bf, DBf, Utldf, Kf, KBf, λf, αf, Dαf, r_f, oMp2);

    plot(r_c[3:ri_max], real(Φc[3:ri_max]), title = "t=$(tc)",  label=L"Φ_c(r)", linewidth = 2,
         ylim = (1.1*minimum(real(Φc[3:ri_max])), 1.1*maximum(real(Φc[3:ri_max]))),
         #color = "royalblue",
         frame = true)
    plot!(r_m[3:2*ri_max], real(Φm[3:2*ri_max]), label=L"Φ_m(r)", linewidth = 2.5,
          ylim = (1.1*minimum(real(Φm[3:2*ri_max])), 1.1*maximum(real(Φm[3:2*ri_max]))),
          style= :dash,
          #color = "coral",
          frame = true)
    plot!(r_f[3:4*ri_max], real(Φf[3:4*ri_max]), label=L"Φ_f(r)", linewidth = 3, 
          ylim =(1.1*minimum(real(Φf[3:4*ri_max])), 1.1*maximum(real(Φf[3:4*ri_max]))),
          style = :dot,
          #color = "darkgreen",
          frame = true)
    p0 = plot!([AH], seriestype="vline", label=L"AH", linewidth=1, color = "black", frame = true)
    
    plot(r_c[3:ri_max], abs.(real(Hc[3:ri_max])), label=L"|H_c(r)|", linewidth=2, 
         ylim = (0, 1.1*maximum(abs.(real(Hc[1:ri_max-2])))),
         #color = "royalblue"
         )
    plot!(r_m[3:2*ri_max], abs.(real(Hm[3:2*ri_max])), label=L"|H_m(r)|", linewidth=2.5, 
          ylim = (0, 1.1*maximum(abs.(real(Hc[3:ri_max])))),
          style = :dash,
          #color = "coral"
          )
    plot!(r_f[3:4*ri_max], abs.(real(Hf[3:4*ri_max])), label=L"|H_f(r)|", linewidth=3.0,
          ylim = (0, 1.1*maximum(abs.(real(Hc[3:ri_max])))),
          style = :dot,
          #color = "darkgreen",
          frame = true)
    p1 = plot!([AH], seriestype="vline", label=L"AH", linewidth=1, color = "black")
    
    plot(r_c[3:ri_max], abs.(real(Pc[3:ri_max])), label=L"|P_c(r)|", linewidth=2,
         ylim = (0, 1.1*maximum(abs.(real(Pc[3:ri_max])))),
         #color = "royalblue"
         )
    plot!(r_m[3:2*ri_max], abs.(real(Pm[3:2*ri_max])), label=L"|P_m(r)|", linewidth=2.5,
          ylim = (0, 1.1*maximum(abs.(real(Pc[3:ri_max])))),
          style = :dash,
          #color = "coral"
          )
    plot!(r_f[3:4*ri_max], abs.(real(Pf[3:4*ri_max])), label=L"|P_f(r)|", linewidth=3.0,
          ylim =(0, 1.1*maximum(abs.(real(Pc[3:ri_max])))),
          style = :dot,
          #color = "darkgreen",
          frame = true)
    p2 = plot!([AH], seriestype="vline", label=L"AH", linewidth=1, color = "black") 
    
    plot(r_c[3:ri_max], real(αc[3:ri_max]), label=L"α_c(r)", linewidth=2, 
         #ylim =(1.1*minimum(real(αc[3:ri_max])), 1.1*maximum(real(αc[3:ri_max]))),
         #color = "royalblue"
         )
    plot!(r_m[3:2*ri_max], real(αm[3:2*ri_max]), label=L"α_m(r)", linewidth=2.5, 
        #ylim =(1.1*minimum(real(αm[3:2*ri_max])), 1.1*maximum(real(αm[3:2*ri_max]))),
          style = :dash,
          #color = "coral"
          )
    plot!(r_f[3:4*ri_max], real(αf[3:4*ri_max]), label=L"α_f(r)", linewidth=3.0,
          ylim =(0.9*minimum(real(αf[3:4*ri_max])),1.1*maximum(real(αf[3:4*ri_max]))),
          style = :dot,
          #color = "darkgreen",
          frame = true, legend =:bottomright)
    p3 = plot!([AH], seriestype="vline", label=L"AH", linewidth=1, color = "black") 
    
    plt = plot(p0, p1, p3, p2, layout = grid(2, 2), wsize = (800,600))#heights=[0.5 ,0.5, 0.5, 0.5]), wsize = (800,600))
end

function cmf_conv(i::Int, ri_min, ri_max)
    it_c = its_c[i]
    it_m = its_m[i]
    it_f = its_f[i] 
    it_str_c  = lpad(it_c, 4, "0")
    it_str_m  = lpad(it_m, 4, "0")
    it_str_f  = lpad(it_f, 4, "0")

    vc = h5read(your_dir*"/data_$((Nr-3)*2^n + 3)/data_$(it_str_c).h5","v")
    vm = h5read(your_dir*"/data_$((Nr-3)*2^(n+1) + 3)/data_$(it_str_m).h5","v")
    vf = h5read(your_dir*"/data_$((Nr-3)*2^(n+2) + 3)/data_$(it_str_f).h5","v")

    tc = h5readattr(your_dir*"/data_$((Nr-3)*2^n + 3)/data_$(it_str_c).h5", "/")["time"]
    tm = h5readattr(your_dir*"/data_$((Nr-3)*2^(n+1) + 3)/data_$(it_str_m).h5", "/")["time"]
    tf = h5readattr(your_dir*"/data_$((Nr-3)*2^(n+2) + 3)/data_$(it_str_f).h5", "/")["time"]

    Φc    = vc[:,1]
    Πc    = vc[:,2]
    Ψc    = vc[:,3]
    Ac    = vc[:,4]
    Bc    = vc[:,5]
    DBc   = vc[:,6]
    Utldc = vc[:,7]
    Kc    = vc[:,8]
    KBc   = vc[:,9]
    λc    = vc[:,10]
    αc    = vc[:,11]
    Dαc   = vc[:,12]

    Φm    = vm[:,1]
    Πm    = vm[:,2]
    Ψm    = vm[:,3]
    Am    = vm[:,4]
    Bm    = vm[:,5]
    DBm   = vm[:,6]
    Utldm = vm[:,7]
    Km    = vm[:,8]
    KBm   = vm[:,9]
    λm    = vm[:,10]
    αm    = vm[:,11]
    Dαm   = vm[:,12]

    Φf    = vf[:,1]
    Πf    = vf[:,2]
    Ψf    = vf[:,3]
    Af    = vf[:,4]
    Bf    = vf[:,5]
    DBf   = vf[:,6]
    Utldf = vf[:,7]
    Kf    = vf[:,8]
    KBf   = vf[:,9]
    λf    = vf[:,10]
    αf    = vf[:,11]
    Dαf   = vf[:,12]

    AH = find_AH(r_f, Bf, Af, KBf)

    Φcm = Φc[3:end] - Φm[3:2:end]
    Φmf = Φm[3:2:end] - Φf[3:4:end]

    αcm = αc[3:end] - αm[3:2:end]
    αmf = αm[3:2:end] - αf[3:4:end]

    Hc, Pc  = constraints_sbp21(Φc, Πc, Ψc, Ac, Bc, DBc, Utldc, Kc, KBc, λc, αc, Dαc, r_c, oMp2);
    Hm, Pm  = constraints_sbp21(Φm, Πm, Ψm, Am, Bm, DBm, Utldm, Km, KBm, λm, αm, Dαm, r_m, oMp2);
    Hf, Pf  = constraints_sbp21(Φf, Πf, Ψf, Af, Bf, DBf, Utldf, Kf, KBf, λf, αf, Dαf, r_f, oMp2);

    Hcm = Hc[3:end] - Hm[3:2:end]
    Hmf = Hm[3:2:end] - Hf[3:4:end]

    Pcm = Pc[3:end] - Pm[3:2:end]
    Pmf = Pm[3:2:end] - Pf[3:4:end]

    plot(r[3:ri_max], real(Φcm[1:ri_max-2]), title = "t=$(tc)",
         label=L"Φ_{cm}(r)", linewidth = 2.5,
         ylim = (1.1*minimum(Φcm[1:ri_max-2]), 1.1*maximum(Φcm[1:ri_max-2])),
         #color = "royalblue",
         frame = true,  xaxis=nothing,
         legendfontsize=10)
    plot!(r[3:ri_max], 4*real(Φmf[1:ri_max-2]), title = "t=$(tc)",  label=L"4\, Φ_{mf}(r)",
          linewidth = 3,
          ylim =(1.1*minimum(4*real(Φmf[1:ri_max-2])), 1.1*maximum(4*real(Φmf[1:ri_max-2]))),
          style = :dash,
          #color = "coral",
          frame = true)
    p0 = plot!([AH], seriestype="vline", label=L"AH", linewidth=1, color = "black", frame = true,
               ytickfont=10)
    
    plot(r[3:ri_max], abs.(real(Hcm[1:ri_max-2])), label=L"|H_{cm}(r)|", linewidth=2.5, 
         ylim = (0, 1.1*maximum(abs.(real(Hcm[1:ri_max-2])))),
         legendfontsize=10
         #color = "royalblue"
         )
    plot!(r[3:ri_max], 4*abs.(real(Hmf[1:ri_max-2])), label=L"4\, |H_{mf}(r)|", linewidth=3.0,
          ylim = (0, 1.1*maximum(4*abs.(real(Hmf[1:ri_max-2])))),
          style = :dash,
          #color = "coral",
          frame = true,  xaxis=nothing)
    p1 = plot!([AH], seriestype="vline", label=L"AH", linewidth=1, color = "black",
               ytickfont=10)
    
    plot(r[3:ri_max], abs.(real(Pcm[1:ri_max-2])), label=L"|P_{cm}(r)|", linewidth=2.5,
         ylim = (0, 1.1*maximum(abs.(real(Pcm[1:ri_max-2])))),
         #color = "royalblue"
         )
    plot!(r[3:ri_max], 4*abs.(real(Pmf[1:ri_max-2])), label=L"4\, |P_{mf}(r)|", linewidth=3.0,
          ylim =(0, 1.1*maximum(4*abs.(real(Pmf[1:ri_max-2])))),
          style = :dash,
          #color = "coral",
          frame = true, legendfontsize=10)
    p2 = plot!([AH], seriestype="vline", label=L"AH", linewidth=1, color = "black",
               xtickfont=10, ytickfont=10,
               xlabel="r") 

    plot(r[3:ri_max], real(αcm[1:ri_max-2]), label=L"α_{cm}(r)", linewidth=2.5, 
         ylim =(1.1*minimum(real(αcm[1:ri_max-2])), 1.1*maximum(real(αcm[1:ri_max-2]))),
         #color = "royalblue"
         )
    plot!(r[3:ri_max], 4*real(αmf[1:ri_max-2]), label=L"4\, α_{mf}(r)", linewidth=3.0,
          ylim =(1.1*minimum(4*real(αmf[1:ri_max-2])),1.1*maximum(4*real(αmf[1:ri_max-2]))),
          style = :dash,
          #color = "coral",
          frame = true, legend =:topright, legendfontsize=10)
    p3 = plot!([AH], seriestype="vline", label=L"AH", linewidth=1, color = "black",
               xtickfont=10, ytickfont=10,
               xlabel="r") 

    plt = plot(p0, p1, p3, p2, layout = grid(2, 2), wsize = (800,600))
end

# make plots and gifs and save

println("remember to choose the function you want by tuning fi")

# gifs

# println("$(v_classic_labels[fi]) evolution gif")
# anim = @animate for i ∈ ProgressBar(1:length(its_c))
#     f_evol(i, fi, 3, 4*ri_max_causal)
# end
# gif(anim, out_dir*"/v_$(fi)_evol.gif", fps = 5)

# println("$(v_classic_labels[fi]) convergence evolution gif")
# anim = @animate for i ∈ ProgressBar(1:length(its_c))
#     f_conv(i, fi, 3, ri_max_causal)
# end
# gif(anim, out_dir*"/v_$(fi)_conv.gif", fps = 5)

# println("Φ, α, H, P evolution gif")
# anim = @animate for i ∈ ProgressBar(1:length(its_c))
#     plot_cmf(i, 3, ri_max_causal)
# end
# gif(anim, out_dir*"/Phi_alpha_H_P_cmf.gif", fps = 5)

# println("Φ, α, H, P convergence evolution gif")
# anim = @animate for i ∈ ProgressBar(1:length(its_c))
#     cmf_conv(i, 3, ri_max_causal)
# end
# gif(anim, out_dir*"/Phi_alpha_H_P_cmf_conv.gif", fps = 5)

# # plots

# println("$(v_classic_labels[fi]) evolution pictures in pdf")
# for i ∈ ProgressBar(1:length(its_c))
#     plt = f_evol(i, fi, 3, 4*ri_max_causal)
#     savefig(plt, out_dir*"/v_$(fi)_evol-$(i).pdf")
# end

# println("Φ, α, H, P evolution pictures in pdf")
# for i in ProgressBar(1:1:length(its_c))
#     plt = plot_cmf(i, 3, ri_max_causal)
#     savefig(plt, out_dir*"/Phi_alpha_H_P_cmf-$(i).pdf")
# end

println("Φ, α, H, P convergence evolution pictures in pdf")
for i in ProgressBar(1:1:length(its_c))
    plt = cmf_conv(i, 3, ri_max_causal)
    savefig(plt, out_dir*"/Phi_alpha_H_P_cmf_conv-$(i).pdf")
end
println("end")
