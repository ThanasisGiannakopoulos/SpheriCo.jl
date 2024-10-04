# needed packages
using HDF5
using DelimitedFiles

using SpheriCo
import Base.Threads.@threads

# you can parallelize
println("Running with number of threads = ", Threads.nthreads())

# give the directory where the data from all the runs are saved
dir = "../examples/quantum_runs/"
par = "a0.4_b0.0_c4.0_rmax15.0_tmax7.0_cfl0.0625_sigma0.0_overMp2_25.132741228718345_reg_true_backreact_false_mPV200.0_dk_denom_15_kmax10.0_lmax30.0"
your_dir = dir*par

####################################################################
# set manually according to your simulation
hbar = 1
c = 1
dk = π/15
kmax = 10
lmax = 30

D = 3
Nr = 128*2^D + 3 # the overal course graining
####################################################################

####################################################################
# set manually; maybe first check UV correlators to see how many NU,
# NV is a good choice
NU = 400 # index i, labels rows, where mat[i,j] is a matrix
NV = 300 # index j, labels columns, where mat[i,j] is a matrix
####################################################################

# load the r grid
r = h5read(your_dir*"/data_$(Nr)/r.h5","r")

# list all available iterations (and corresponding files)
(its, all_filenames) = list_h5_files(your_dir*"/data_$(Nr)", prefix="quantum_");

# index of last timestep saved (in list of its, not iteration number)
ti_max = length(its)
# initiate matrices to save a function f in t,r
# U,V,f,g, are needed to perform the coordinate transformation to
# double-null
U_tr = zeros(ti_max, length(r));
V_tr = zeros(ti_max, length(r));
f_tr = zeros(ti_max, length(r));
g_tr = zeros(ti_max, length(r));

# the same as above, for all the quantum modes the 4 at the end is for
# the massive ghost fields that do the regularization
uq_tr = zeros(ComplexF64, (ti_max, length(r), lmax+1, kmax, 4) )

t1 = time() # to measure the time passed in calculations
ichecked = zeros(ti_max) # to measure progress in loops below
# load and save function in t,r
for i in 1:ti_max
    ichecked[i] = 1
    print("\rloading uq_tr: $(round((sum(ichecked)/ti_max)*100.0, digits=2)) %")
    it = its[i]
    it_str  = lpad(it, 4, "0")
    
    U   =  h5read(your_dir*"/data_$(Nr)/data_$(it_str).h5","v")[:,17]
    V   =  h5read(your_dir*"/data_$(Nr)/data_$(it_str).h5","v")[:,18]
    f   =  h5read(your_dir*"/data_$(Nr)/data_$(it_str).h5","v")[:,15]
    g   =  h5read(your_dir*"/data_$(Nr)/data_$(it_str).h5","v")[:,16]

    U_tr[i,:] = U
    V_tr[i,:] = V
    f_tr[i,:] = f
    g_tr[i,:] = g

    uq_tr[i,:,:,:,:] = h5read(your_dir*"/data_$(Nr)/quantum_$(it_str).h5","v")[:,:,:,:,1]
    
end
println()
println("Δt = ", time() - t1) # time lapsed in calculation above
GC.gc() # garbage collector; might free up memory space

# needed for algorithm to transform a t,r function in U,V
#take min max of U(t,r), and V(t,r)
Umin = minimum(U_tr)
Umax = maximum(U_tr)
Vmin = minimum(V_tr)
Vmax = maximum(V_tr)

U_axis = zeros(NU)
for i in 1:NU
    U_axis[i] = Umin + (i-1)*(Umax - Umin)/(NU-1)
end

V_axis = zeros(NV)
for j in 1:NV
    V_axis[j] = Vmin + (j-1)*(Vmax - Vmin)/(NV-1)
end

# initiate rescaled quantum modes as nans (not a number)
utld_UV = NaN* zeros(ComplexF64, (NU, NV, lmax+1, kmax, 4) ) # 4 for PV

Umin = minimum(U_tr)
Umax = maximum(U_tr)
Vmin = minimum(V_tr)
Vmax = maximum(V_tr)

dU = (Umax - Umin)/(NU-1)
dV = (Vmax - Vmin)/(NV-1)

println()
t1 = time()
ichecked = zeros(ti_max)
# coordinate transformation for rescaled quantum modes to U,V
@threads for i in 1:ti_max
    ichecked[i] = 1
    print("\rcalculating utld_UV: $(round((sum(ichecked)/ti_max)*100.0, digits=2)) %")
    for j in 1:length(r)

        Utr    = U_tr[i,j]
        Vtr    = V_tr[i,j]
        idown  = Int( round( 1 + (NU - 1)*(Utr - Umin)/(Umax - Umin) ) )
        jleft  = Int( round( 1 + (NV - 1)*(Vtr - Vmin)/(Vmax - Vmin) ) )

        if jleft < NV
            jright = jleft + 1
        else
            jright = jleft
        end
        if idown < NU
            iup =  idown + 1
        else
            iup =  idown
        end
        utld_UV[idown, jleft,1,:,:] = uq_tr[i,j,1,:,:]
        utld_UV[iup, jleft,1,:,:] = uq_tr[i,j,1,:,:]
        utld_UV[idown, jright,1,:,:] = uq_tr[i,j,1,:,:]
        utld_UV[iup, jright,1,:,:] = uq_tr[i,j,1,:,:]
        for l in 2:lmax
            utld_UV[idown, jleft,l,:,:]  = r[j]^(l-1)*uq_tr[i,j,l,:,:]
            utld_UV[iup, jleft,l,:,:]    = r[j]^(l-1)*uq_tr[i,j,l,:,:]
            utld_UV[idown, jright,l,:,:] = r[j]^(l-1)*uq_tr[i,j,l,:,:]
            utld_UV[iup, jright,l,:,:]   = r[j]^(l-1)*uq_tr[i,j,l,:,:]
        end
    end
end
println()
println("Δt = ", time() - t1)
GC.gc()

println()
# correlators for fixed V, across U, takes times
UU_crlt = zeros(ComplexF64, (NV, NU, NU))
t1 = time()
vichecked = zeros(NV)
@threads for vi in 1:NV
    vichecked[vi] = 1
    print("\rcalculating U-U correlators: $(round((sum(vichecked)/NV)*100.0, digits=2)) %")
    for ui in 1:NU
        for uj in 1:NU
            for l in 0:lmax-1
                #for k in 1:kmax
                    UU_crlt[vi, ui, uj] += dk*(2.0*l+1.0)*sum(
                        utld_UV[ui, vi, l+1, :, 1].*conj(utld_UV[uj, vi, l+1, :, 1]) .-
                        2.0*utld_UV[ui, vi, l+1, :, 2].*conj(utld_UV[uj, vi, l+1, :, 2]) .+
                        2.0*utld_UV[ui, vi, l+1, :, 3].*conj(utld_UV[uj, vi, l+1, :, 3]) .-
                        utld_UV[ui, vi, l+1, :, 4].*conj(utld_UV[uj, vi, l+1, :, 4]) )
               # end
            end
        end
    end
end
println()
println("Δt = ", time() - t1)

h5write(your_dir*"/data_$(Nr)/UU_crlt_NU$(NU).h5", "crlt", UU_crlt)

"""
If you want to calculate the correlators for fixed U, across V,
then uncomment bellow
"""

# println()
# # correlators for fixed U
# VV_crlt = zeros(ComplexF64, (NU, NV, NV))
# t1 = time()
# uichecked = zeros(NU)
# @threads for ui in 1:NU
#     uichecked[ui] = 1
#     print("\rcalculating V-V correlators: $(round((sum(uichecked)/NU)*100.0, digits=2)) %")
#     for vi in 1:NV
#         for vj in 1:NV
#             for l in 0:lmax-1
#                 #for k in 1:kmax
#                     VV_crlt[ui, vi, vj] += dk*(2.0*l+1.0)*sum(
#                         utld_UV[ui, vi, l+1, :, 1].*conj(utld_UV[ui, vj, l+1, :, 1]) .-
#                         2.0*utld_UV[ui, vi, l+1, :, 2].*conj(utld_UV[ui, vj, l+1, :, 2]) .+
#                         2.0*utld_UV[ui, vi, l+1, :, 3].*conj(utld_UV[ui, vj, l+1, :, 3]) .-
#                         utld_UV[ui, vi, l+1, :, 4].*conj(utld_UV[ui, vj, l+1, :, 4]) )
#                # end
#             end
#         end
#     end
# end
# println()
# println("Δt = ", time() - t1)

# h5write(your_dir*"/data_$(Nr)/VV_crlt.h5", "crlt", VV_crlt)
# GC.gc()
