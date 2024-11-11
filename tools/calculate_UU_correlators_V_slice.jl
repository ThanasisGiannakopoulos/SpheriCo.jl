# copied from Jecco.jl/src/AdS5_3_1/AdS5_3_1.jl
# needed to set OPENBLAS_NUM_THREADS=1 automatically
using LinearAlgebra

# always set the number of BLAS threads to 1 upon loading the module. by default
# it uses a bunch of them and we don't want that since they trample over each
# other when solving the nested systems equations. it's much better to thread
# over the loop. see also the discussion here:
# https://github.com/JuliaLang/julia/issues/33409
#
# this saves us the tedious task of always setting OMP_NUM_THREADS=1 before
# launching julia.
function __init__()
    LinearAlgebra.BLAS.set_num_threads(1)
    nothing
end

# needed packages
using HDF5
using DelimitedFiles

using SpheriCo
import Base.Threads.@threads

# you can parallelize
println("Running with number of threads = ", Threads.nthreads())

# give the directory where the data from all the runs are saved
dir = "../examples/quantum_runs/"
par = "a1.25_b0.0_c1.0_rmax30.0_tmax8.0_cfl0.0625_sigma0.02_overMp2_1.0_reg_true_backreact_false_mPV100.0_dk_denom_30_kmax20.0_lmax60.0"
your_dir = dir*par

####################################################################
# set manually according to your simulation
hbar = 1
c = 1
dk = π/30
kmax = 20
lmax = 60
D = 3
Nr = 128*2^D + 3 # the overal course graining
####################################################################

####################################################################
# set manually; maybe first check UV correlators to see how many NU,
# NV is a good choice
NU = 400 # index i, labels rows, where mat[i,j] is a matrix
NV = 400 # index j, labels columns, where mat[i,j] is a matrix
# the V slice; for UU correlator
vi = 50
####################################################################

# load the r grid
r = h5read(your_dir*"/data_$(Nr)/r.h5","r")

# list all available iterations (and corresponding files)
(its, all_filenames) = list_h5_files(your_dir*"/data_$(Nr)", prefix="quantum_");

# index of last timestep saved (in list of its, not iteration number)
ti_max = 69#length(its)
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

# needed for algorithm to transform a t,r function in U,V
#take min max of U(t,r), and V(t,r)
Umin = minimum(U_tr)
Umax = maximum(U_tr)
Vmin = minimum(V_tr)
Vmax = maximum(V_tr)
dU = (Umax - Umin)/(NU-1)
dV = (Vmax - Vmin)/(NV-1)
U_axis = zeros(NU)
for i in 1:NU
    U_axis[i] = Umin + (i-1)*(Umax - Umin)/(NU-1)
end
V_axis = zeros(NV)
for j in 1:NV
    V_axis[j] = Vmin + (j-1)*(Vmax - Vmin)/(NV-1)
end

# initiate UU correlators
UU_crlt = zeros(ComplexF64, (NU, NU))

println()
t1 = time()
ichecked = zeros(lmax) # to measure progress in loops below
# calculate U correlator for vi slice
# buffers for sums to make calculation threads safe
buffers = zeros(ComplexF64, (NU, NU, Threads.nthreads()))
@threads for l in 1:lmax
    # thread id
    id = Threads.threadid()
    ichecked[l] = 1
    print("\rcalculating U-U correlator: $(round((sum(ichecked)/lmax)*100.0, digits=2)) %")
    for k in 1:kmax

        # UV coord transf
        uq_UV = NaN* zeros(ComplexF64, (NU, NV, 4))
        for i in 1:ti_max#length(tlist)
            for j in 1:length(r)
                idown  = Int( round( 1 + (NU - 1)*(U_tr[i,j] - Umin)/(Umax - Umin) ) )
                jleft  = Int( round( 1 + (NV - 1)*(V_tr[i,j] - Vmin)/(Vmax - Vmin) ) )

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

                uq_UV[idown, jleft,:]  = (2*l-1)*r[j]^(l-1)*uq_tr[i,j,l,k,:]
                uq_UV[iup, jleft,:]    = (2*l-1)*r[j]^(l-1)*uq_tr[i,j,l,k,:]
                uq_UV[idown, jright,:] = (2*l-1)*r[j]^(l-1)*uq_tr[i,j,l,k,:]
                uq_UV[iup, jright,:]   = (2*l-1)*r[j]^(l-1)*uq_tr[i,j,l,k,:]
            end # r loop
        end # time loop
        # end UV coord transf

        buffers[:,:,id] +=  dk*(
            uq_UV[:,vi,1].*transpose(conj.(uq_UV[:,vi,1])) .-
            2.0*uq_UV[:,vi,2].*transpose(conj.(uq_UV[:,vi,2])) .+
            2.0*uq_UV[:,vi,3].*transpose(conj.(uq_UV[:,vi,3])) .-
            uq_UV[:,vi,4].*transpose(conj.(uq_UV[:,vi,4]))
        )
     end # end k loops
end # end l loops; has @threads

UU_crlt[:,:] .= sum(buffers[:,:,:], dims=3)

println()
println("Δt = ", time() - t1)
GC.gc()

h5write(your_dir*"/data_$(Nr)/UU_crlt_NU$(NU)_vi$(vi).h5", "crlt", UU_crlt)
