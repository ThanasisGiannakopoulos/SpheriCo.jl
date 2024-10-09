
using SpheriCo
using SpheriCo.quantum

# amplitude
"""
if overMp2 = 8.0*π  and amp = A
then the same is 
overMp2 = 1 and amp = A*sqrt(8.0*π)
"""
a = 0.4
# position of center
b = 0.0
#width
c = 4.0

# change D for number of points
D = 3
Nr = (128)*2^D + 1 + 2 # +1 for last point, + 2 for ghosts left of r=0
# noise_amplitude_drop = 0.25
sigma = 0.0 # for KO diss
cfl_denom = 16.0
rmax  = 15.0
r_cut = 12.0
tmax  = 7.0
# infalling r_max
infalling_rmax = true
# regularization
reg = false
# backreaction; needs non-zero Cosmological Constant (CC)
backreact = false
# Pauli-Villars mass
mPV = 100.0
# denominator of dk
dk_denom = 15

# convention for G Newton
# 1/Mp^2 = 1 or 8π if there is no backreaction
#overMp2 = 8.0*π
"""
 if there is backreaction
 1/Mp^2|simulation = 1/(Mp^2|physical - mPV^2*ln(2^4/3^3)/(12*4*π^2))
 choose appropriately
"""
PhysMp2 = 1.0 #1/(8*π)
overMp2 = 1/(PhysMp2 - mPV^2*log(2^4/3^3)/(12*4*π^2))

# Cosmological Constant = -ln(3^9/2^16)*mPV^4/(8*(2pi)^2)
# NEEDED for backreaction
CC = 0.0 #-log(3^9/2^16)*(mPV^4)/(8*(2*π)^2)
# above in CC remember to devide by the Mpl^2 (use the same convention as for the rest)
# quantum modes
kmax = 10.0
lmax = 30.0
if backreact==false
    root_dir = "./quantum_runs/a$(a)_b$(b)_c$(c)_rmax$(rmax)_tmax$(tmax)_cfl$(1.0/cfl_denom)_sigma$(sigma)_overMp2_$(overMp2)_reg_$(reg)_backreact_$(backreact)_mPV$(mPV)_dk_denom_$(dk_denom)_kmax$(kmax)_lmax$(lmax)/"
end
if backreact==true
    if CC==0.0
        println("tune CC correctly")
        exit()
    end
    root_dir = "./quantum_runs/a$(a)_b$(b)_c$(c)_rmax$(rmax)_tmax$(tmax)_cfl$(1.0/cfl_denom)_sigma$(sigma)_overMp2_$(overMp2)_reg_$(reg)_backreact_$(backreact)_rcut$(r_cut)_mPV$(mPV)_dk_denom_$(dk_denom)_kmax$(kmax)_lmax$(lmax)/"
end

# create the folders where data are saved
out_dir = joinpath(root_dir, "data_$(Nr)/")
mkpath(out_dir)
cp("./evol_quantum.jl", out_dir*"/evol_quantum.jl", force=true)

# radial grid
g = Grid(
    # discretization parameters
    Nr         = Nr,
    r_max      = rmax
)

# parameters to be passed in the model
p = Param(
    # time of simulation
    t_max      = tmax,
    # directory to save data
    out_dir    = out_dir,
    #CFL
    cfl        = 1.0/cfl_denom,
    # KO diss
    sigma      = sigma,
    # constraint violation
    damping    = 0.0, # 1.0 (there is constraint damping), or 0.0 (no damping)
    κ1         = 0.0,
    κ2         = 0.0,
    # convention: 1/Mp^2 = 1.0 or = 8*π
    # the 8pi convention does not work well with backreaction; related to CC reg?
    overMp2    = overMp2,
    # cosmological constant; non-zero for backreaction in quantum case
    CC         = CC,
    # for Gaussian
    amp        = a,
    rc         = b,
    width      = c,
    # infalling_rmax
    infalling_rmax = infalling_rmax,
    # exit the code if an Apparent horizon is found
    AH = false,
    # how often to save data
    save_data   = true,
    data_every  = 8*2^D,
    # how often to save data
    save_data_r0   = false,
    data_r0_every  = 1*2^D,
    # how often to save data for checkpoint
    save_checkpoint  = true,
    checkpoint_every = 1.0, # this is given in hours
    ##########################################################
    # quantum
    hbar = 1.0, # default = 1.0
    # for filter based on tanh, to cut backreaction in causally disconnecter r-part
    # the combo of steepness and r_cut has to give a filter with filter[end]=0
    # to preserve the flat asymptotics of the problem
    # otherwise the BC are wrong for the PDE setup and noise comes from rmax
    # needed for stability
    steepness = 1.5,
    # also for filter; the r where the transition from 1->0 happens
    r_cut = r_cut,
    # number of quantum modes
    kmax = kmax,
    lmax = lmax,
    # for quantum modes
    dk = π/dk_denom, # for the k in the Bessel functions for the quantum ID
    # chose quantum version (regularized or non-regularized)
    PV_reg = reg, # false=non-reg., true=reg.
    # masses for PV regularization are in principle [m0, m1, m2, m3, m4, m5]
    # m0 = 0.0
    # m1 = m3, m2 = m4 = sqrt(3)*m1, m5 = 2.0*m1
    # m1 = mPV = 1.0 as default. It can be changed in the example
    # to avoid repetition we have mlist = [m0, m1, m2, m5]
    mlist = [0.0, 1.0*mPV, sqrt(3.0)*mPV, 2.0*mPV],
    # backreaction
    backreaction = backreact,
    # how often to save bilinears and correlators
    save_quantum      = true,
    quantum_every     = 8*2^D,
    save_quantum_r0   = false,
    quantum_r0_every  = 32*2^D,
    save_bilinears    = true,
    bilinears_every   = 8*2^D,
    save_correlators  = true,
    correlators_every = 8*2^D

)

@time run_quantum(g, p)
