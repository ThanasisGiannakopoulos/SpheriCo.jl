#!/bin/bash

# Script that searches for the critical point of the scalar field
# collapse in spherical symmetry Initial scalar field: Gaussian

# Specify the output directory directly
inidir="./bisection"
aupper=1.0 # Upper bound of the interval
alower=0.03125 # Lower bound of the interval
damin=0.01 # Minimum interval size for stopping the bisection
# above if I write 1e-2 it give syntax error: TODO fix

# Parameters
bcscale=7      # Scale for the basic calculator; if damin=1e-5, bcscale >=5, etc
R=5.0          # center of Gaussian (classical scalar ID)
s=1.0          # width of Gaussian
sigma=0.02     # for KO dissipation
cfl_denom=8.0  # cfl = 1/cfl_denom for time integration
damp=0.0       # damp=0 means no damping of H, P constraints
#kk1=0.02      # when damp=1, kk1, kk2 control the amount of H, P constraint damping
#kk2=0.0
#infalling_rmax=true
r_max=12.0
t_max=11.5

# Create the output directory if it does not exist
mkdir -p "${inidir}"

# Create the interval file if it does not exist
intervalfile="${inidir}/AHinterval"
if [ ! -e "${intervalfile}" ]; then
    echo "Creating interval file at ${intervalfile}"
    echo "aupper=${aupper}" > "${intervalfile}"
    echo "alower=${alower}" >> "${intervalfile}"
else
    echo "Reading interval variables from ${intervalfile}"
    source "${intervalfile}"
fi

# Function to run the simulation
run_simulation() {
    local a=$1
    local rundir=$2
    echo "------------------------------------------------------------------------------" >&1
    echo "------------------------------------------------------------------------------" >&2
    echo "Creating directories and parameter files"
    tag="rmax${r_max}_tmax${t_max}_cfl1o${cfl_denom}_sigma${sigma}_damping${damp}_amp${a}_width${s}_rc${R}"

    tagdir="${inidir}/${rundir}/${tag}"
    mkdir -p "${tagdir}"
    parfile="${tagdir}/${tag}.jl"
    
    # Set Nr dynamically
    D=4
    Nr=$(( (128)*2**D + 1 + 2 ))

    ### Write par-file ####
    cat <<EOF > "$parfile"
###############################################################################
# $tag.jl
#
# Evolve scalar field
###############################################################################

using SpheriCo
using SpheriCo.classical

# r=rc centered gaussian
# position of center
rc = ${R}
# amplitude
amp = ${a}
#width
width = ${s}

# change D for number of points
D = ${D}
Nr = ${Nr} # +1 for last point, + 2 for ghosts left of r=0
# noise_amplitude_drop = 0.25
sigma = ${sigma} # for KO diss
cfl_denom = ${cfl_denom}
damp = ${damp}
kk1 = 0.02
kk2 = 0.0
r_max  = ${r_max}
t_max  = ${t_max}
infalling_rmax = true

root_dir = "${tagdir}"

# create the folders where data are saved
out_dir = joinpath(root_dir, "data_${Nr}")
mkpath(out_dir)
cp("${parfile}", joinpath(out_dir, "${tag}.jl"), force=true)

# parameters to be passed in the model
p = Param(
    # discretization parameters
    Nr         = Nr,
    r_max      = r_max,
    t_max      = t_max,
    # directory to save data
    out_dir    = out_dir,
    #CFL
    cfl        = 1.0/cfl_denom,
    # KO diss
    sigma      = sigma,
    # constraint violation
    damping    = damp, # 1.0 (there is constraint damping), or 0.0 (no damping)
    κ1         = kk1,
    κ2         = kk2,
    # convention: 1/Mp^2 = 1.0 or = 8*π
    overMp2    = 8.0*π,
    # cosmological constant; non-zero for backreaction in quantum case
    CC         = 0.0,
    # for Gaussian
    amp        = amp,
    width      = width,
    rc         = rc,
    # infalling_rmax
    infalling_rmax = infalling_rmax,
    # to exit the code if an Apparent horizon is found
    AH = true,
    # how often to save data
    save_data   = true,
    data_every  = 32*2^D,
    # how often to save data
    save_data_r0   = true,
    data_r0_every  = 1*2^D,
    # how often to save data for checkpoint
    save_checkpoint  = true,
    checkpoint_every = 1.0 # this is given in hours   
)

run_classical(p)
EOF

    #### Run the created SpheriCo example ####

    echo "Running for Gaussian initial data with amp=${a} width=${s} rc=${R}"
    Julialog="${tagdir}/${tag}.log"
    echo "julia ${parfile}"
    julia "${parfile}" > "${Julialog}"

    #### Check if file r_AH.txt exists ####
    if [ -e "${tagdir}/data_${Nr}/r_AH.txt" ]; then
        echo "Apparent horizon found"
        aupper="$a"
    else
        echo "No apparent horizon found"
        alower="$a"
    fi

    # Update the interval file
    echo "aupper=${aupper}" > "${intervalfile}"
    echo "alower=${alower}" >> "${intervalfile}"
}

# Initial run for aupper
run_simulation "$aupper"

if [ ! -e "${inidir}/${tag}/data_${Nr}/r_AH.txt" ]; then
    echo "No apparent horizon found for initial aupper. Exiting..."
    exit 1
fi

# Bisection loop for the critical point search
while true; do
    diff=$(echo "${aupper} - ${alower}" | bc -l)
    if (( $(echo "${diff} <= ${damin}" | bc -l) )); then
        break
    fi

    a=$(echo "scale=${bcscale};(${alower} + ${aupper})/2.0" | bc -l)
    run_simulation "$a"
done

echo "Bisection completed. Critical point found between ${alower} and ${aupper}."
