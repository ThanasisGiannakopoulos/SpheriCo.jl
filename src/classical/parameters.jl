@with_kw struct Param
    # max evolution time
    t_max      :: Float64 = 10.0
    # directory to save data
    out_dir    :: String
    # CFL
    cfl        :: Float64 = 0.125
    # KO diss
    sigma      :: Float64 = 0.02
    # constraint violation damping
    damping    :: Float64 = 0.0 # 1.0 (there is constraint damping), or 0.0 (no damping)
    κ1         :: Float64 = 0.02
    κ2         :: Float64 = 0.0
    # convention: 1/Mp^2 = 1.0 or = 8*π
    overMp2    :: Float64 = 8*π
    # cosmological constant, always 0 here (asymptotically flat)
    CC         :: Float64 = 0.0
    # classical Gaussian profile (amplitude, width, location of center)
    amp        :: Float64
    width      :: Float64
    rc         :: Float64
    # random data amplitude (for robust stability test)
    rand   :: Bool = false
    A_rand :: Float64 = 1.0
    # infalling_rmax
    infalling_rmax :: Bool = false
    # to exit the code if an Apparent horizon is found
    AH :: Bool = false
    # how often to save data
    save_data   :: Bool = true
    data_every  :: Int
    # how often to save data at r0; useful for critical collapse
    save_data_r0   :: Bool = false
    data_r0_every  :: Int
    # how often to save checkpoint
    save_checkpoint  :: Bool = true
    checkpoint_every :: Float64 = 1.0 # this is given in hours

end
