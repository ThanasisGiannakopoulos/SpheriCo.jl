# start operators

# FD2 center
function  Dr_FD2!(f_r::Vector, f::Vector, sys::System)
    #Nr = length(f)
    ohr2 = 0.5 / sys.hr

    #r[1], r[2] are ghosts; they will be overwritten
    f_r[1] = 0.0
    f_r[2] = 0.0

    # r=0 is at i=3:
    @inbounds for i in 3:length(f)-1
        f_r[i] = (f[i+1] - f[i-1]) * ohr2
    end

    # end point: 1st order accurate (backward derivarive)
    f_r[end] = (f[end] - f[end-1]) / (sys.hr)

    f_r
end
function Dr_FD2(f::Vector, sys::System)
    f_r = similar(f)
    Dr_FD2!(f_r, f, sys)
end
# multiple dispatch
# for tools
function  Dr_FD2!(f_r::Vector, f::Vector, hr::Float64)
    #Nr = length(f)
    ohr2 = 0.5 / hr

    #r[1], r[2] are ghosts; they will be overwritten
    f_r[1] = 0.0
    f_r[2] = 0.0

    # r=0 is at i=3:
    @inbounds for i in 3:length(f)-1
        f_r[i] = (f[i+1] - f[i-1]) * ohr2
    end

    # end point: 1st order accurate (backward derivarive)
    f_r[end] = (f[end] - f[end-1]) / hr

    f_r
end
function Dr_FD2(f::Vector, hr::Float64)
    f_r = similar(f)
    Dr_FD2!(f_r, f, hr)
end

function  Drr_FD2!(f_rr::Vector, f::Vector, sys::System)
    Nr = length(f)
    odr2 = 1.0 / sys.hr^2

    #r[1], r[2] are ghosts; they will be overwritten
    f_rr[1] = 0.0
    f_rr[2] = 0.0

    @inbounds for i in 3:Nr-1
        f_rr[i] = (f[i+1] -2*f[i] +f[i-1]) *odr2
    end

    # 1st order accurate
    f_rr[end] = (2*f[end] -5*f[end-1] +4*f[end-2] -f[end-3]) *odr2

    f_rr
end
function Drr_FD2(f::Vector, sys::System)
    f_rr = similar(f)
    Drr_FD2!(f_rr, f, sys)
end

function setup_wbar!(wbar::Array{Float64}, L::Int)

    Nr = length(wbar)
    # r[3] = 0.0
    # r[1] = -2h, r[2] = -h, are ghosts
    # L!/2^L is the normalization
    wbar[3] = Float64(1.0)
    for i in 1:L
       wbar[3] = wbar[3]*i
    end
    #wbar0 =
    wbar[3] = wbar[3]* (0.5^L)
    #wbar1 =
    wbar[4] = wbar[3]*(1.0 + L)
    # recursive relation
    @inbounds for i in 5:Nr
        xi = i - 3
        term1 = (2.0* (L+ 1.0)/xi )* ( (1.0 - 1.0/xi )^L)* wbar[i-1]
        term2 = ( (1.0 - 2.0/xi )^(L+1) )* wbar[i-2]
        wbar[i] = term1 + term2
    end

    # the following does not seem needed in fact
    # Ghost points: for L odd, wbar is odd so w is even.
    if mod(L,2) == 0
        wbar[2] = wbar[4]
        wbar[1] = wbar[5]
    else
        wbar[2] = - wbar[4]
        wbar[1] = - wbar[5]
    end

    wbar
end
function setup_wbar(wbar::Array{Float64}, L::Int)
    setup_wbar!(wbar, L)
end

# for real valued
# for odd functions; it is assumed for the l' Hospital on r=0 of the centered grid
function  Dr_SBP2!(f_r::Array{Float64}, f::Array{Float64},
                   sys::System,
                   wbar::Array{Float64}, L::Int)

    Nr = length(f)
    r = sys.r
    ohr = 1.0/sys.hr

    # ghosts for i=1,2 ; they will be overwritten
    f_r[1] = 0.0
    f_r[2] = 0.0

    # at r = 0; that is i=3 here
    # for odd function (vanishes at r=0), l' Hospital for r=0 and centered FD2
    #f_r[3]  = (1.0 + L)*f[4]*ohr # true only for odd
    # use the standard FD stencil
    f_r[3]  = (1.0 + L)*(f[4] - f[2])*ohr*0.5

    @inbounds for i in 4:Nr-1
        xi = i - 3
        f_r[i] = ohr*0.5*(wbar[i+1]/wbar[i])*((1+1/xi)^L)*f[i+1] -
        ohr*0.5*(wbar[i-1]/wbar[i])*((1-1/xi)^L)*f[i-1]
    end

    # end point
    # 1st oder accurate
    xi = Nr-3
    f_r[end]   = ohr*f[end] - ohr*f[end-1]*(wbar[end-1]/wbar[end])*(1-1/xi)^L

    f_r
end

# for complex
function  Dr_SBP2!(f_r::Array{ComplexF64}, f::Array{ComplexF64},
                   sys::System,
                   wbar::Array{Float64}, L::Int)

    Nr = length(f)
    r = sys.r
    ohr = 1.0/sys.hr

    # ghosts for i=1,2 ; they will be overwritten
    f_r[1] = 0.0
    f_r[2] = 0.0

    # at r = 0; that is i=3 here
    # for odd function (vanishes at r=0), l' Hospital for r=0 and centered FD2
    f_r[3]  = (1.0 + L)*f[4]*ohr

    @inbounds for i in 4:Nr-1
        xi = i - 3
        f_r[i] = ohr*0.5*(wbar[i+1]/wbar[i])*((1+1/xi)^L)*f[i+1] -
        ohr*0.5*(wbar[i-1]/wbar[i])*((1-1/xi)^L)*f[i-1]
    end

    # end point
    # 1st oder accurate
    xi = Nr-3
    f_r[end]   = ohr*f[end] - ohr*f[end-1]*(wbar[end-1]/wbar[end])*(1-1/xi)^L

    f_r
end
function Dr_SBP2(f, sys::System, L::Int)

    f_r = similar(f)

    Nr = length(f)
    wbar = zeros(Float64, Nr)
    setup_wbar(wbar, L)

    Dr_SBP2!(f_r, f, sys, wbar, L)
end

# for tools
function  Dr_SBP2!(f_r::Vector, f::Vector,
                   hr::Float64,
                   wbar::Array{Float64}, L::Int)

    Nr = length(f)
    ohr = 1.0/hr

    # ghosts for i=1,2 ; they will be overwritten
    f_r[1] = 0.0
    f_r[2] = 0.0

    # at r = 0; that is i=3 here
    # for odd function (vanishes at r=0), l' Hospital for r=0 and centered FD2
    f_r[3]  = (1.0 + L)*f[4]*ohr

    @inbounds for i in 4:Nr-1
        xi = i - 3
        f_r[i] = ohr*0.5*(wbar[i+1]/wbar[i])*((1+1/xi)^L)*f[i+1] -
        ohr*0.5*(wbar[i-1]/wbar[i])*((1-1/xi)^L)*f[i-1]
    end

    # end point
    # 1st oder accurate
    xi = Nr-3
    f_r[end]   = ohr*f[end] - ohr*f[end-1]*(wbar[end-1]/wbar[end])*(1-1/xi)^L

    f_r
end
function Dr_SBP2(f::Vector, hr::Float64, L::Int)

    f_r = similar(f)

    Nr = length(f)
    wbar = zeros(Float64, Nr)
    setup_wbar(wbar, L)

    Dr_SBP2!(f_r, f, hr, wbar, L)
end

#KO diss for FD2
function diss_FD2!(f_r4::Vector, f::Vector, sys::System)
    Nr = length(f)

    #r[1], r[2] are ghosts; they will be overwritten
    f_r4[1] = 0.0
    f_r4[2] = 0.0

    @inbounds for i in 3:Nr-2
        f_r4[i] = (f[i-2] - 4.0* f[i-1] + 6.0* f[i] - 4.0* f[i+1] + f[i+2]) / (16.0* sys.hr)
    end

    # f_r4[end-1] = (f[end-4] - 4.0* f[end-3] + 6.0* f[end-2] - 4.0* f[end-1] + f[end]) / (16.0* sys.hr) # human error?
    f_r4[end-1] = (f[end-5] - 4.0* f[end-4] + 6.0* f[end-3] - 4.0* f[end-2] + f[end-1]) / (16.0* sys.hr)
    f_r4[end]   = (f[end-4] - 4.0* f[end-3] + 6.0* f[end-2] - 4.0* f[end-1] + f[end]) / (16.0* sys.hr)

    f_r4
end
function KO_FD2(f, sys::System)
    f_r4 = similar(f)
    diss_FD2!(f_r4, f, sys)
end

#KO diss for FD2 acted as a convolution (dot via dt)
function low_pass!(f_r4::Vector, f::Vector)
    Nr = length(f)

    #r[1], r[2] are ghosts; they will be overwritten
    f_r4[1] = 0.0
    f_r4[2] = 0.0

    @inbounds for i in 3:Nr-2
        f_r4[i] = (f[i-2] - 4.0* f[i-1] + 6.0* f[i] - 4.0* f[i+1] + f[i+2]) / 16.0
    end

    # f_r4[end-1] = (f[end-4] - 4.0* f[end-3] + 6.0* f[end-2] - 4.0* f[end-1] + f[end]) / 16.0 # human error
    f_r4[end-1]   = (f[end-5] - 4.0* f[end-4] + 6.0* f[end-3] - 4.0* f[end-2] + f[end-1]) / 16.0
    f_r4[end]   = (f[end-4] - 4.0* f[end-3] + 6.0* f[end-2] - 4.0* f[end-1] + f[end]) / 16.0

    f_r4
end
function low_pass(f)
    f_r4 = similar(f)
    low_pass!(f_r4, f)
end

# projections on a smaller grid using interpolation
# needs old r grid, new r grid and old grid function; returns new grid funtion (old projected on new grid)
function project_on_smaller_grid(r_old::Array, r_new::Array, f_old::Array)

    #interpolation
    f_itp = cubic_spline_interpolation(range(r_old[1], r_old[end], length=length(r_old)), f_old[:])
    # return the projection onto the new grid
    f_itp.(r_new)

end

### not good at rmax; check again

# end operators
