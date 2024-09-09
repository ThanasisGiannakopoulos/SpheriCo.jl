# start time integrators

# Runge-Kutta 4th order (RK4)

#mutate v
function RK4!(t::Float64, dt::Float64,
              #
              v::Array,
              v1::Array,
              v2::Array,
              v3::Array,
              v4::Array,
              sys::System,
              p::Param)

    #v1 =
    F!(t, v1, v, sys, p)

    #v2 =
    F!(t, v2, v +0.5*dt*v1, sys, p)
    
    #v3 =
    F!(t, v3, v +0.5*dt*v2, sys, p)
    
    #v4 =
    F!(t, v4, v +dt*v3, sys, p)
        
    v .= v + dt*(0.166666667*v1 + 0.333333333*v2 + 0.333333333*v3 + 0.166666667*v4)

    nothing
end
function RK4(t::Float64, dt::Float64,
             v::Array,
             sys::System,
             p::Param)

    v1 = similar(v)
    v2 = similar(v)
    v3 = similar(v)
    v4 = similar(v)
    RK4!(t, dt,
         v, v1, v2, v3, v4,
         sys, p)
end

# Adams-Bashford 3rd order (AB3) for the rest of the time evolution
# faster than RK4 but needs smaller dt

#mutate v
function AB3!(
    v::Array,
    F2::Array,
    F1::Array,
    F0::Array,
    dt::Float64)
    
    # build v3
    v .= v .+ dt*(1.916666667*F2 .-1.333333333*F1 .+0.416666667*F0)

    nothing
end

# for a state vector split in v_classic ad v_quantum:
# Runge-Kutta 4th order (RK4)

#mutates v_classic and v_quantum
function RK4!(t::Float64, dt::Float64,
              #
              v_classic::Array,
              v1_classic::Array,
              v2_classic::Array,
              v3_classic::Array,
              v4_classic::Array,
              v_quantum::Array,
              v1_quantum::Array,
              v2_quantum::Array,
              v3_quantum::Array,
              v4_quantum::Array,
              sys::System,
              p::Param)

    #v1 =
    F!(t, v1_classic, v_classic, v1_quantum, v_quantum, sys, p)

    #v2 =
    F!(t, v2_classic, v_classic +0.5*dt*v1_classic,
       v2_quantum, v_quantum +0.5*dt*v1_quantum, sys, p)
    
    #v3 =
    F!(t, v3_classic, v_classic +0.5*dt*v2_classic,
       v3_quantum, v_quantum +0.5*dt*v2_quantum, sys, p)
    
    #v4 =
    F!(t, v4_classic, v_classic +dt*v3_classic,
       v4_quantum, v_quantum +dt*v3_quantum, sys, p)
        
    v_classic .= v_classic + dt*(0.166666667*v1_classic + 0.333333333*v2_classic + 0.333333333*v3_classic + 0.166666667*v4_classic)
    v_quantum .= v_quantum + dt*(0.166666667*v1_quantum + 0.333333333*v2_quantum + 0.333333333*v3_quantum + 0.166666667*v4_quantum)

    nothing
end
function RK4(t::Float64, dt::Float64,
             v_classic::Array,
             v_quantum::Array,
             sys::System,
             p::Param)

    v1_classic = similar(v_classic)
    v2_classic = similar(v_classic)
    v3_classic = similar(v_classic)
    v4_classic = similar(v_classic)
    v1_quantum = similar(v_quantum)
    v2_quantum = similar(v_quantum)
    v3_quantum = similar(v_quantum)
    v4_quantum = similar(v_quantum)
    RK4!(t, dt,
         v_classic, v1_classic, v2_classic, v3_classic, v4_classic,
         v_quantum, v1_quantum, v2_quantum, v3_quantum, v4_quantum,
         sys, p)
end

# Adams-Bashford 3rd order (AB3) for the rest of the time evolution
# faster than RK4 but needs smaller dt

#mutates v_classic and v_quantum
function AB3!(
    v_classic::Array,
    F2_classic::Array,
    F1_classic::Array,
    F0_classic::Array,
    v_quantum::Array,
    F2_quantum::Array,
    F1_quantum::Array,
    F0_quantum::Array,
    dt::Float64)
    
    # build v3
    v_classic .= v_classic .+ dt*(1.916666667*F2_classic .-1.333333333*F1_classic .+0.416666667*F0_classic)
    v_quantum .= v_quantum .+ dt*(1.916666667*F2_quantum .-1.333333333*F1_quantum .+0.416666667*F0_quantum)

    nothing
end

# end time integrators
