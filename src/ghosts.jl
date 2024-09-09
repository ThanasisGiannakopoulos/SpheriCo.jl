# start ghosts (for centered grid)
# populate ghost points r[1]=-2h, r[2]=-h

# even | f[3] = f(r=0)
# f[1] = f[5] = f(r=2h)
# f[2] = f[4] = f(r=h)

function even_ghosts(f::Array)
    f[1] = f[5]
    f[2] = f[4]
    f
end

# odd  |  f[3] = f(r=0)
# f[1] = -f[5] = -f(r=2h)
# f[2] = -f[4] = -f(r=h)

function odd_ghosts(f::Array)
    f[1] = -f[5]
    f[2] = -f[4]
    f
end

# end ghosts
