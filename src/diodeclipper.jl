module DiodeClipper

using Interpolations
using NLsolve: nlsolve
using NumericalIntegration: cumul_integrate

export runsim, runsimaaout, runsimaa

# component values
const R = 1e3
const C = 33e-9

# nonlinear function mapping p = x+u to y (depending on used sampling rate)
function nlfunction(fS, p)
    Is = 1e-15
    Vt = 25e-3
    function f!(F, x)
        vD=x[1]
        F[1] = p - (2*R*C*fS + 1)*vD - 2*R*Is*sinh(vD/Vt)
    end
    return nlsolve(f!, [0.0]).zero[1]
end

# precompute antiderivate table for given sampling rate
function makeantiderivtable(fS)
    pmax = 1.5 * (4 + 8R*C*fS) # guesstimate of relevant range
    ps = range(-pmax, pmax, length=1024)
    zs = [nlfunction(fS, p) for p in ps]
    zsint = cumul_integrate(ps, zs)
    return scale(interpolate(zsint, BSpline(Cubic(Line(OnGrid())))), ps)
end

# simulation without any antialiasing
function runsim(fS, u)
    y = similar(u) # preallocate output
    x = 0.0 # initial state
    for n in eachindex(u)
        p = x + u[n]
        y[n] = nlfunction(fS, p)
        x = -x + 4*C*R*fS*y[n]
    end
    return y
end

# simulation with antialiasing of the output only
function runsimaaout(fS, u)
    y = similar(u)
    adtab = makeantiderivtable(fS)
    x = 0.0
    last_p = 0.0
    last_ad = adtab(last_p)
    for n in eachindex(u)
        p = R*x + u[n]
        ad = adtab(p)
        if abs(p-last_p) > 1e-15
            y[n] = (ad - last_ad) / (p - last_p)
        else
            y[n] = nlfunction(fS, 0.5 * (p + last_p))
        end
        last_p = p
        last_ad = ad
        x = -x + 4*C*fS*nlfunction(fS, p)
    end
    return y
end

# simulation with full antialiasing
function runsimaa(fS, u)
    fS = fS/1.5
    y = similar(u)
    adtab = makeantiderivtable(fS)
    x = 0.0
    last_p = 0.0
    last_ad = adtab(last_p)
    prev_x = 0.0
    for n in eachindex(u)
        p = R*x + u[n]
        ad = adtab(p)
        if abs(p-last_p) > 1e-15
            y[n] = (ad - last_ad) / (p - last_p)
        else
            y[n] = nlfunction(fS, 0.5 * (p + last_p))
        end
        last_p = p
        last_ad = ad
        prev_x, x = x, -0.5*(x + prev_x) + 4fS*C*y[n]
    end
    return y
end

end
