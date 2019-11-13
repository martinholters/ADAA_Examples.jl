module TubeScreamer

using ACME: ACME, DiscreteModel, @circuit, capacitor, diode, opamp, resistor, run!,
    voltagesource, voltageprobe
using Interpolations: BSpline, Cubic, Line, OnGrid, interpolate, scale
using NumericalIntegration: cumul_integrate
using StaticArrays: SVector

export runsim, runsimaaout, runsimaa

# circuit description
const circ = @circuit begin
    # power supply
    j3 = voltagesource(9), [+] ⟷ vcc, [-] ⟷ gnd # 9V power supply

    # simplified input stage
    j1 = voltagesource(), [-] ⟷ gnd # input
    c2 = capacitor(18e-9), [1] ⟷ j1[+]
    r5 = resistor(100e3), [1] ⟷ c2[2], [2] ⟷ gnd

    # distortion stage
    ic1a = opamp(), ["in+"] ⟷ c2[2], ["out-"] ⟷ gnd
    d1 = diode(is=4e-9, η=2), [-] ⟷ ic1a["out+"], [+] ⟷ ic1a["in-"] # 1N914
    d2 = diode(is=3e-9, η=2), [-] ⟷ ic1a["in-"] # 1N914
    d3 = diode(is=5e-9, η=2), [+] ⟷ ic1a["out+"], [-] ⟷ d2[+] # 1N914
    r6 = resistor(500e3), [1] ⟷ ic1a["in-"], [2] ⟷ ic1a["out+"]
    c99 = capacitor(51e-12), [1] ⟷ ic1a["in-"], [2] ⟷ ic1a["out+"]
    c4 = capacitor(47e-9), [1] ⟷ ic1a["in-"]
    r7 = resistor(4.7e3), [1] ⟷ c4[2], [2] ⟷ gnd

    # output stage
    c10 = capacitor(1e-6), [1] ⟷ ic1a["out+"]
    r16 = resistor(100e3), [1] ⟷ c10[2], [2] ⟷ gnd
    j2 = voltageprobe(), [+] ⟷ c10[2], [-] ⟷ gnd # output
end

# precompute antiderivate table
function makeantiderivtable(model)
    pmax = 1e-3 # guesstimate of relevant range
    ps = range(-pmax, pmax, length=1024)
    zs = [SVector{ACME.nn(model,1)}(ACME.solve(model.solvers[1], [p])) for p in ps]
    zsint = cumul_integrate(ps, zs)
    return scale(interpolate(zsint, BSpline(Cubic(Line(OnGrid())))), ps)
end

# simulation without any antialiasing
function runsim(fs, u)
    model = DiscreteModel(circ, 1//fs)
    permutedims(run!(model, permutedims(u)))
end

# simulation with antialiasing of the output only
function runsimaaout(fs, u)
    model = DiscreteModel(circ, 1//fs)
    y = similar(u)
    adtab = makeantiderivtable(model)
    last_x = copy(model.x)
    last_p = 0.0
    last_u = 0.0
    last_ad = adtab(last_p)
    for n in eachindex(u)
        p = (model.dqs[1]*model.x + model.eqs[1]*[u[n]])[1]
        ad = adtab(p)
        if abs(p-last_p) > 1e-15
            z = (ad - last_ad) / (p - last_p)
        else
            z = ACME.solve(model.solvers[1], [0.5*(p+last_p)])
        end
        y[n] = (0.5*model.dy*(model.x+last_x) + 0.5*model.ey*[u[n]+last_u] + model.fy*z + model.y0)[1]
        last_p = p
        last_ad = ad
        last_u = u[n]
        copyto!(last_x, model.x)
        copyto!(model.x, model.a * model.x + model.b * [u[n]] + model.c * ACME.solve(model.solvers[1], [p]) + model.x0)
    end
    return y
end

# simulation with full antialiasing
function runsimaa(fs, u)
    fs = fs/1.5
    model = DiscreteModel(circ, 1/fs)
    y = similar(u)
    adtab = makeantiderivtable(model)
    last_x = copy(model.x)
    old_x = copy(model.x)
    last_p = 0.0
    last_u = 0.0
    last_ad = adtab(last_p)
    for n in eachindex(u)
        p = (model.dqs[1]*model.x + model.eqs[1]*[u[n]])[1]
        ad = adtab(p)
        if abs(p-last_p) > 1e-15
            z = (ad - last_ad) / (p - last_p)
        else
            z = ACME.solve(model.solvers[1], [0.5*(p+last_p)])
        end
        y[n] = (0.5*model.dy*(model.x+last_x) + 0.5*model.ey*[u[n]+last_u] + model.fy*z + model.y0)[1]
        copyto!(old_x, model.x)
        copyto!(model.x, 0.5*model.a*(model.x+last_x) + 0.5*model.b*[u[n]+last_u] + model.c*z + model.x0)
        copyto!(last_x, old_x)
        last_p = p
        last_ad = ad
        last_u = u[n]
    end
    return y
end

end
