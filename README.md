# Example code for "Antiderivative Antialiasing for Stateful Systems"

The example code is provided as a Julia package which can be installed
with
```julia
import Pkg; Pkg.add(Pkg.PackageSpec(url="https://github.com/martinholters/ADAA_Examples.jl.git"))
```
Then, load the package with
```julia
using ADAA_Examples
```
which brings two sub-modules, `DiodeClipper` and `TubeScreamer` into scope.

They both provide three functions `rumsim`, `runsimaaout`, and `runsimaa` to
run the respective model without antialiasing, with antialiasing of the output
only, or with full antialiasing, respectively. They all take two arguments, the
sampling rate (in Hz) and a vector of input samples.

To e.g. recreate the spectra of the tube screamer at 88.2kHz, first install
packages for DSP (for window functions), FFT, and plotting:
```julia
import Pkg; Pkg.add(["DSP", "FFTW", "Plots"])
```
Then do
```julia
using ADAA_Examples, DSP, FFTW, Plots

f0 = 1000
fS = 88100
u = [sin(2π*f0/fS*n) for n in 1:round(Int, fS)]
w = hanning(length(u))
w ./= sum(w)/2
f = range(0, fS/2, length=length(u) ÷ 2 + 1)
kmax = floor(Int, 20000/(fS/2)*length(f))
f = f[1:kmax]

y = TubeScreamer.runsim(fS, u)
Y = rfft(w.*y)[1:kmax]
p1 = plot(f, amp2db.(abs.(Y)); legend=false, title="no aliasing mitigation", ylim=(-130, 10));

y = TubeScreamer.runsimaaout(fS, u)
Y = rfft(w.*y)[1:kmax]
p2 = plot(f, amp2db.(abs.(Y)); legend=false, title="output aliasing mitigation",ylim=(-130, 10));

y = TubeScreamer.runsimaa(fS, u)
Y = rfft(w.*y)[1:kmax]
p3 = plot(f, amp2db.(abs.(Y)); legend=false, title="full aliasing mitigation", ylim=(-130, 10));

plot(p1, p2, p3, layout=(3, 1), size=(600, 800))
```

Note that these functions setup the model and initialize its states on every
invocation, so they are not suitable for successive processing of signal
frames.
