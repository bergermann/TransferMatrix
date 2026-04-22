
using Plots
using BenchmarkTools

include("src/transfer_matrix.jl")




freqs = range(1e9,30,100);
distances = (1.45:0.05:1.55).*1e-3
tilts = range(-deg2rad(0.01),deg2rad(0.01),3)

coords = Coordinates(1,0.02; diskR=0.15);
modes = Modes(coords,4,3);

tiltsx = [0.]
tiltsy = [0.]

@time gpm = GPM(freqs,distances,tilts,modes,coords);

@time RB = transfer_matrix_3d(Dist,dists,tiltsx,tiltsy,gpm,freqs;);


E0 = modes2field(ax,modes)[:,:,1]
showField(E0,coords)

k0 = 2π*freqs[1]/c0*sqrt(1)
propagate!(E0,k0,coords,10e-3,deg2rad(0.5),0)
showField(E0,coords)
c1 = field2modes(E0,modes);

dists = [
    7.005317,
    7.161926,
    7.436722,
    7.144421,
    7.185010,
    7.209110,
    7.278833,
    7.169816,
    7.250541,
    7.214103,
    7.170475,
    7.245183,
    7.241939,
    7.191030,
    7.208307,
    7.300933,
    7.203299,
    7.265450,
    6.785361,
    7.310886,
]*1e-3

