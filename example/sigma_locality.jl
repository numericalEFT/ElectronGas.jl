"""
Check how the sigma function changes with the external momentum within the G0W0 approximation.

Both the z-factor and band mass ratio are calculated.

It is found that that the plasma approximation for the polarization didn't change the result much.

The locality of the z-factor and band mass seems to be protected by the condition that ω_p >> E_F
"""

using ElectronGas
using Plots
using LaTeXStrings
using GreenFunc

rs = 5.0
beta = 1000.0
mass2 = 1e-6

para = Parameter.rydbergUnit(1.0 / beta, rs, 3; Λs=mass2)

# create the benchmark self-energy with the standard polarization
sigma0 = SelfEnergy.G0W0(para; pifunc=Polarization.Polarization0_3dZeroTemp);
kgrid = [k for k in sigma0.spaceGrid.grid if k < para.kF * 2]
z0 = [SelfEnergy.zfactor(para, sigma0; kamp=k)[1] for k in kgrid]
mass0 = [SelfEnergy.bandmassratio(para, sigma0; kamp=k)[1] for k in kgrid]


factorList = [0.1875, 0.75, 3.0, 18.74, 300.0, 30000.0]
# factorList = [1.0, 3.0]

z = zeros(length(factorList), length(kgrid))
mass = zeros(length(factorList), length(kgrid))

plasmafreq(factor) = round(1 / sqrt(factor / 3), sigdigits=3) # fake plasma freq vesus physical one
ef(factor) = round(plasmafreq(factor) * para.ωp / para.EF, sigdigits=3) # fake plasma freq vesus physical EF
# plasmafreq(factor) = round(factor, sigdigits=3)

Threads.@threads for fi in eachindex(factorList)
    factor = factorList[fi]
    println("factor $factor")
    sigma = SelfEnergy.G0W0(para, Euv=1000.0 * para.EF, Nk=12, order=8; pifunc=Polarization.Polarization0_3dZeroTemp_Plasma, factor=factor)
    # sigma = SelfEnergy.G0W0(para; pifunc=Polarization.Polarization0_3dZeroTemp_Plasma, factor=factor)
    for (i, k) in enumerate(kgrid)
        z[fi, i] = SelfEnergy.zfactor(para, sigma; kamp=k)[1]
        # mass[fi, i] = SelfEnergy.bandmassratio(para, sigma; kamp=k)
        mass[fi, i] = SelfEnergy.bandmassratio(para, sigma; kamp=k)[1]
        # mass[fi, i] = sigma.dynamic[]
    end
end

let
    p = plot(xlabel=L"$k/k_F$", ylabel=L"$z_k$", legend=:bottomright)
    plot!(p, kgrid ./ para.kF, z0, label="physical")
    for (fi, factor) in enumerate(factorList[1:end-1])
        plot!(p, kgrid ./ para.kF, z[fi, :], label=L"$\tilde{\omega}_p = %$(plasmafreq(factor))\omega_p = %$(ef(factor))E_F$")
    end
    savefig(p, "sigma_locality_z.pdf")
    display(p)
end

let
    p = plot(xlabel=L"$k/k_F$", ylabel=L"$m^*_k/m_0$", legend=:topright, xlim=[0.1, kgrid[end] / para.kF], ylim=[0.75, 1.401])
    plot!(p, kgrid ./ para.kF, mass0, label="physical")
    for (fi, factor) in enumerate(factorList[1:end-1])
        plot!(p, kgrid ./ para.kF, mass[fi, :], label=L"$\tilde{\omega}_p = %$(plasmafreq(factor))\omega_p = %$(ef(factor))E_F$")
    end
    savefig(p, "sigma_locality_massratio.pdf")
    display(p)
end