module Interaction

# using Parameters, GreenFunc
# include(srcdir*"/parameter.jl")
# using .Parameter
# include(srcdir*"/convention.jl")
# using .Convention
# include(srcdir*"/polarization.jl")
# using .Polarization

using ..Parameter, ..Convention, ..Polarization
using ..Parameters, ..CompositeGrids, ..GreenFunc

srcdir = "."
rundir = isempty(ARGS) ? pwd() : (pwd()*"/"*ARGS[1])

export RPA, KO, RPAwrapped, KOwrapped, coulomb

# if !@isdefined Para
#     include(rundir*"/para.jl")
#     using .Para
# end

# println(Parameter.Param)
# @unpack me, kF, rs, e0, beta , Λs, ϵ0= Parameter.Param

function inf_sum(q,n)
    # Calculate a series sum for Takada anzats
    # See Takada(doi:10.1103/PhysRevB.47.5202)(Eq.2.16).
    a=q*q
    sum=1.0
    i=0
    j=1.0
    k=2.0
    for i in 1:n
        sum = sum+a/j/k
        a=a*q*q
        j=j*(i+1.0)
        k=k*(i+2.0)
    end
    return 1.0/sum/sum
end

"""
    function coulomb(q,param)

Bare interaction in momentum space. Coulomb interaction if Λs=0, Yukawa otherwise.

#Arguments:
 - q: momentum
 - param: other system parameters
"""
function coulomb(q, param)
    @unpack me, kF, rs, e0s, e0a, beta, Λs, Λa, ϵ0 = param
    return e0s^2/ϵ0/(q^2+Λs), e0a^2/ϵ0/(q^2+Λa)
end

function bubbledyson(V, V0, G, Π, n)
    # V:bare interaction
    # V0:symmetric bare interaction
    # G:G^{+-} is local field factor,0 for RPA
    # Π:Polarization. 2*Polarization0 for spin 1/2
    # n:matfreq. special case for n=0
    # comparing to previous convention, an additional V is multiplied
    K = 0
    if n == 0
        K = - V0*Π*(V/V0-G)^2/( 1.0/V0  + Π*(V/V0-G))
    else
        K = - (Π)*(V-V0*G)^2/( 1.0  + (Π)*(V-V0*G))
    end
    return K
end

function bubblecorrection(q, n, param;
            pifunc=Polarization0_ZeroTemp, gfactorfunc=localFieldFactorTakada, V_Bare=coulomb, type=:s)
    Gs::Float64, Ga::Float64 = gfactorfunc(q, n, param)
    Ks::Float64, Ka::Float64 = 0.0, 0.0
    Vs::Float64, Va::Float64 = V_Bare(q, param)
    V0 = (type==:s) ? Vs : Va
    if abs(q) > EPS 
        Π::Float64 = 2*pifunc(q, n, param)
        Ks = bubbledyson(Vs,V0,Gs,Π,n)
        Ka = bubbledyson(Va,V0,Ga,Π,n)
    else
        Ks, Ka = 0.0, 0.0
    end

    return Ks, Ka
end

"""
    function RPA(q, n, param)

Dynamic part of RPA interaction, with polarization approximated by zero temperature Π0.

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function RPA(q, n, param; pifunc=Polarization0_ZeroTemp, V_Bare=coulomb)
    return bubblecorrection(q,n,param;pifunc=pifunc,gfactorfunc=localFieldFactor0,V_Bare=V_Bare)
end

function RPAwrapped(Euv, rtol, sgrid::SGT, param;
                   pifunc=Polarization0_ZeroTemp,gfactorfunc=localFieldFactorTakada, V_Bare=coulomb) where{SGT}
    @unpack me, kF, rs, e0, beta , Λs, ϵ0 = param
    gs = GreenFunc.Green2DLR{Float64}(:rpa,GreenFunc.IMFREQ,beta,false,Euv,sgrid,1; timeSymmetry=:ph,rtol=rtol)
    ga = GreenFunc.Green2DLR{Float64}(:rpa,GreenFunc.IMFREQ,beta,false,Euv,sgrid,1; timeSymmetry=:ph,rtol=rtol)
    green_dyn_s = zeros(Float64, (gs.color, gs.color, gs.spaceGrid.size, gs.timeGrid.size))
    green_ins_s = zeros(Float64, (gs.color, gs.color, gs.spaceGrid.size))
    green_dyn_a = zeros(Float64, (ga.color, ga.color, ga.spaceGrid.size, ga.timeGrid.size))
    green_ins_a = zeros(Float64, (ga.color, ga.color, ga.spaceGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(gs.dlrGrid.n)
            green_dyn_s[1,1,ki,ni],green_dyn_a[1,1,ki,ni] = RPA(k, n, param; pifunc=pifunc, V_Bare=V_Bare)
        end
        green_ins_s[1,1,ki],green_ins_a[1,1,ki] = V_Bare(k, param)
    end
    gs.dynamic=green_dyn_s
    gs.instant=green_ins_s
    ga.dynamic=green_dyn_a
    ga.instant=green_ins_a
    return gs, ga
end

"""
    function localFieldFactorTakada(q, n, param)

G factor with Takada's anzats. See Takada(doi:10.1103/PhysRevB.47.5202)(Eq.2.13-2.16).

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function localFieldFactorTakada(q, n, param)
    @unpack me, kF, rs, e0, beta , Λs, ϵ0= param
    r_s_dl=sqrt(4*0.521*rs/ π );
    C1=1-r_s_dl*r_s_dl/4.0*(1+0.07671*r_s_dl*r_s_dl*((1+12.05*r_s_dl)*(1+12.05*r_s_dl)+4.0*4.254/3.0*r_s_dl*r_s_dl*(1+7.0/8.0*12.05*r_s_dl)+1.5*1.363*r_s_dl*r_s_dl*r_s_dl*(1+8.0/9.0*12.05*r_s_dl))/(1+12.05*r_s_dl+4.254*r_s_dl*r_s_dl+1.363*r_s_dl*r_s_dl*r_s_dl)/(1+12.05*r_s_dl+4.254*r_s_dl*r_s_dl+1.363*r_s_dl*r_s_dl*r_s_dl));
    C2=1-r_s_dl*r_s_dl/4.0*(1+r_s_dl*r_s_dl/8.0*(log(r_s_dl*r_s_dl/(r_s_dl*r_s_dl+0.990))-(1.122+1.222*r_s_dl*r_s_dl)/(1+0.533*r_s_dl*r_s_dl+0.184*r_s_dl*r_s_dl*r_s_dl*r_s_dl)));
    D=inf_sum(r_s_dl,100);
    A1=(2.0-C1-C2)/4.0/e0^2*π ;
    A2=(C2-C1)/4.0/e0^2*π ;
    B1=6*A1/(D+1.0);
    B2=2*A2/(1.0-D);
    G_s=A1*q^2/(1.0+B1*q^2)+A2*q^2/(1.0+B2*q^2);
    G_a=A1*q^2/(1.0+B1*q^2)-A2*q^2/(1.0+B2*q^2);
    return G_s, G_a
end

@inline function localFieldFactor0(q, n, param)
    return 0.0, 0.0
end

"""
    function KO(q, n, param)

Dynamic part of KO interaction, with polarization approximated by zero temperature Π0.
Returns the spin symmetric part and asymmetric part separately.

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function KO(q, n, param; pifunc=Polarization0_ZeroTemp, gfactorfunc=localFieldFactorTakada, V_Bare=coulomb)
    return bubblecorrection(q,n,param;pifunc=pifunc,gfactorfunc=gfactorfunc,V_Bare=coulomb)
end

function KOwrapped(Euv, rtol, sgrid::SGT, param;
                   pifunc=Polarization0_ZeroTemp,gfactorfunc=localFieldFactorTakada, V_Bare=coulomb) where{SGT}
    @unpack me, kF, rs, e0, beta , Λs, ϵ0 = param
    gs = GreenFunc.Green2DLR{Float64}(:ko,GreenFunc.IMFREQ,beta,false,Euv,sgrid,1; timeSymmetry=:ph,rtol=rtol)
    ga = GreenFunc.Green2DLR{Float64}(:ko,GreenFunc.IMFREQ,beta,false,Euv,sgrid,1; timeSymmetry=:ph,rtol=rtol)
    green_dyn_s = zeros(Float64, (gs.color, gs.color, gs.spaceGrid.size, gs.timeGrid.size))
    green_ins_s = zeros(Float64, (gs.color, gs.color, gs.spaceGrid.size))
    green_dyn_a = zeros(Float64, (ga.color, ga.color, ga.spaceGrid.size, ga.timeGrid.size))
    green_ins_a = zeros(Float64, (ga.color, ga.color, ga.spaceGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(gs.dlrGrid.n)
            green_dyn_s[1,1,ki,ni],green_dyn_a[1,1,ki,ni] = KO(k, n, param; pifunc=pifunc, gfactorfunc=gfactorfunc,V_Bare=V_Bare)
        end
        green_ins_s[1,1,ki],green_ins_a[1,1,ki] = V_Bare(k, param)
    end
    gs.dynamic=green_dyn_s
    gs.instant=green_ins_s
    ga.dynamic=green_dyn_a
    ga.instant=green_ins_a
    return gs, ga
end

end
