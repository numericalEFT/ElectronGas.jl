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

export RPA, KO, RPAwrapped, KOwrapped, V_Bare

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
    function V_Bare(q,param)

Bare interaction in momentum space. Coulomb interaction if Λs=0, Yukawa otherwise.

#Arguments:
 - q: momentum
 - param: other system parameters
"""
function V_Bare(q, param)
    @unpack me, kF, rs, e0, beta, Λs, ϵ0 = param
    return e0^2/ϵ0/(q^2+Λs)
end


"""
    function RPA(q, n, param)

Dynamic part of RPA interaction, with polarization approximated by zero temperature Π0.

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function RPA(q, n, param; pifunc=Polarization0_ZeroTemp, V_Bare=V_Bare)
    ks = 0.0
    Vs = V_Bare(q, param)
    if abs(q) > EPS 
        Π = pifunc(q, n, param)
        if n == 0
            ks = - Π/( 1.0/Vs + Π )
        else
            ks = - (Π*Vs)/( 1.0  + Π*Vs )
        end
    else
        ks = 0.0
    end
    return ks
end

function RPAwrapped(Euv, rtol, sgrid::SGT, param; pifunc=Polarization0_ZeroTemp, V_Bare=V_Bare) where{SGT}
    @unpack beta = param

    gs = GreenFunc.Green2DLR{Float64}(:rpa,GreenFunc.IMFREQ,beta,false,Euv,sgrid,1; timeSymmetry=:ph,rtol=rtol)
    green_dyn_s = zeros(Float64, (gs.color, gs.color, gs.spaceGrid.size, gs.timeGrid.size))
    green_ins_s = zeros(Float64, (gs.color, gs.color, gs.spaceGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(gs.dlrGrid.n)
            green_dyn_s[1,1,ki,ni] = RPA(k, n, param; pifunc=pifunc,V_Bare=V_Bare)
        end
        green_ins_s[1,1,ki] = V_Bare(k, param)
    end
    gs.dynamic=green_dyn_s
    gs.instant=green_ins_s
    return gs
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

"""
    function KO(q, n, param)

Dynamic part of KO interaction, with polarization approximated by zero temperature Π0.
Returns the spin symmetric part and asymmetric part separately.

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function KO(q, n, param; pifunc=Polarization0_ZeroTemp, gfactorfunc=localFieldFactorTakada, V_Bare=V_Bare)

    G_s::Float64, G_a::Float64 = gfactorfunc(q, n, param)
    Ks::Float64, Ka::Float64 = 0.0, 0.0
    Vs::Float64 = V_Bare(q, param)
    if abs(q) > EPS 
        Π::Float64 = pifunc(q, n, param)
        if n == 0
            Ks = - Π*(1-G_s)^2/( 1.0/Vs  + Π*(1-G_s))
            Ka = -Π*(-G_a)^2/( 1.0/Vs + Π*(-G_a))
        else
            Ks = - (Π*Vs)*(1-G_s)^2/( 1.0  + (Π*Vs)*(1-G_s))
            Ka = -(Π*Vs)*(-G_a)^2/(1.0 + (Π*Vs)*(-G_a))
        end
    else
        Ks, Ka = 0.0, 0.0
    end

    return Ks, Ka
end

function KOwrapped(Euv, rtol, sgrid::SGT, param;
                   pifunc=Polarization0_ZeroTemp,gfactorfunc=localFieldFactorTakada, V_Bare=V_Bare) where{SGT}
    @unpack me, kF, rs, e0, beta , Λs, ϵ0 = param
    gs = GreenFunc.Green2DLR{Float64}(:rpa,GreenFunc.IMFREQ,beta,false,Euv,sgrid,1; timeSymmetry=:ph,rtol=rtol)
    ga = GreenFunc.Green2DLR{Float64}(:rpa,GreenFunc.IMFREQ,beta,false,Euv,sgrid,1; timeSymmetry=:ph,rtol=rtol)
    green_dyn_s = zeros(Float64, (gs.color, gs.color, gs.spaceGrid.size, gs.timeGrid.size))
    green_ins_s = zeros(Float64, (gs.color, gs.color, gs.spaceGrid.size))
    green_dyn_a = zeros(Float64, (ga.color, ga.color, ga.spaceGrid.size, ga.timeGrid.size))
    green_ins_a = zeros(Float64, (ga.color, ga.color, ga.spaceGrid.size))
    for (ki, k) in enumerate(sgrid)
        for (ni, n) in enumerate(gs.dlrGrid.n)
            green_dyn_s[1,1,ki,ni],green_dyn_a[1,1,ki,ni] = KO(k, n, param; pifunc=pifunc, gfactorfunc=gfactorfunc,V_Bare=V_Bare)
        end
        green_ins_s[1,1,ki],green_ins_a[1,1,ki] = V_Bare(k, param), V_Bare(k, param)
    end
    gs.dynamic=green_dyn_s
    gs.instant=green_ins_s
    ga.dynamic=green_dyn_a
    ga.instant=green_ins_a
    return gs, ga
end

end
