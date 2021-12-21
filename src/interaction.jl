module Interaction

using Parameters


srcdir = "."
rundir = isempty(ARGS) ? pwd() : (pwd()*"/"*ARGS[1])

include(srcdir*"/parameter.jl")
using .Parameter

include(srcdir*"/convention.jl")
using .Convention

include(srcdir*"/polarization.jl")
using .Polarization

export RPA, KO

srcdir = "."
rundir = isempty(ARGS) ? pwd() : (pwd()*"/"*ARGS[1])

# if !@isdefined Para
#     include(rundir*"/para.jl")
#     using .Para
# end

# println(Parameter.Param)
# @unpack me, kF, rs, e0, beta , mass2, ϵ0= Parameter.Param

function inf_sum(q,n)
    a=q*q
    sum=1.0
    i=0
    j=1.0
    k=2.0
    for i in 1:n
        sum = sum+a/j/k
        a=a*q*q
        j=j*(i+2.0)
        k=k*(i+3.0)
    end
    return 1.0/sum/sum
end




"""
    function V_Bare(q,param)

Bare interaction in momentum space. Coulomb interaction if mass2=0, Yukawa otherwise.

#Arguments:
 - q: momentum
 - param: other system parameters
"""
function V_Bare(q, param)
    @unpack me, kF, rs, e0, beta , mass2, ϵ0 = param
    e0^2/ϵ0/(q^2+mass2)
end


"""
    function RPA(q, n, param)

Dynamic part of RPA interaction, with polarization approximated by zero temperature Π0.

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function RPA(q, n, param)
    @unpack me, kF, rs, e0, beta , mass2, ϵ0= param
    kernel = 0.0
    if abs(q) > EPS 
        Π = Polarization0_ZeroTemp(q, n, param)
        if n == 0
            kernel = - Π/( 1.0/V_Bare(q,param)  + Π )
        else
            Π0 = Π * V_Bare(q,param)
            kernel = - Π0/( 1.0  + Π0 )
        end
    else
        kernel = 0
    end

    return kernel
end

"""
    function GFactorTakada(q, n, param)

G factor with Takada's anzats.

#Arguments:
 - q: momentum
 - n: matsubara frequency given in integer s.t. ωn=2πTn
 - param: other system parameters
"""
function GFactorTakada(q, n, param)
    @unpack me, kF, rs, e0, beta , mass2, ϵ0= param
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
function KO(q, n, param)
    @unpack me, kF, rs, e0, beta , mass2, ϵ0= param

    G_s, G_a = GFactorTakada(q, n, param)
    Ks, Ka = 0.0, 0.0

    if abs(q) > EPS 
        Π = Polarization0_ZeroTemp(q, n, param)
        if n == 0
            Ks = - Π*(1-G_s)^2/( 1.0/V_Bare(q, param)  + Π*(1-G_s))
            Ka = -Π*(-G_a)^2/( 1.0/V_Bare(q,param) + Π*(-G_a))
        else
            Π0 = Π * V_Bare(q,param)
            Ks = - Π0*(1-G_s)^2/( 1.0  + Π0*(1-G_s))
            Ka = -Π0*(-G_a)^2/(1.0 + Π0*(-G_a))
        end
    else
        Ks, Ka = 0.0, 0.0
    end

    return Ks, Ka
end

end


if abspath(PROGRAM_FILE) == @__FILE__
    println(Interaction.RPA(1.0, 1, Interaction.Parameter.Param))
    println(Interaction.KO(1.0, 1, Interaction.Parameter.Param))
end
