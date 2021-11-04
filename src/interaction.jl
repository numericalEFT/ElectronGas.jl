module Interaction

using ElectronGas.Convention
using Parameters

export RPA, KO

srcdir = "."
rundir = isempty(ARGS) ? "." : (pwd()*"/"*ARGS[1])

# println("Loading parameter from:"*rundir*"/para.jl")
include(rundir*"/para.jl")
using .Para

@unpack me, kF, rs, e0, β , mass2, ϵ0= Para.Param

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

r_s_dl=sqrt(4*0.521*rs/ π );
C1=1-r_s_dl*r_s_dl/4.0*(1+0.07671*r_s_dl*r_s_dl*((1+12.05*r_s_dl)*(1+12.05*r_s_dl)+4.0*4.254/3.0*r_s_dl*r_s_dl*(1+7.0/8.0*12.05*r_s_dl)+1.5*1.363*r_s_dl*r_s_dl*r_s_dl*(1+8.0/9.0*12.05*r_s_dl))/(1+12.05*r_s_dl+4.254*r_s_dl*r_s_dl+1.363*r_s_dl*r_s_dl*r_s_dl)/(1+12.05*r_s_dl+4.254*r_s_dl*r_s_dl+1.363*r_s_dl*r_s_dl*r_s_dl));
C2=1-r_s_dl*r_s_dl/4.0*(1+r_s_dl*r_s_dl/8.0*(log(r_s_dl*r_s_dl/(r_s_dl*r_s_dl+0.990))-(1.122+1.222*r_s_dl*r_s_dl)/(1+0.533*r_s_dl*r_s_dl+0.184*r_s_dl*r_s_dl*r_s_dl*r_s_dl)));
D=inf_sum(r_s_dl,100);
A1=(2.0-C1-C2)/4.0/e0^2*π ;
A2=(C2-C1)/4.0/e0^2*π ;
B1=6*A1/(D+1.0);
B2=2*A2/(1.0-D);

function V_Bare(q)
    e0^2/ϵ0/(q^2+mass2)
end

function Polarization0_ZeroTemp(q, n, β=β)
    Π = 0.0
    x = q/2/kF
    ω_n = 2*π*n/β
    y = me*ω_n/q/kF

    if n == 0
        if abs(q - 2*kF) > EPS
            Π = me*kF/2/π^2*(1 + (1 -x^2)*log1p(4*x/((1-x)^2))/4/x)
        else
            Π = me*kF/2/π^2
        end
    else
        if abs(q - 2*kF) > EPS
            if y^2 < 1e-4/EPS                    
                theta = atan( 2*y/(y^2+x^2-1) )
                if theta < 0
                    theta = theta + π
                end
                @assert theta >= 0 && theta<= π
                Π = me*kF/2/π^2 * (1 + (1 -x^2 + y^2)*log1p(4*x/((1-x)^2+y^2))/4/x - y*theta)
            else
                Π = me*kF/2/π^2 * (2.0/3.0/y^2  - 2.0/5.0/y^4) #+ (6.0 - 14.0*(ω_n/4.0)^2)/21.0/y^6)
            end
        else
            theta = atan( 2/y )
            if theta < 0
                theta = theta + π
            end
            @assert theta >= 0 && theta<= π
            Π = me*kF/2/π^2*(1 + y^2*log1p(4/(y^2))/4 - y*theta)
        end
    end

    return Π
end

function RPA(q, n, β=β)
    kernel = 0.0
    if abs(q) > EPS 
        Π = Polarization0_ZeroTemp(q, n, β)
        if n == 0
            kernel = - Π/( 1.0/V_Bare(q)  + Π )
        else
            Π0 = Π * V_Bare(q)
            kernel = - Π0/( 1.0  + Π0 )
        end
    else
        kernel = 0
    end

    return kernel
end

function KO(q, n, β=β)
    G_s=A1*q^2/(1.0+B1*q^2)+A2*q^2/(1.0+B2*q^2);
    G_a=A1*q^2/(1.0+B1*q^2)-A2*q^2/(1.0+B2*q^2);

    Ks, Ka = 0.0, 0.0

    if abs(q) > EPS 
        Π = Polarization0_ZeroTemp(q, n, β)
        if n == 0
            Ks = - Π*(1-G_s)^2/( 1.0/V_Bare(q)  + Π*(1-G_s))
            Ka = -Π*(-G_a)^2/( 1.0/V_Bare(q) + Π*(-G_a))
        else
            Π0 = Π * V_Bare(q)
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
    println(Interaction.RPA(1.0, 1))
    println(Interaction.KO(1.0, 1))
    println(Interaction.Polarization0_ZeroTemp(1e-5, 1)*1e10)
    println(Interaction.Polarization0_ZeroTemp(1e-7, 1)*1e14)
    println(Interaction.Polarization0_ZeroTemp(1e-9, 1)*1e18)
end
