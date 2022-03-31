# Analytic Expression of Moroni's $G_{+}$

## Exchange-correlation Energy

The analytic expression of local field factor is based on the parametrization of correlation energy $\epsilon_c(r_s)$ in Vosko's paper doi: 10.1139/p80-159, equation (4.4):

```math
\epsilon_c(r_s) = A\{n\frac{x^2}{X(x)} + \frac{2b}{Q}tan^{-1}\frac{Q}{2x+b}-\frac{bx_0}{X(x_0)}[ln\frac{(x-x_0)^2}{X(x)}+\frac{2(b+2x_0)}{Q}tan^{-1}\frac{Q}{2x+b}] \}Ry
```

Where $x_0$, $b$,  and $c$ are free parameters obtained by fitting with numerical data, and $x=r_s^{1/2}$, $Q = (4c - b )^{1/2}$, $X(x) = x^2 + bx +c$. 
The parametrization we use for paramagnetic situation are $A = 0.0621814$,  $x_0= -0.10498$,  $b = 3.72744$, $c = 12.9352$, in table 5 of Vosko's paper. 
Note that this parametrization is under the Rydberg energy unit.

## Analytic Expression of $G_{+}$


In 1995 Moroni and others used Diffusion Monte Carlo to produce the local field factor $G_{+}$ (doi: 10.1103/PhysRevLett.75.689). Later in 1998 Corradini and others formulated an analytical expression based on Moroni's numerical data (doi: 10.1103/PhysRevB.57.14569).
The expression is based on three parameters:

```math
\begin{align}
&A =\frac{1}{4}-\frac{k_F^2}{4 \pi e^2}\frac{d \mu_c}{dn_0}\\
&B = \frac{(1+a_1x+a_2 x^2)}{(3+b_1 x +b_2 x^2)}\\
&C =\frac{\pi}{2e^2k_F}\frac{-d(r_s\epsilon_c)}{dr_s} 
\end{align}
```

Here $x=r_s^{1/2}$,  $n_0 = 3/(4\pi a_0^3 r_s^3)$ is the density of homogeneous electron gas at $r_s$, and $\mu_c$ is the contribution of correlation $\epsilon_c$ to the chemical potential:
```math
\mu_c = \frac{d(n_0 \epsilon_c)}{dn_0}
```
The parameters are $a_1 = 2.15, a_2 = 0.435, b_1=1.57, b_2=0.409$, valid for $r_s$ in the range 2-10.
With all parameters ready, the final expression of $G_{+}$ is:

```math
G_{+}(q)=CQ^2+\frac{BQ^2}{g+Q^2}+\alpha Q^4 e^{-\beta Q^2}
```

where 

```math
\begin{align}
&g = \frac{B}{A-C}\\
&\alpha = \frac{1.5}{r_s^{1/4}}\frac{A}{Bg}\\
&\beta = \frac{1.2}{Bg}
\end{align}
```

And $Q = q/k_F$.
