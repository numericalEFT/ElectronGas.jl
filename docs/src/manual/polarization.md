
# Polarization of free electron

## Definition

The bare polarization is defined as

```math
\begin{aligned}
\Pi_0(\omega_n, \vec{q})=
N_ST\sum_{m}\int\frac{d^Dk}{{(2\pi)}^D}G_0(\omega_m, \vec{k})G_0(\omega_{m+n}, \vec{k}+\vec{q})
\end{aligned}
```

where ``N_S`` is spin number, ``G_0`` is bare Green's function.
We have

```math
\begin{aligned}
G_0(\omega_m, \vec{k})=\frac{1}{i\omega_m-\varepsilon_{\vec{k}}+\mu}
\end{aligned}
```

with dispersion given by ``\varepsilon_{\vec{k}}`` and chemical potential ``\mu``.


## Free electron in 3D

   From now on we consider free electron in 3D, where ``D=3``, ``N_S=2`` and ``\varepsilon_{\vec{k}}=\frac{k^2}{2m}``.
Then by summing over frequency we have

```math
\begin{aligned}
\Pi_0(\omega_n, \vec{q})&=2T\sum_{m}\int\frac{d^3k}{{(2\pi)}^3}
\frac{1}{(i\omega_m-\varepsilon_{\vec{k}}+\mu)(i\omega_{m+n}-\varepsilon_{\vec{k}+\vec{q}}+\mu)}\\
&=-2\int\frac{d^3k}{{(2\pi)}^3}
\frac{n(\vec{k}+\vec{q})-n(\vec{k})}{i\omega_n-\varepsilon_{\vec{k}+\vec{q}}+\varepsilon_{\vec{k}}}\\
\end{aligned}
```

with

```math
\begin{aligned}
n(\vec{k})=\frac{1}{e^{\beta(\varepsilon_{\vec{k}}-\mu)}+1}.
\end{aligned}
```

Now plug in the dispersion and we have

```math
\begin{aligned}
\Pi_0(\omega_n, \vec{q})
&=-2\int\frac{d^3k}{{(2\pi)}^3}
[\frac{n(\vec{k})}{i\omega_n-\varepsilon_{\vec{k}}+\varepsilon_{\vec{k}+\vec{q}}}
-\frac{n(\vec{k})}{i\omega_n-\varepsilon_{\vec{k}+\vec{q}}+\varepsilon_{\vec{k}}}
]\\
&=-2\int\frac{d^3k}{{(2\pi)}^3}
\frac{4mn(\vec{k})(q^2+2kq\cos(\theta))}
{4m^2\omega_n^2+{(q^2+2kq\cos(\theta))}^2}\\
&=-\int\frac{kmdk}{2\pi^2q}
n(\vec{k})\ln(\frac{4m^2\omega_n^2+{(q^2+2kq)}^2}{4m^2\omega_n^2+{(q^2-2kq)}^2})\\
\end{aligned}
```

which could be handled with one dimensional integral of ``k``.

