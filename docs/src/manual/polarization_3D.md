
# Polarization of free electron

## Generic formalism

The bare polarization is defined as

```math
\begin{aligned}
\Pi_0(\Omega_n, \vec{q})=
-S T\sum_{m}\int\frac{d^d k}{{(2\pi)}^d}G_0(\omega_m, \vec{k})G_0(\omega_m+\Omega_n, \vec{k}+\vec{q})
\end{aligned}
```
in the matsubara representation, where ``S`` is spin number, ``G_0`` is bare Green's function.
We have

```math
\begin{aligned}
G_0(\omega_m, \vec{k})=\frac{1}{i\omega_m-\epsilon_{\vec{k}}}
\end{aligned}
```

with bare electron dispersion given by ``\epsilon_{\vec{k}}``. By summing over frequency we have 

```math
\begin{aligned}
\Pi_0(\Omega_n, \vec{q}) &= -S \int\frac{d^d \vec k}{{(2\pi)}^d}
\frac{n(\epsilon_{\vec{k}+\vec{q}})-n(\epsilon_{\vec{k}})}{i\omega_n-\epsilon_{\vec{k}+\vec{q}}+\epsilon_{\vec{k}}} \\
&=-S\int \frac{d^d \vec k}{{(2\pi)}^d} n(\epsilon_{\vec k}) \left[ \frac{1}{i\Omega+\epsilon_{\vec k}-\epsilon_{\vec k+\vec q}}-\frac{1}{i\Omega+\epsilon_{\vec k-\vec q}-\epsilon_{\vec k}}\right]
\end{aligned}
```
with the fermi distribution function 

```math
n(\epsilon_{\vec k}) =\frac{1}{e^{\beta\varepsilon_{\vec{k}}}+1}
```

## Free electron in 3D

From now on we consider free electron in 3D, where ``d=3`` and dispersion ``\epsilon_{\vec{k}}=k^2/2m-\mu``. We have

```math
\begin{aligned}
\Pi_0(\Omega, \vec{q})&=-S\int_0^{\infty} \frac{k^2dk}{4\pi^2} n(\epsilon_k) \int_{-1}^{1} d(\cos \theta) \left[ \frac{1}{i\Omega+\epsilon_{\vec k}-\epsilon_{\vec k+\vec q}}-\frac{1}{i\Omega+\epsilon_{\vec k-\vec q}-\epsilon_{\vec k}}\right]\\
&=-S\int_0^{\infty} \frac{k^2dk}{4\pi^2} n(\epsilon_k) \frac{m}{kq}\ln\frac{4m^2\Omega^2+(q^2-2kq)^2}{4m^2\Omega^2+(q^2+2kq)^2} \,,
\end{aligned}
```
which could be handled with one dimensional integral of ``k``.

- In the limit ``q^2+2k_F q \ll 2m\Omega_n ``, the intergrand of ``\Pi_0`` is expanded as 
  ```math
  \frac{m}{kq}\ln\frac{4m^2\Omega^2+(q^2-2kq)^2}{4m^2\Omega^2+(q^2+2kq)^2}=-\frac{2q^2}{m\Omega^2}+\frac{2k^2q^4}{m^3\Omega^4}+\frac{(-4k^2+m^2\Omega^2)q^6}{2m^5\Omega^6}+...
  ```
- Zero temperature polarization can be calculated explicitly
  ```math
  \Pi_0(\Omega,q) = -\frac{N_F}{2}\left[1-\frac{1}{8 k_F q}\left\{ \left[\frac{(i2m\Omega-q^2)^2}{q^2}-4 k_F^2\right]\log\left(\frac{i2m\Omega-q^2-2 k_F q}{i2m\Omega-q^2+2 k_F q}\right)+\left[\frac{(i2m\Omega+q^2)^2}{q^2}-4 k_F^2\right]\log\left(\frac{i2m\Omega+q^2+2 k_F q}{i2m\Omega+q^2-2 k_F q}\right)\right\}\right]
  ```
- In the static limit ``\Omega=0``, 
  ```math 
  \Pi_0(0, q) = -N_F F(q/2k_F) \,,
  ```
   where ``N_F=Smk_F/(2\pi^2)`` is the density of states, and ``F(x)=\frac{1}{2}-\frac{x^2-1}{4x}\ln \left|\frac{1+x}{1-x}\right|`` is the Lindhard function. 
   
   The weak logarithmic singularity near ``2k_F`` is the cause of the Friedel oscillation and Kohn-Luttinger superconductivity.
  
## Polarization in the large frequency limit ``\Omega \gg q v_F``
As derived in [[polarization (free electron)#Generic formalism]]
- In the Matsubara representation,
```math
P_{q, \Omega}=S\int \frac{d^D k}{(2\pi)^D} \frac{n(\epsilon_k)-n(\epsilon_{k+q})}{i\Omega+\epsilon_k-\epsilon_{k+q}}
```

```math
P_{q, \Omega}=S\int \frac{d^D k}{(2\pi)^D} n(\epsilon_k) \left[ \frac{1}{i\Omega+\epsilon_k-\epsilon_{k+q}}-\frac{1}{i\Omega+\epsilon_{k-q}-\epsilon_{k}}\right]
```

Consider the limit ``\Omega \gg q v_F``,

```math
P_{q, \Omega}=\frac{S}{i\Omega}\int \frac{d^D k}{(2\pi)^D} n(\epsilon_k) \left[ \frac{1}{1-\Lambda_q/i\Omega}-\frac{1}{1+\Lambda_q/i\Omega}\right]
```

where ``\Lambda_q=\epsilon_{k+q}-\epsilon_k=(2k \cdot q+q^2)/2m``
   
```math
P_{q, \Omega}=\frac{2S}{(i\Omega)^2}\int \frac{d^D k}{(2\pi)^D} n(\epsilon_k) \left[\Lambda_q+\Lambda_q^3/(i\Omega)^2+...\right]=\frac{n}{m}\frac{q^2}{(i\Omega)^2}\left[1+O\left(\frac{q^2}{(i\Omega)^2}\right)\right]
```

where is exact in arbitrary dimensions and the electron density ``n=S\int \frac{d^D k}{(2\pi)^D} n(\epsilon_k)``.
 
- The correction term is ``\left[\frac{3}{5}(q v_F)^2+\epsilon_q^2\right]/(i\Omega)^2`` for 3D.
 
 In real frequency, 

```math
\operatorname{Re} P_{q, \omega} = \frac{n}{m}\frac{q^2}{\omega^2}\left[1+O\left(\frac{q^2}{\omega^2}\right)\right]
```

### Plasmon frequency
Plasmon dispersion is the zeros of the dynamic dielectric function,
```math
\epsilon=1-v_q P_{q, \omega}=0
```
where ``v_q=4\pi e^2/q^2``. 
The dispersion of the plasmon is,
```math
\omega_p^2=\frac{4\pi e^2n}{m}\left[1+O\left(\frac{q^2}{\omega^2}\right)\right]
```

Plasmon emerges only if ``\epsilon_q \cdot v_q \sim \text{constant}``.
