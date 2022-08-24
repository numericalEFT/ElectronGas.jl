# Polarization Approximation

In this note, we discuss serveral approximations of the polarization

## Linearized Dispersion
The original polarization is given by,
```math
\begin{aligned}
\Pi_0(\Omega_n, \vec{q}) &= -S \int\frac{d^d \vec k}{{(2\pi)}^d}
\frac{n(\epsilon_{\vec{k}+\vec{q}})-n(\epsilon_{\vec{k}})}{i\omega_n-\epsilon_{\vec{k}+\vec{q}}+\epsilon_{\vec{k}}} \\
&=-S\int \frac{d^d \vec k}{{(2\pi)}^d} n(\epsilon_{\vec k}) \left[ \frac{1}{i\Omega+\epsilon_{\vec k}-\epsilon_{\vec k+\vec q}}-\frac{1}{i\Omega+\epsilon_{\vec k-\vec q}-\epsilon_{\vec k}}\right].
\end{aligned}
```
One possible approximation is to replace the kinetic energy with a dispersion linearized near the Fermi surface, 
```math
\xi_{\mathbf{p}+\mathbf{q}}-\xi_{\mathbf{p}}=(1 / m) \mathbf{p} \cdot \mathbf{q}+\mathcal{O}\left(q^{2}\right)
```
so that,
```math
n_{\mathrm{F}}\left(\epsilon_{\mathbf{p}+\mathbf{q}}\right)-n_{\mathrm{F}}\left(\epsilon_{\mathbf{p}}\right) \simeq \partial_{\epsilon_p} n_{\mathrm{F}}\left(\epsilon_{p}\right)(1 / m) \mathbf{p} \cdot \mathbf{q} \simeq-\delta\left(\epsilon_{p}-\mu\right)(1 / m) \mathbf{p} \cdot \mathbf{q}
```
where, in the zero-temperature limit, the last equality becomes exact. Converting the momentum sum into an integral, we thus obtain
```math
\Pi_0(\mathbf{q}, \omega_{m})=-S \int \frac{d^{3} p}{(2 \pi)^{3}} \delta\left(\epsilon_{p}-\mu\right) \frac{\frac{1}{m} \mathbf{p} \cdot \mathbf{q}}{i \omega_{m}+\frac{1}{m} \mathbf{p} \cdot \mathbf{q}}.
```
Evaluate the integral gives
```math
\begin{aligned}
\Pi_0(\mathbf{q}, \omega_{m}) &=-\frac{S}{(2 \pi)^{3}} \int d p p^{2} \int d \delta\left(\epsilon_{p}-\mu\right) \frac{v_{\mathrm{F}} \mathbf{n} \cdot \mathbf{q}}{i \omega_{m}+v_{\mathrm{F}} \mathbf{n} \cdot \mathbf{q}} \\
&=-\underbrace{\frac{S}{(2 \pi)^{3}} \int d p p^{2} \int d \delta\left(\epsilon_{p}-\mu\right)}_{N_F} \frac{1}{\int d\Omega} \int d\Omega \frac{v_{\mathrm{F}} \mathbf{n} \cdot \mathbf{q}}{i \omega_{m}+v_{\mathrm{F}} \mathbf{n} \cdot \mathbf{q}} \\
&=-\frac{N_F}{2} \int_{-1}^{1} d x \frac{v_{\mathrm{F}} x q}{i \omega_{m}+v_{\mathrm{F}} x q}=-N_F\left[1-\frac{i \omega_{m}}{2 v_{\mathrm{F}} q} \ln \left(\frac{i \omega_{m}+v_{\mathrm{F}} q}{i \omega_{m}-v_{\mathrm{F}} q}\right)\right] .
\end{aligned}
```
The above derivation is adapted from the A. Altland and B. Simons' book "Condensed Matter Field Theory" Chapter 5.2, Eq. (5.30).
### Two limits:

For the exact free-electron polarization, we expect
In the limit ``q ≫ ω_m``,
```math
Π_0(q, iω_m) \rightarrow -N_F \left(1-\frac{π}{2}\frac{|ω_m|}{v_{\mathrm{F}} q}\right)
```
where we use the Taylor expansion for ``\text{Log}\left[\frac{1+i x}{-1+i x}\right]`` where ``x=\omega_m/(v_{\mathrm{F}} q)``,
```math
\begin{array}{cc}
 \{ & 
\begin{array}{cc}
 -i \pi +2 i x-\frac{2 i x^3}{3}+O\left(x^4\right) & x \ge 0 \\
 i \pi +2 i x-\frac{2 i x^3}{3}+O\left(x^4\right) & x<0 \\
\end{array}
 \\
\end{array}
```

and in the limit ``q ≪ ω_m``, 
```math
Π_0(q, iω_m) \rightarrow -\frac{N_F}{3}\left(\frac{v_{\mathrm{F}} q}{ω_m}\right)^2 = -N_F \left(\frac{q}{q_{\mathrm{TF}}}\frac{\omega_p}{\omega_m}\right)^2
```
where the plasma-frequency and the Thomas-Fermi screening momentum is related by ``ω_p=v_F q_{\mathrm{TF}}/\sqrt{3}``.


## Plasma Approximation 

It is sometimes convenient to approximate the polarization with the plasma poles,
```math
Π_0(q, iω_m) \approx -N_F \frac{(q/q_{\mathrm{TF}})^2}{(q/q_{\mathrm{TF}})^2+(\omega_m/\omega_p)^2}
```