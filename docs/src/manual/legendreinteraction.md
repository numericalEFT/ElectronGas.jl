
# Legendre Decomposition of Interaction


## Why decompose?

In many cases we need to calculate integrals of the following form:
```math
   \begin{aligned}
\Delta(\vec{k}) = \int \frac{d^dp}{{2\pi}^d} W(|\vec{k}-\vec{p}|) F(\vec{p}),
\tag{1}
   \end{aligned}
```
where only momentum dependency is of our concern. In such cases we assume that the interaction depends only on the magnitude of the momentum difference, and the other part of the integrand has some sort of space symmetry.

Two typical cases are the calculation of self-energy and gap-function equation. In the case of self-energy, the other part of the integrand depends only on the magnitude of the momentum, thus we can integrate out the angular dependence to simplify the calculation. In the case of gap-function equation, we do not know priorly the symmetry of anomalous correlator, but we can always decompose the interaction and anomalous correlator into angular momentum channels. In that sense, the self-energy case is actually a special case where we only care about s-channel.


## How?
We first express the interaction in Legendre polynomials:
```math
   \begin{aligned}
W(|\vec{k}-\vec{p}|)&=\sum_{\ell=0}^{\infty}\frac{N(d,\ell)}{2} w_l(k, p) P_{l}(\hat{kp}) \,,\\
w_{\ell}(k, p) &= \int_{-1}^{1}d\chi P_l(\chi) W(\sqrt{k^2+p^2-2kp\chi}) \,, \\
   \end{aligned}
```
where 
```math
N(d, \ell)=\frac{2 \ell+d-2}{\ell}\left(
\begin{array}{c}
\ell+d-3 \\
\ell-1
\end{array}\right)
```
denotes the number of linearly independent Legendre polynomials of degree ``\ell`` in ``d`` dimensions.
The Legendre polynomials of a scalar product of unit vectors can be expanded using the addition theorem for spherical harmonics
```math
P_{\ell}(\hat{k p})=\frac{\Omega_{d}}{N(d,\ell)} \sum_{m=1}^{N(d,\ell)} Y_{\ell m}(\hat{k}) Y_{\ell m}^{*}(\hat{p})\,,
```
where ``\Omega_{d}=(2\pi)^{\frac d 2}/\Gamma(\frac d 2)`` is the solid angle in ``d`` dimensions. The spherical harmonics are orthonormal as 
```math
\int {\rm d}\hat k Y_{\ell m}(\hat k)  Y^*_{\ell^\prime m^\prime}(\hat k)  = \delta_{\ell \ell^\prime} \delta_{mm^\prime} \,.
```
Hence, the ``W(|\vec{k}-\vec{p}|)`` function can expressed further on as
```math
W(|\vec{k}-\vec{p}|)=\frac{\Omega_d}{2} \sum_{\ell m}  w_l(k, p) Y_{lm}(\hat{k})Y_{lm}^{*}(\hat{p}).
```

The other part of the integrand could also be decomposed as
```math
   \begin{aligned}   
   F(\vec p) &= \sum_{\ell m} f_{\ell m}(p) Y_{\ell m}(\hat p)\,, \\
f_{lm}(p) &= \int d\hat{k}Y^*_{lm}(\hat{p})F(\vec{p})
   \end{aligned}
```

Projecting Eq.(1), on the spherical harmonic ``Y_{\ell m}(\hat k)``, we have
```math
\begin{aligned}
\sum_{\ell m} \Delta_{\ell m }(k)Y_{\ell m}(\hat k) &= \frac{\Omega_d}{2} \int \frac{{\rm d}\vec p}{(2\pi)^d} \sum_{\ell m} w_{\ell}(k,p)   Y_{\ell m}(\hat{k}) Y_{\ell m}^{*}(\hat{p} ) \sum_{\ell^\prime m^\prime}f_{\ell^\prime m^\prime}(p) Y_{\ell^\prime m^\prime}(\hat p)  \\
&=\frac{\Omega_d}{2} \sum_{\ell m} \int \frac{p^{d-1}dp}{(2\pi)^d}  w_{\ell}(k,p,\tau) f_{\ell m}(p, \tau) Y_{\ell m}(\hat k) \,,
\end{aligned}
```
which leads to the decoupled equations with channels.

For gap-function equation, we have
```math
   \begin{aligned}
\Delta_l(k) = \frac{\Omega_d}{2}  \int \frac{p^{d-1} {\rm d}p}{(2\pi)^d} w_l(k, p) f_l(p) \,.
   \end{aligned}
```
For ``GW`` self-energy symmertric with ``\hat k``, we have
```math
   \begin{aligned}
\Sigma(k) = \frac{\Omega_d}{2}  \int \frac{p^{d-1} {\rm d}p}{(2\pi)^d} w_0(k, p) G(p).
   \end{aligned}
```

## Helper function

We can do the decomposition with helper function:
```math
\begin{aligned}
H_n(y) = \int_{0}^{y} z^n W(z)dz.
\end{aligned}
```
Then
```math
\begin{aligned}
\int_{|k-p|}^{k+p} z^n W(z)dz = H_n(k+p)-H_n(|k-p|)= \delta H_n(k, p).
\end{aligned}
```
Changing variables in Eq.(1) to ``z^2=k^2+p^2-2kp\chi``, we obatin 
```math
w_{\ell}(k, p)=\frac{1}{k p} \int_{|k-p|}^{k+p} z d z P_{\ell}\left(\frac{k^{2}+p^{2}-z^{2}}{2 k p}\right) W(z)
```
Since ``P_{\ell}(x)`` is a polynomial in ``x``, for any ``k``,``p``, and integer ``\ell``, ``w_l`` can be expressed as the combination of helper functions.

For electron gas, the ``W`` function contains two terms, bare interaction ``V`` and generic ``W(q,\tau)``. We only tabulated helper functions for the second term. Helper functions for the first term ``h_n`` can be done analytically.

## Three dimensions
For 3D, ``N(3,\ell)=2\ell +1``, and ``Y_{\ell m}`` is the standard sphereical harmonic function. The ``GW`` self energy is
```math
   \begin{aligned}
\Sigma(k) =   \int \frac{p^2 {\rm d}p}{4\pi^2} w_0(k, p) G(p) \,.
   \end{aligned}
```
For 3D electron gas with bare interaction ``V=\frac{4\pi e^2}{q^2} \delta(\tau)\, \left( V(r)=\frac{e^2}{r} \right)``, the helper functions for the first term have
```math
\delta h_{1}(k, p)=4 \pi e^{2} \ln \frac{k+p}{|k-p|}, \quad \delta h_{3}(k, p)=4 \pi e^{2}[2 k p], \quad \delta h_{5}(k, p)=4 \pi e^{2}\left[2 k p\left(k^{2}+p^{2}\right)\right] .
```
3D Yukawa interaction has
```math
\begin{aligned}
V(r) &= \frac{e^2}{r}e^{-mr}, \quad V(q)=\frac{4\pi e^2}{q^2+m^2} ,\\
\delta h_{1}(k, p)&=2 \pi e^{2} \ln\left[\frac{(k+p)^2+m^2}{(k-p)^2+m^2}\right], \\
\delta h_{3}(k, p)&=4 \pi e^{2} \left[2kp - \frac{m^2}{2} \ln \frac{(k+p)^2+m^2}{(k-p)^2+m^2} \right], \\
\delta h_{5}(k, p)&=4 \pi e^{2}\left[2 k p\left(k^2+p^2 -m^2\right) +\frac{m^4}{2}\ln \frac{(k+p)^2+m^2}{(k-p)^2+m^2} \right].
\end{aligned}
```

``w_l`` can be expressed as
```math
   \begin{aligned}
w_0(k,p) &= \frac{1}{kp} \delta H_1(k, p), \\
w_1(k,p) &= \frac{1}{2{(kp)}^2} {[(k^2+p^2)\delta H_1(k, p)-\delta H_3(k, p)]}, \\
w_2(k,p) &= \frac{1}{{(2kp)}^3}
{\{[3{(k^2+p^2)}^2-4k^2p^2]\delta H_1(k, p)-6(k^2+p^2)\delta H_3(k, p)+3\delta H_5(k, p)\} }.
   \end{aligned}
```

## Two dimensions
For 2D, the sphereical harmonic is just a Fourier series as
```math
\begin{aligned}
 P_{\ell}(\hat{kp}) &=  \pi \cos[\ell(\theta_{\hat k} - \theta_{\hat p})] , \quad N(2,\ell) = 2\\
Y_{\ell 1}(\hat k) &= \cos(\ell \theta),\;  Y_{\ell 2}(\hat k) = \sin(\ell \theta).
\end{aligned}
```
The ``GW`` self energy is
```math
   \begin{aligned}
\Sigma(k) =   \int \frac{p {\rm d}p}{4\pi} w_0(k, p) G(p) \,.
   \end{aligned}
```
For 2D electron gas with 2D coulomb potential ``V=\frac{2\pi e^2}{q} \delta(\tau)\, \left( V(r)=\frac{e^2}{r} \right)``, 
the helper functions for the first term have
```math
\delta h_{1}(k, p)=2\pi e^{2} (k+p-|k-p|), \quad \delta h_{3}(k, p)=\frac 2 3 \pi e^{2}\left [(k+p)^3-|k-p|^3\right ].
```
For ``V=\frac{4\pi e^2}{q^2} \delta(\tau)\,\left(V(r)=-\ln \frac{r}{L}\right)``, ``\delta h_n`` has the same form as in 3D electron gas. 

2D Yukawa interaction has
```math
\begin{aligned}
V(r)&=\frac{e^2}{r} e^{-mr} ,\\
V(q)&=\int {\rm d}^2\vec r V(r) e^{i\vec q\cdot \vec r} \\
    & = \int r dr \frac{e^2}{r} e^{-mr} \int_0^{2\pi} d\theta e^{iqr\cos \theta} \\
    & = \int dr e^2 e^{-mr} 2\pi J_0(qr) \\
    & = \frac{2\pi e^2}{\sqrt{q^2+m^2}} \,,
\end{aligned}
```
and the helper functions for the first term have
```math
\begin{aligned}
\delta h_{1}(k, p)&=2\pi e^{2} (\sqrt{(k+p)^2+m^2}-\sqrt{(k-p)^2+m^2}), \\
\delta h_{3}(k, p)&=\frac 2 3 \pi e^{2}\left ([(k+p)^2-2m^2]\sqrt{(k+p)^2+m^2} - [(k-p)^2-2m^2]\sqrt{(k-p)^2+m^2} \right ).
\end{aligned}
```

**Reference**:
1. Christopher Frye and Costas J. Efthimiou, Spherical Harmonics in ``p`` Dimensions, arXiv 1205.3548