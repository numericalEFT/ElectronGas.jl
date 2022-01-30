
# Legendre Decomposition of Interaction


## Why decompose?

In many cases we need to calculate integrals of the following form:
```math
   \begin{aligned}
\Delta(\vec{k}) = \int \frac{d^dp}{{2\pi}^d} W(|\vec{k}-\vec{p}|) F(\vec{p}),
   \end{aligned}
```
where only momentum dependency is of our concern.
In such cases we assume that the interaction depends only on the magnitude of
the momentum difference, and the other part of the integrand has some sort of
space symmetry.

Two typical cases are the calculation of self-energy and gap-function equation.
In the case of self-energy, the other part of the integrand depends only on the
magnitude of the momentum, thus we can integrate out the angular dependence to
simplify the calculation. In the case of gap-function equation, we do not know
priorly the symmetry of anomalous correlator, but we can always decompose the
interaction and anomalous correlator into angular momentum channels. In that sense
the self-energy case is actually a special case where we only care about s-channel.


## How?

From now on we illustrate 3D case.

We first express the interaction in Legendre polynomials:
```math
   \begin{aligned}
W(|\vec{k}-\vec{p}|)&=\sum_{l=0}^{\infty}\frac{2l+1}{2} w_l(k, p) p_{l}(\hat{kp}),\\
&=2\pi\sum_{l,m} w_l(k, p) Y_{lm}(\hat{k}/k)Y_{lm}^{*}(\hat{p}/p).
   \end{aligned}
```
with
```math
   \begin{aligned}
w_l(k, p) = \int_{-1}^{1}d\chi P_l(\chi) W(\sqrt{k^2+p^2-2kp\chi}).
   \end{aligned}
```
The other part of the integrand could also be decomposed as
```math
   \begin{aligned}
f_{lm}(k) = \int d\hat{k}Y_{lm}(\hat{k}/k)F(\vec{k})
   \end{aligned}
```
In the cases mentioned above there's no dependency on ``m``, thus the whole calculation
could be decomposed with ``l``.

For gap-function equation we have
```math
   \begin{aligned}
\Delta_l(k) = \int \frac{p^2dp}{{4\pi}^2} w_l(k, p) f_l(p),
   \end{aligned}
```
and for self-energy we have
```math
   \begin{aligned}
\Sigma(k) = \int \frac{p^2dp}{{4\pi}^2} w_0(k, p) G(p).
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
Thus all integral for ``w_l`` could be expressed as combination of helper functions:
```math
   \begin{aligned}
w_0(k,p) &= \frac{1}{kp} \delta H_1(k, p), \\
w_1(k,p) &= \frac{1}{2{(kp)}^2} {[(k^2+p^2)\delta H_1(k, p)-\delta H_3(k, p)]}, \\
w_2(k,p) &= \frac{1}{{(2kp)}^3}
{\{[3{(k^2+p^2)}^2-4k^2p^2]\delta H_1(k, p)-6(k^2+p^2)\delta H_3(k, p)+3\delta H_5(k, p)\} }.
   \end{aligned}
```
