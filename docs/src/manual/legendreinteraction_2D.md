# Decomposition of Interaction

In GW-approximation, we calculation self-energy as
$$
\Sigma(\mathbf{k},\omega_n)=-T\int \frac{{\rm d}^d \mathbf{q}}{(2\pi)^d} \sum_m G(\mathbf{p},\omega_m)W(\mathbf{k-p},\omega_n-\omega_m) \,,
\tag{1} 
$$
where $G$ is the Green's function and W is the effective interaction. Here, we suppress spin index.

## Spherical harmonic representation
We first express the $W(q,\tau)$ function as an expansion in Legendre polynomials $P_\ell(\chi)$
$$
\begin{gathered}
W(|\mathbf{k}-\mathbf{p}|, \tau)=\sum_{\ell=0}^{\infty} \bar{w}_{\ell}(k, p, \tau) P_{\ell}(\hat{k p}) \,, \\
\bar{w}_{\ell}(k, p, \tau)=\frac{N(d,\ell)}{2} \int_{-1}^{1} d \chi P_{\ell}(\chi) W\left(\sqrt{k^{2}+p^{2}-2 k p \chi} ,\tau\right)\,.
\end{gathered}
$$
Since the Legendre polynomials of a scalar product of unit vectors can be expanded with spherical harmonics using
$$
P_{\ell}(\hat{k p})=\frac{\Omega_{d}}{N(d,\ell)} \sum_{m=1}^{N(d,\ell)} Y_{\ell m}(\hat{k}) Y_{\ell m}^{*}(\hat{p})\,,
$$
where $\Omega_{d}$ is the solid angle in $d$ dimensions, and 
$$
N(d, \ell)=\frac{2 \ell+d-2}{\ell}\left(
\begin{array}{c}
\ell+d-3 \\
\ell-1
\end{array}\right)
$$
denotes the number of linearly independent homogeneous harmonic polynomials of degree $\ell$ in $d$ dimensions. The spherical harmonics are orthonormal as 
$$
\int {\rm d}\Omega_{\hat k} Y_{\ell m}(\hat k)  Y_{\ell^\prime m^\prime}(\hat k)  = \delta_{\ell \ell^\prime} \delta_{mm^\prime} \,.
$$
Hence, the $W(q,\tau)$ function can expressed further on as
$$
W(|\mathbf{k}-\mathbf{p}|, \tau)=\sum_{\ell} \frac{\Omega_{d}}{N(d,\ell)}  \bar{w}_{\ell}(k, p, \tau) \sum_{m} Y_{\ell m}(\hat{k}) Y_{\ell m}^{*}(\hat{p})
$$

or 
$$
W(|\mathbf{k}-\mathbf{p}|, \tau)= \frac{\Omega_{d}}{2}\sum_{\ell}  w_{\ell}(k, p, \tau) \sum_{m} Y_{\ell m}(\hat{k}) Y_{\ell m}^{*}(\hat{p} )\,.
$$
with 
$$
w_{\ell}(k, p, \tau)=\int_{-1}^{1} d \chi P_{\ell}(\chi) W\left(\sqrt{k^{2}+p^{2}-2 k p \chi} ,\tau\right)
$$
In addition, the Green's function $G(\mathbf p, \tau)$ has 
$$
\begin{aligned}
	G(\mathbf p, \tau)= \sum_{\ell=0}^{\infty}\sum_{m=1}^{N(d,\ell)} G_{\ell m}(p,\tau) Y_{\ell m}(\hat p)\,, \\
	G_{\ell m}(p,\tau) = \int {\rm d}\Omega_{\hat k} G(\mathbf p,\tau) Y^*_{\ell m}(\hat p)
\end{aligned}
$$
### Decouple with channels
By the sphereical harmonic expansion of Eq.(1), the self-energy 

$$
\begin{aligned}
\sum_{\ell m} \Sigma_{\ell m }(k,\tau)Y_{\ell m}(\hat k) &= \frac{\Omega_d}{2} \int \frac{{\rm d}\mathbf p}{(2\pi)^d} \sum_{\ell m}G_{\ell m}(p, \tau) Y_{\ell m}(\hat p) \sum_{\ell^\prime m^\prime} w_{\ell^\prime}(k,p,\tau)   Y_{\ell^\prime m^\prime}(\hat{k}) Y_{\ell^\prime m^\prime}^{*}(\hat{p} ) \\
&=\frac{\Omega_d}{2} \sum_{\ell m} \int \frac{p^{d-1}dp}{(2\pi)^d} G_{\ell m}(p, \tau) w_{\ell}(k,p,\tau) Y_{\ell m}(\hat k)
\end{aligned}
$$
Since self-energy is symmertric with $\hat k$, we just need to project on the s-wave channel, namely
$$
\Sigma(k,\tau) =\frac{\Omega_d}{2}  \int \frac{p^{d-1}dp}{(2\pi)^d} G(p,\tau) w_0(k,p,\tau)
$$
###  Two dimensions
$$
\begin{aligned}
N(2,\ell) &= 2 \\
Y_{\ell 1}(\hat k) &= \cos(\ell \theta),\;  Y_{\ell 2}(\hat k) = \sin(\ell \theta) \\
P_{\ell}(\hat{kp}) &=  \pi \cos[\ell(\theta_{\hat k} - \theta_{\hat p})]
\end{aligned}
$$
Hence, the self-energy is 
$$
\Sigma(k,\tau) = \int \frac{p^{d-1}dp}{4\pi} G(p,\tau)w_0(k,p,\tau)
$$

## Helper function  

We can do the decomposition with helper function:
$$
H_n(y) = \int_{0}^{y} z^n W(z)dz\,.
$$
Then
$$
\begin{aligned}
\int_{|k-p|}^{k+p} z^n W(z)dz = H_n(k+p)-H_n(|k-p|)= \delta H_n(k, p)\,,\\
w_{\ell}(k, p)=\frac{1}{k p} \int_{|k-p|}^{k+p} z d z P_{\ell}\left(\frac{k^{2}+p^{2}-z^{2}}{2 k p}\right) W(z)
\end{aligned}
$$
Thus all integral for $w_{\ell}$ could be expressed as combination of helper functions:
$$
\begin{aligned}
w_0(k,p) &= \frac{1}{kp} \delta H_1(k, p), \\
w_1(k,p) &= \frac{1}{2{(kp)}^2} {[(k^2+p^2)\delta H_1(k, p)-\delta H_3(k, p)]}, \\
w_2(k,p) &= \frac{1}{{(2kp)}^3}
{\{[3{(k^2+p^2)}^2-4k^2p^2]\delta H_1(k, p)-6(k^2+p^2)\delta H_3(k, p)+3\delta H_5(k, p)\} }.
 \end{aligned}
 $$