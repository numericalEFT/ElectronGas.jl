# Polarization of free electron in two dimensions

We consider the polarziation of free electron gas in 2D, 
$$
\begin{aligned}
\Pi_0(q, \Omega_n)&=-S T\sum_{m}\int\frac{d^2 \vec k}{{(2\pi)}^2}
\frac{1}{(i\omega_m-\varepsilon_{\vec{k}})(i\omega_m+i\Omega_n-\varepsilon_{\vec{k}+\vec{q}})}\\
&=-S\int\frac{d^2 \vec k}{{(2\pi)}^2}
\frac{n(\vec{k})-n(\vec{k}+\vec{q})}{i\Omega_n+\varepsilon_{\vec{k}}-\varepsilon_{\vec{k}+\vec{q}}}\\
\end{aligned}
$$
where $n(\vec k) = 1/(e^{\beta(\epsilon_{\vec k}-\mu)}+1)$, $\epsilon_{\vec k}=k^2/(2m_e)-\mu$, and $S$ is the spin factor. It is expressed as
$$
\Pi_0(q, \Omega_n)=-S\int_0^{\infty} \frac{m_ekdk}{2\pi^2} n(\epsilon_k) \int_{0}^{2\pi} d\theta \left[ \frac{1}{i2m_e\Omega_n-2kq \cos\theta-q^2}-\frac{1}{i2m_e\Omega_n-2kq\cos \theta+q^2}\right]
$$
For static polarization,
	$$
	\Pi_0(q, \Omega_n=0)= -S\int_0^{q/2}\frac{m_ekdk}{2\pi^2} n(\epsilon_k) \frac{4\pi}{\sqrt{q^4-4k^2q^2}}
	$$
- Zero temperature

$$
\Pi_0(q, \Omega_n)=-\frac{Sm_e}{2\pi} \left[1-\frac{2k_F}{q}{\rm Re}\sqrt{\left(\frac{q}{2k_F}+ i \frac{m_e\Omega_n}{qk_F}\right)^2-1} \right]
$$