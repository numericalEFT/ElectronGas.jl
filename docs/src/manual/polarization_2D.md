# Polarization of free electron in two dimensions

We consider the polarziation of free electron gas in 2D, 
```math
\Pi_0(q, \Omega_n)
=-S\int\frac{d^2 \vec k}{{(2\pi)}^2}
\frac{n(\epsilon_{\vec{k}})-n(\epsilon_{\vec{k}+\vec{q}})}{i\Omega_n+\epsilon_{\vec{k}}-\epsilon_{\vec{k}+\vec{q}}}\\
```
where ``n(\vec k) = 1/(e^{\beta(\epsilon_{\vec k}-\mu)}+1)``, ``\epsilon_{\vec k}=k^2/(2m)-\mu``, and ``S`` is the spin factor. It is expressed as
```math
\Pi_0(q, \Omega_n)=-S\int_0^{\infty} \frac{mkdk}{2\pi^2} n(\epsilon_k) \int_{0}^{2\pi} d\theta \left[ \frac{1}{i2m\Omega_n-2kq \cos\theta-q^2}-\frac{1}{i2m\Omega_n-2kq\cos \theta+q^2}\right]
```

## Static limit ``\Omega_n=0`` 

For real ``a,b``, the integral ``\int_0^{2\pi}\frac{1}{a-b\cos \theta}=\frac{2\pi}{\sqrt{a^2-b^2}}`` and ``\int_0^{2\pi}\frac{1}{a+b\cos \theta}=C\frac{2\pi}{\sqrt{a^2-b^2}}`` where ``C=1`` for ``a>b`` and ``C=-1`` for ``a<b``. Hence, we have 	
```math
\Pi_0(q, 0)= -S\int_0^{q/2}\frac{mkdk}{2\pi^2} n(\epsilon_k) \frac{4\pi}{\sqrt{q^4-4k^2q^2}}
```

At zero temperature, 
- for ``q<2k_F``, 
  ```math
  \Pi_0(q, 0)=-S\int_0^{q/2} \frac{m}{\pi q} \frac{dk^2}{\sqrt{q^2-4k^2}}=-\int_0^1\frac{mS}{4\pi}\frac{dx}{\sqrt{1-x}}=-\frac{mS}{2\pi} \,;
  ```
- for  ``q>2k_F``, 
  ```math 
  \Pi_0(q, 0)=-S\int_0^{k_F} \frac{m}{\pi q} \frac{dk^2}{\sqrt{q^2-4k^2}}=-\int_0^{4k_F^2/q^2}\frac{mS}{4\pi}\frac{dx}{\sqrt{1-x}}=-\frac{mS}{2\pi}\left( 1-\sqrt{1-\frac{4k_F^2}{q^2}}\right) \,.
  ```

## Zero temperature
```math
\Pi_0(q, \Omega_n)=-\frac{mS}{2\pi} \left[1-\frac{2k_F}{q}{\rm Re}\sqrt{\left(\frac{q}{2k_F}+ i \frac{m\Omega_n}{qk_F}\right)^2-1} \right] \,.
```
For ``q\to 0`` and ``\Omega_n \neq 0``, ``\Pi_0(q, \Omega_n) \to 0``.