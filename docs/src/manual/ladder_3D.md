# Ladder (Particle-particle bubble) of free electrons

```math
\begin{split}
 &\int \frac{d^3 \vec{p}}{\left(2\pi^3\right)} T \sum_{i \omega_n} \frac{1}{i \omega_n+i \Omega_n-\frac{(\vec{k}+\vec{p})^2}{2 m}+\mu} \frac{1}{-i \omega_n-\frac{p^2}{2 m}+\mu} \\
 =&\int \frac{d^3 \vec{p}}{\left(2\pi^3\right)}  \frac{f\left(\frac{(\vec{k}+\vec{p})^2}{2 m}-\mu\right)-f\left(-\frac{p^2}{2 m}+\mu\right)}{i \Omega_n-\frac{(k+\vec{p})^2}{2 m}-\frac{p^2}{2 m}+2 \mu}  
 \end{split}
```
- In the limit ``k>2k_F``

Define ``\vec{k}+\vec{p}=-\vec{p}' ``, then rename it as ``\vec{p}``,

```math
\begin{split}
=&  \int \frac{d^3 \vec{p}}{\left(2\pi^3\right)} \frac{f\left(\frac{p^2}{2 m}-\mu\right)-f\left(-\frac{p^2}{2 m}+\mu\right)}{i \Omega_n-\frac{p^2}{2 m}-\frac{(\vec{p}+\vec{k})^2}{2 m}+2 \mu} \\
=&  \int \frac{d^3 \vec{p}}{\left(2\pi^3\right)} \frac{2f\left(\frac{p^2}{2 m}-\mu\right)-1}{i \Omega_n-\frac{p^2}{2 m}-\frac{(\vec{p}+\vec{k})^2}{2 m}+2 \mu} \\
=& -\int \frac{d^3 \vec{p}}{\left(2\pi^3\right)} \frac{1}{i \Omega_n-\frac{p^2}{2 m}-\frac{(\vec{p}+\vec{k})^2}{2 m}+2 \mu} +\int_{|\vec{p}|< k_F} \frac{d^3 \vec{p}}{(2 \pi)^3} \frac{2}{i \Omega_n-\frac{p^2}{2 m}-\frac{\left(\vec{p}+k^2\right)}{2 m}+2 \mu}
\end{split}
```

Define ``\vec{p}'=\vec{p}+\vec{k}/2``, the first term becomes
```math
\frac{1}{i \Omega_n-\frac{\left(\vec{p}’+\vec{k}/2\right)^2}{2 m}-\frac{\left(\vec{p}^{\prime}-\vec{k}/{2}\right)^2}{2 m}+2\mu}=\frac{1}{i \Omega_n-\frac{\vec{p}^{\prime 2}}{m}-\frac{k^2}{4 m}+2 \mu}
```
```math
\begin{split}
=& -\int \frac{d^3 \vec{p}}{\left(2\pi^3\right)} \frac{1}{i \Omega_n-\frac{p^2}{ m}-\frac{\vec{k}^2}{4 m}+2 \mu}\\
=&4 \pi \int \frac{d p}{(2 \pi)^3} \frac{p^2}{\frac{p^2}{m}+\frac{k^2}{4 m}-i \Omega_n-2 \mu} \\
=&4 \pi m^{3/2} \int \frac{d p/m^{1/2}}{(2 \pi)^3} \frac{p^2/m}{\frac{p^2}{m}+\frac{k^2}{4 m}-i \Omega_n-2 \mu} \\
=&4 \pi m^{3/2} \int \frac{d x}{(2 \pi)^3} \frac{x^2}{x^2+\frac{k^2}{4 m}-i \Omega_n-2 \mu}
\end{split}
```
Using ``\int \frac{x^2 d x}{x^2+a}=x-\sqrt{a} \tan ^{-1}\left(\frac{x}{\sqrt{a}}\right)``,
```math
\begin{split}
=&4 \pi m^{3/2} \int \frac{d x}{(2 \pi)^3} \frac{x^2}{x^2+\frac{k^2}{4 m}-i \Omega_n-2 \mu} \\
=& \frac{4\pi m^{3/2}}{(2\pi)^3}\frac{\Lambda}{\sqrt{m}}-\frac{4\pi m^{3/2}}{(2\pi)^3}\sqrt{\frac{k^2}{4 m}-i \Omega_n-2 \mu}\cdot \left.\tan^{-1}\left(\frac{x}{\sqrt{\frac{k^2}{4 m}-i \Omega_n-2 \mu}}\right)\right|^{\Lambda/\sqrt{m}}_0 \\
=& \frac{m}{2\pi^2}\Lambda-\frac{m^{3/2}}{4\pi}\sqrt{\frac{k^2}{4 m}-i \Omega_n-2 \mu} 
\end{split}
```
We conclude
```math
\begin{split}
 &\int \frac{d^3 \vec{p}}{\left(2\pi^3\right)} T \sum_{i \omega_n} \frac{1}{i \omega_n+i \Omega_n-\frac{(\vec{k}+\vec{p})^2}{2 m}+\mu} \frac{1}{-i \omega_n-\frac{p^2}{2 m}+\mu} \\
=&\frac{m \Lambda}{2\pi^2} -\frac{m^{3 / 2}}{4\pi} \sqrt{-i \Omega_n+\frac{k^2}{4 m}-2 \mu} 
+\int_{|\vec{p}|< k_F} \frac{d^3 \vec{p}}{(2 \pi)^3} \frac{2}{i \Omega_n-\frac{p^2}{2 m}-\frac{\left(\vec{p}+\vec{k}\right)^2}{2 m}+2 \mu}
\end{split}
```

- Generic ``k``
```math
\begin{split}
=&  \int \frac{d^3 \vec{p}}{\left(2\pi^3\right)} \frac{2f\left(\frac{p^2}{2 m}-\mu\right)-1}{i \Omega_n-\frac{p^2}{2 m}-\frac{(\vec{p}+\vec{k})^2}{2 m}+2 \mu} \\
=& \frac{1}{(2\pi)^2} \int_0^\pi \sin(\theta)d\theta \int_0^\infty p^2 dp \frac{2f\left(\frac{p^2}{2 m}-\mu\right)-1}{i \Omega_n-\frac{p^2}{m}-\frac{k^2}{2m}-\frac{kp\cos(\theta)}{m}+2 \mu}
\end{split}
```
The angle can be integrated explicitly,
```math
-\int_0^\pi \frac{d\cos\theta}{a+b\cos\theta} = \int_{-1}^{1} \frac{dx}{a+bx} = \frac{1}{b} \ln\frac{a+b}{a-b}
```
where ``a`` is not real-valued as long as ``\Omega_n`` doesn't vanish, and ``\ln(a/b) \equiv \ln(a)-\ln(b)``.

Therefore, for generic Matsubara-frequency,
```math
= \frac{m}{(2\pi)^2} \int_0^\infty dp \left[\frac{p}{k}\ln \frac{i\Omega_n -\frac{p^2}{m}-\frac{k^2}{2m}+\frac{pk}{m}+2\mu}{i\Omega_n -\frac{p^2}{m}-\frac{k^2}{2m}-\frac{pk}{m}+2\mu}\left(2f(\frac{p^2}{2m}-\mu)-1\right)-2\right]+\frac{m\Lambda}{2\pi^2}
```

In the small momentum limit, it is simplifies to
```math
= \frac{m}{(2\pi)^2} \int_0^\infty dp \left[\frac{\frac{2p^2}{m}}{i\Omega_n -\frac{p^2}{m}+2\mu}\left(2f(\frac{p^2}{2m}-\mu)-1\right)-2\right]+\frac{m\Lambda}{2\pi^2}
```