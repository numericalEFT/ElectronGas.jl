
# Quasiparticle properties of electon gas

## Renormalization factor
 The renormalization constant $Z$ gives the strength of the quasiparticle pole, and can be obtained from the frequency dependence of the self-energy as
 ```math
Z=\frac{1}{1-\left.\frac{1}{\hbar} \frac{\partial \operatorname{Im} \Sigma(k, i\omega_n)}{\partial \omega_n}\right|_{k=k_{F}, \omega_n=0^+}}
 ```

 ## Effective mass
 ```math
\frac{m^{*}}{m}= \frac{Z^{-1}}{1+\frac{m}{\hbar^{2} k_{F}} \left. \frac{\partial \operatorname{Re} \Sigma(k, i\omega_n)}{\partial k}\right|_{k=k_{F}, \omega_n=0^+}} 
 ```

 ## Benchmark 
### 2D UEG
| $r_s$ | $Z$ (RPA) | $Z$ ($G_0W_0$ [1]) | $m^*/m$ (RPA) | $m^*/m$ ($G_0W_0$ [1]) |
| :---: | :-------: | :----------------: | :-----------: | :--------------------: |
|  0.5  |   0.787   |       0.786        |               |         0.981          |
|  1.0  |   0.662   |       0.662        |               |         1.020          |
|  2.0  |   0.519   |       0.519        |               |         1.078          |
|  3.0  |   0.437   |       0.437        |               |         1.117          |
|  4.0  |   0.383   |       0.383        |               |         1.143          |
|  5.0  |   0.344   |       0.344        |               |         1.162          |
|  8.0  |   0.271   |       0.270        |               |         1.196          |
| 10.0  |   0.240   |       0.240        |               |         1.209          |

 ### 3D UEG
| $r_s$ | $Z$ (RPA) |       $Z$ ($G_0W_0$)        | $m^*/m$ (RPA) | $m^*/m$ ($G_0W_0$ [2]) |
| :---: | :-------: | :-------------------------: | :-----------: | :--------------------: |
|  1.0  |  0.8601   |        0.859 [**3**]        |               |         0.970          |
|  2.0  |  0.7642   | 0.768 [**3**] 0.764 [**4**] |               |         0.992          |
|  3.0  |  0.6927   |        0.700 [**3**]        |               |         1.016          |
|  4.0  |  0.6367   | 0.646 [**3**] 0.645 [**4**] |               |         1.039          |
|  5.0  |  0.5913   |        0.602 [**3**]        |               |         1.059          |
|  6.0  |  0.5535   |        0.568 [**3**]        |               |         1.078          |


 
[**References**]
1. [H.-J. Schulze, P. Schuck, and N. Van Giai, Two-dimensional electron gas in the random-phase approximation with exchange and self-energy corrections. *Phys. Rev. B 61, 8026* (2000).](https://link.aps.org/doi/10.1103/PhysRevB.61.8026)
2. [Simion, G. E. & Giuliani, G. F., Many-body local fields theory of quasiparticle properties in a three-dimensional electron liquid. *Phys. Rev. B 77, 035131* (2008).](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.77.035131)
3. G. D Mahan, *Many-Particle Physics* (Plenum, New York, 1991), Chap. 5. 
4. [B. Holm and U. von Barth, Fully self-consistent GW self-energy of the electron gas. *Phys. Rev. B 57, 2108* (1998).](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.57.2108)
