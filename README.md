# Computer-assisted proofs of periodic solutions of state-dependent delay equations beyond perturbative regimes



Table of contents:


* [Introduction](#introduction)
* [Equations under study](#equations-under-study)
   * [A periodically forced toy model](#a-periodically-forced-toy-model)
   * [An autonomous SDDE with two state-dependent delays](#an-autonomous-sdde-with-two-state-dependent-delays)
   * [A gene regulatory model with threshold delay](#a-gene-regulatory-model-with-threshold-delay)
* [Utilisation and References](#utilisation-and-references)
* [License and Citation](#license-and-citation)
* [Contact](#contact)



# Introduction

This Matlab code is a complement to the article

#### [[1]](https://arxiv.org/abs/2610.36216) : "Computer-assisted proofs of periodic solutions of state-dependent delay equations beyond perturbative regimes", J. B. van den Berg, M. Breden, M. Cadiot, K. E. M. Church and J.-P. Lessard, [ArXiv Link](https://arxiv.org/abs/2610.36216).

in which we develop a computer-assisted framework for proving the existence of periodic solutions of state-dependent delay differential equations (SDDEs) directly, outside perturbative regimes. The approximate solution is a trigonometric polynomial, and a true periodic solution in its vicinity is obtained by a Newton-Kantorovich argument (Theorem 3) in the space of $C^1$ periodic functions. The Fourier coefficients of the state-dependent terms, such as $\bar{u}(t-\bar{\tau}(t))$, are rigorously enclosed by combining a geometric bound on a complex strip with trigonometric interpolation and a control of the aliasing error (Section 2.5).

This code provides the necessary rigorous computations of the bounds presented along the paper. The computations are performed using [Intlab](https://www.tuhh.de/ti3/intlab/) for interval arithmetic, and the Fourier coefficients are enclosed with interval-valued FFTs. The second and third applications are taken from the following papers :

#### [[2]](https://arxiv.org/abs/1607.02683) Resonance phenomena in a scalar delay differential equation with two state-dependent delays, R. C. Calleja, A. R. Humphries and B. Krauskopf. SIAM Journal on Applied Dynamical Systems, 16(3):1474–1513, 2017.
#### [[3]](https://arxiv.org/abs/2410.13092) Dynamics of a state-dependent delay-differential equation, T. Gedeon, A. R. Humphries, M. C. Mackey, H.-O. Walther and Z. Wang. Electronic Journal of Qualitative Theory of Differential Equations, Paper No. 37, 2025.


# Equations under study

We present some numerical details for the three applications of Sections 2, 3 and 4 of [[1]](https://arxiv.org/abs/2610.36216). Each application is given in a separate folder, in which we prove constructively the existence of a periodic solution in a neighborhood of a given approximate solution. In each case, the code computes rigorously the bounds $Y$, $Z_1$ and $Z_2$ of Theorem 3 and validates (or not) the computer-assisted proof. If the proof succeeds, the smallest radius $r_0$ for which Theorem 3 applies is displayed, and the code verifies that Theorem 3 holds with the radius $r_0$ stated in the paper.

## A periodically forced toy model

In Section 2 of [[1]](https://arxiv.org/abs/2610.36216), we introduce the main components of our approach on the toy model

$$u'(t) = -u(t) - \varepsilon_1 u(t-1-\varepsilon_2 u(t)) + \cos(t) + \sin(t),$$

for which we look for $2\pi$-periodic solutions. For $\varepsilon_1 = 0$, $u(t) = \sin(t)$ is an isolated $2\pi$-periodic solution, from which the parameters can be continued away from zero. The code "proof_non_autonomous.m" in the folder "1. Non-autonomous" first computes an approximate solution $\bar{u}$ for $\varepsilon_1 = 0.5$ and $\varepsilon_2 = 2$, given as a trigonometric polynomial of order 150, using a Newton method starting from $\sin(t)$. Then, choosing the truncation size $N = 660$, $N_{\text{FFT}} = 2^{12}$ and $\bar{\nu} = 1.1$, it computes the bounds $Y$, $Z_1$ and $Z_2$ of Lemmas 4, 9 and 5, using Lemmas 10, 11 and 13 and Proposition 8 for the constant $\theta_N$. We obtain

- $Y \leq 1.21 \times 10^{-7}$
- $Z_1 \leq 0.410$
- $Z_2 \leq 158.9$

and Theorem 3 provides a $2\pi$-periodic solution $\tilde{u}$ such that $\|\tilde{u}-\bar{u}\|_{C^1} \leq r_0 = 2.06 \times 10^{-7}$ (Theorem 14).

## An autonomous SDDE with two state-dependent delays

In Section 3 of [[1]](https://arxiv.org/abs/2610.36216), we consider the scalar SDDE with two state-dependent delays

$$u'(t) = -\gamma u(t) - \kappa_1 u(t-a_1-c_1u(t)) - \kappa_2 u(t-a_2-c_2u(t))$$

studied in [[2]](https://arxiv.org/abs/1607.02683), with $\gamma = 4.75$, $a_1 = 1.3$, $a_2 = 6$, $c_1 = c_2 = 1$ and $\kappa_2 = 3$. As $\kappa_1$ increases, the zero solution loses stability in a Hopf bifurcation at $\kappa_1 \approx 3.206$, and we prove the existence of a periodic orbit on the bifurcating branch for $\kappa_1 = 3.4$. Since the period is unknown, we use the time and amplitude rescalings (71)-(72), so that the unknowns are $\beta_1 = \omega$, $\beta_2 = u(-\pi)$ and a $2\pi$-periodic function $u$, and we add the phase condition (76).

The approximate solution $\bar{x} = (\bar{\beta}_1,\bar{\beta}_2,\bar{u})$ of Figure 2 in [[1]](https://arxiv.org/abs/2610.36216) is given in the file "sol1.mat", where $\bar{u}$ is a trigonometric polynomial of order 30 and $2\pi\bar{\beta}_1 \approx 2.578$ is the approximate period. The code "proof_autonomous.m" in the folder "2. Autonomous" computes the bounds of Lemmas 4, 9 (using Lemma 16 and (83)) and 17, with $N = 800$, $N_{\text{FFT}} = 2^{12}$, $\bar{\nu} = 1.1$ and $r_* = 10^{-3}$. We obtain

- $Y \leq 1.6 \times 10^{-11}$
- $Z_1 \leq 0.834$
- $Z_2(r_*) \leq 3.3 \times 10^{4}$

and Theorem 3 provides $\tilde{\beta}_1 > 0$ and $\tilde{u}$ such that $|\tilde{\beta}_1-\bar{\beta}_1| \leq r_0$ and $\|\tilde{u}-\bar{u}\|_{C^1} \leq r_0$, where $r_0 = 9.7 \times 10^{-11}$. In particular, $t \mapsto \tilde{\beta}_1\tilde{u}(t/\tilde{\beta}_1)$ is a $2\pi\tilde{\beta}_1$-periodic solution of the above equation (Theorem 18).

## A gene regulatory model with threshold delay

In Section 4 of [[1]](https://arxiv.org/abs/2610.36216), we consider the threshold-type SDDE studied in [[3]](https://arxiv.org/abs/2410.13092)

$$w'(t) = b e^{-\mu\tau(t)}\frac{v(w(t))}{v(w(t-\tau(t)))}q(w(t-\tau(t))) - \gamma w(t), \qquad a = \int_{-\tau(t)}^0 v(w(s+t))ds,$$

where $q$ and $v$ are the Hill functions

$$q(x) = \frac{q^{-}\theta_q^n + q^{+}x^n}{\theta_q^n + x^n} \quad \text{and} \quad v(x) = \frac{v^{-}\theta_v^m + v^{+}x^m}{\theta_v^m + x^m}.$$

Here, $w$ is the effector concentration, $q$ the production flux of messenger RNA and $v$ the transcription/translation velocity, so that the delay $\tau(t)$ is the time needed to process a strand of fixed length $a$. The model is reformulated as a system of SDDEs for $(w,\tau)$, with an integral constraint and an unfolding parameter (see (102)), the period being an unknown as in Section 3.

The approximate solution $\bar{x}$ of Figure 3 in [[1]](https://arxiv.org/abs/2610.36216), for $a = 3$, $b = 3$, $\mu = 0.3$, $\gamma = 0.13$, $m = 56$, $v^{-} = 0.5$, $v^{+} = 0.2$, $\theta_v = 1$ and $q^{-} = q^{+} = 1$ (the parameters of Figure 5.7 in [[3]](https://arxiv.org/abs/2410.13092)), is given in the file "pt2.mat". The function $\bar{w}$ and the rescaled delay $\bar{\tau}$ are trigonometric polynomials of order 70, and the approximate period is $2\pi\bar{\beta}_2 \approx 11.22$. Since $q$ and $v$ are rational functions, the code "proof_threshold.m" in the folder "3. Threshold" first verifies, using Lemma 23, that the terms of the vector field are analytic on the strip $|\text{Im}(z)| \leq \bar{\rho}$, where $\bar{\rho} = \ln(1.05)$. That is, the images of the strip by $\bar{w}$ and $\bar{w}(\cdot-\bar{\tau})$ do not contain the roots (117) of the denominators of the Hill functions, which is verified on the boundary of the strip thanks to the half-plane condition (118). Then, it computes the bounds of Lemmas 4, 26 and 27, with $N = 1000$, $N_2 = 2N$, $N_{\text{FFT}} = 2^{14}$ and $r_* = 10^{-3}$. We obtain

- $Y \leq 8.54 \times 10^{-10}$
- $Z_1 \leq 0.569$
- $Z_2(r_*) \leq 3.36 \times 10^{5}$

and Theorem 3 provides $\tilde{\beta}_2 > 0$ and $\tilde{u} = (\tilde{u}_1,\tilde{u}_2)$ such that $|\tilde{\beta}_2-\bar{\beta}_2| \leq r_0$ and $\|\tilde{u}-\bar{u}\|_{C^1} \leq r_0$, where $r_0 = 2 \times 10^{-9}$. In particular, $(w,\tau) : s \mapsto (\tilde{u}_1(s/\tilde{\beta}_2), \tilde{\beta}_2\tilde{u}_2(s/\tilde{\beta}_2))$ is a $2\pi\tilde{\beta}_2$-periodic solution of the threshold model (Theorem 28).


# Utilisation and References

Each folder is self-contained, and the main codes "proof_non_autonomous.m", "proof_autonomous.m" and "proof_threshold.m" can be run directly from their folder, once Intlab has been started (startintlab). The approximate solution of the toy model is computed on the run, while the approximate solutions of the two other applications are stored in the files "sol1.mat" and "pt2.mat", which need to be in the same folder as the main code. The folders contain the following functions :

- "computation_bounds.m" computes the bounds $Y$, $Z_1$ and $Z_2$, and is called by the main code.
- "F_SDDE_rig.m" and "DF_SDDE_rig.m" (resp. "F_threshold_rig.m" and "DF_threshold_rig.m") provide rigorous enclosures of the zero-finding map $F$ and of the Galerkin projection of its Jacobian at the approximate solution, in Fourier coefficients.
- "C_rho_contour_SDDE.m" and "C_rho_contour_operator.m" compute, on the complex strip, the constants of the geometric bound (34) on the Fourier coefficients, which are needed for the tails and for the aliasing error (38).
- "comp_mat_norm.m" and "DFT_matrix.m" compute operator norms of finite-dimensional operators on $H_0 = \mathbb{R}^m \times C^0$, following Section 2.6.6 and Lemma 13.
- "kolmogorov_theta_N.m" provides a rigorous enclosure of the constant $\theta_N$ of Lemma 6, which is the optimal constant in $\|\pi^{>N}h\|_{C^0} \leq \theta_N\|h'\|_{C^0}$, using Proposition 8.
- "DF_SDDE.m", "DF_autonomous.m" and "DF_threshold.m" compute the Jacobian in floating point arithmetic, in order to build the approximate inverse $A^{\leq N}$. Moreover, in the folder "1. Non-autonomous", "F_SDDE.m" and "newton.m" provide the Newton method for the toy model, and "plot_periodic.m" plots the approximate solution.
- In the folder "3. Threshold", "check_analyticity.m" and "avoidpole.m" verify the analyticity of the vector field on the strip (Lemma 23), and "Gderivs_threshold.m" encloses the second derivatives needed in Lemma 27.
- In the folders "1. Non-autonomous" and "2. Autonomous", "multiplicationmatrix.m", "tensortomatrix.m", "findfulltype.m" and "findnumberofmodes.m" build the matrices of multiplication operators in Fourier coefficients.

The comments in the codes refer to the equations, lemmas and sections of [[1]](https://arxiv.org/abs/2610.36216), and follow its notations.

The code is built using the following software :
 - [Matlab](https://www.mathworks.com/products/matlab.html)
 - [Intlab](https://www.tuhh.de/ti3/intlab/), for the interval arithmetic computations.


# License and Citation

This code is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).

If you wish to use this code in your publication, research, teaching, or other activities, please cite it using the following BibTeX template:

```
@software{SDDEsPeriodicOrbits.mat,
  author = {Jan Bouwe van den Berg and Maxime Breden and Matthieu Cadiot and Kevin Church and Jean-Philippe Lessard},
  title  = {{SDDEsPeriodicOrbits.mat}},
  url    = {https://github.com/matthieucadiot/SDDEsPeriodicOrbits.mat},
  note = {\url{ https://github.com/matthieucadiot/SDDEsPeriodicOrbits.mat}},
  year   = {2026}
}
```


# Contact

You can contact me at :

matthieu.cadiot@polytechnique.edu
