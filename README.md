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

in which we prove the existence of periodic solutions of state-dependent delay differential equations (SDDEs) outside perturbative regimes. A true periodic solution is obtained in a neighborhood of an approximate solution, given as a trigonometric polynomial, using a Newton-Kantorovich argument in the space of $`C^1`$ periodic functions.

This code provides the necessary rigorous computations of the bounds presented along the paper. The computations are performed using [Intlab](https://www.tuhh.de/ti3/intlab/) for interval arithmetic. The second and third applications are taken from the following papers :

#### [[2]](https://arxiv.org/abs/1607.02683) Resonance phenomena in a scalar delay differential equation with two state-dependent delays, R. C. Calleja, A. R. Humphries and B. Krauskopf. SIAM Journal on Applied Dynamical Systems, 16(3):1474–1513, 2017.
#### [[3]](https://arxiv.org/abs/2410.13092) Dynamics of a state-dependent delay-differential equation, T. Gedeon, A. R. Humphries, M. C. Mackey, H.-O. Walther and Z. Wang. Electronic Journal of Qualitative Theory of Differential Equations, Paper No. 37, 2025.


# Equations under study

Each application of [[1]](https://arxiv.org/abs/2610.36216) is given in a separate folder. In each case, the main code computes rigorously the bounds $`Y`$, $`Z_1`$ and $`Z_2`$ of Theorem 3 in [[1]](https://arxiv.org/abs/2610.36216) and validates (or not) the computer-assisted proof. If the proof succeeds, the radius $`r_0`$ of the ball containing the true solution is displayed.

## A periodically forced toy model

In Section 2 of [[1]](https://arxiv.org/abs/2610.36216), we consider the toy model

```math
u'(t) = -u(t) - \varepsilon_1 u(t-1-\varepsilon_2 u(t)) + \cos(t) + \sin(t)
```

and look for $`2\pi`$-periodic solutions. The code "proof_non_autonomous.m" in the folder "1. Non-autonomous" first computes an approximate solution $`\bar{u}`$ for $`\varepsilon_1 = 0.5`$ and $`\varepsilon_2 = 2`$, using a Newton method starting from $`\sin(t)`$ (which is the solution for $`\varepsilon_1 = 0`$). Then, it computes the bounds of the proof and we obtain

- $`Y \leq 1.21 \times 10^{-7}`$
- $`Z_1 \leq 0.410`$
- $`Z_2 \leq 158.9`$

This proves the existence of a $`2\pi`$-periodic solution $`\tilde{u}`$ such that $`\|\tilde{u}-\bar{u}\|_{C^1} \leq 2.06 \times 10^{-7}`$ (Theorem 14).

## An autonomous SDDE with two state-dependent delays

In Section 3 of [[1]](https://arxiv.org/abs/2610.36216), we consider the scalar SDDE with two state-dependent delays

```math
u'(t) = -\gamma u(t) - \kappa_1 u(t-a_1-c_1u(t)) - \kappa_2 u(t-a_2-c_2u(t))
```

studied in [[2]](https://arxiv.org/abs/1607.02683), with $`\gamma = 4.75`$, $`a_1 = 1.3`$, $`a_2 = 6`$, $`c_1 = c_2 = 1`$, $`\kappa_1 = 3.4`$ and $`\kappa_2 = 3`$. Since the equation is autonomous, the period is an additional unknown. The approximate solution is stored in the file "sol1.mat", and the code "proof_autonomous.m" in the folder "2. Autonomous" computes the bounds of the proof. We obtain

- $`Y \leq 1.6 \times 10^{-11}`$
- $`Z_1 \leq 0.834`$
- $`Z_2 \leq 3.3 \times 10^{4}`$

This proves the existence of a periodic solution, with period close to $`2.578`$, at distance at most $`r_0 = 9.7 \times 10^{-11}`$ of the approximate solution (Theorem 18).

## A gene regulatory model with threshold delay

In Section 4 of [[1]](https://arxiv.org/abs/2610.36216), we consider the threshold-type SDDE studied in [[3]](https://arxiv.org/abs/2410.13092)

```math
w'(t) = b e^{-\mu\tau(t)}\frac{v(w(t))}{v(w(t-\tau(t)))}q(w(t-\tau(t))) - \gamma w(t), \qquad a = \int_{-\tau(t)}^0 v(w(s+t))\,ds,
```

where $`q`$ and $`v`$ are the Hill functions

```math
q(x) = \frac{q^{-}\theta_q^n + q^{+}x^n}{\theta_q^n + x^n}, \qquad v(x) = \frac{v^{-}\theta_v^m + v^{+}x^m}{\theta_v^m + x^m}.
```

The approximate solution, for $`a = 3`$, $`b = 3`$, $`\mu = 0.3`$, $`\gamma = 0.13`$, $`m = 56`$, $`v^{-} = 0.5`$, $`v^{+} = 0.2`$, $`\theta_v = 1`$ and $`q^{-} = q^{+} = 1`$, is stored in the file "pt2.mat". The code "proof_threshold.m" in the folder "3. Threshold" first verifies that the vector field is analytic on a complex strip, and then computes the bounds of the proof. We obtain

- $`Y \leq 8.54 \times 10^{-10}`$
- $`Z_1 \leq 0.569`$
- $`Z_2 \leq 3.36 \times 10^{5}`$

This proves the existence of a periodic solution, with period close to $`11.22`$, at distance at most $`r_0 = 2 \times 10^{-9}`$ of the approximate solution (Theorem 28).


# Utilisation and References

The codes require [Matlab](https://www.mathworks.com/products/matlab.html) and [Intlab](https://www.tuhh.de/ti3/intlab/). To run a proof, start Intlab in Matlab (startintlab), go to the folder of the example and run its main code :

- "proof_non_autonomous.m" in the folder "1. Non-autonomous" (Theorem 14),
- "proof_autonomous.m" in the folder "2. Autonomous" (Theorem 18),
- "proof_threshold.m" in the folder "3. Threshold" (Theorem 28).

Each folder is self-contained : it contains all the functions and data needed for the corresponding proof. The comments in the codes refer to the equations, lemmas and sections of [[1]](https://arxiv.org/abs/2610.36216).


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
