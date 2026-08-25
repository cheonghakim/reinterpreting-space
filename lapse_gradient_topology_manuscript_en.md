# Lapse-Gradient Topology of Static Multi-Black-Hole Spacetimes: Critical-Point Bifurcations in Majumdar-Papapetrou Geometry

## Abstract

This study treats the lapse function \(N=d\tau/dt\) of a static spacetime as a scalar field encoding gravitational redshift and the proper-acceleration structure of static observers, and proposes a gradient-topological framework for classifying the global organization of multi-source gravitational fields.

As an exact benchmark, we use Majumdar-Papapetrou (MP) multi-black-hole spacetimes, for which the lapse is determined by a harmonic function \(U\). We define a descriptor based on finite critical points of the lapse, their Morse indices, gradient separatrices, and heteroclinic structure.

For an equilateral three-center family with masses

$$
(M_1,M_2,M_3)=(1,1,\lambda),
$$

we identify two critical mass ratios,

$$
\lambda_-
=
0.673877474470738\ldots,
$$

and

$$
\lambda_+
=
1.136252210664681\ldots,
$$

at which the lapse critical-point structure bifurcates.

The first bifurcation is a symmetry-breaking pitchfork-type bifurcation and the second is a fold-type saddle-node annihilation. The finite critical structure changes as

$$
2S_1
\longrightarrow
3S_1+S_2
\longrightarrow
2S_1.
$$

These results are supported by analytic normal-form reduction, interval/Krawczyk root certification at representative parameter values, Brouwer-index consistency, and independent validated numerics using outward-rounded interval arithmetic and arbitrary-precision ball arithmetic.

In particular, complete critical sets and Morse indices are certified at \(\lambda=0.5\), \(0.8\), and \(1.2\). Local existence and uniqueness of the off-axis pitchfork branches and fold nondegeneracy are also certified.

For the asymmetric family

$$
(M_1,M_2,M_3)
=
(1+\varepsilon,1-\varepsilon,\lambda),
$$

we obtain a cusp unfolding with the scaling law

$$
\lambda_{\mathrm{fold}}(\varepsilon)-\lambda_-
\sim
C|\varepsilon|^{2/3},
\qquad
C=1.58535311676\ldots.
$$

Finally, we compare the proposed gradient-skeleton representation with contour/Reeb-type representations and argue that the distinguishing information is not level-set connectivity alone, but the combined organization of Morse indices, stable and unstable manifolds, and heteroclinic connections.

**Keywords:** lapse function, gravitational redshift, Majumdar-Papapetrou spacetime, Morse bifurcation, gradient skeleton, validated numerics, interval arithmetic, multi-black-hole topology

---

## 1. Introduction

In general relativity, the lapse function \(N\) of a static spacetime is not merely a coordinate bookkeeping device. For static observers, it gives the ratio between proper time and coordinate time,

$$
N=\frac{d\tau}{dt},
$$

and therefore directly represents the local gravitational redshift factor.

The original motivation of this study was to reinterpret ordinary three-dimensional space by lifting the position-dependent rate of proper time into an additional graph coordinate. In that sense one may visualize a point as

$$
(x,y,z,N),
$$

where the fourth coordinate is not interpreted as a new physical spatial dimension, but as an auxiliary representation of the redshift field.

The physically meaningful research question is therefore:

> Can the lapse/redshift field of a static gravitational system be classified globally through its finite critical points, Morse indices, gradient flows, stable and unstable manifolds, and heteroclinic connections?

This question goes beyond a simple visualization of equal-lapse surfaces. The four-acceleration of a static observer satisfies

$$
a_i=c^2D_i\ln N,
$$

so the spatial gradient of \(N\) is directly related to the proper acceleration required to remain static.

The present work studies this structure in exact Majumdar-Papapetrou multi-black-hole spacetimes. The MP family is particularly suitable because the lapse is

$$
N=U^{-1},
$$

with

$$
U(\mathbf{x})
=
1+\sum_a\frac{M_a}{|\mathbf{x}-\mathbf{x}_a|}.
$$

This makes the critical-point structure analytically tractable while retaining a fully relativistic exact solution.

The main contributions of this work are:

1. We define a **gravitational gradient skeleton** based on lapse critical points, Morse indices, separatrices, heteroclinic connectivity, and source/horizon basins.
2. We derive two exact bifurcation values for the equilateral MP three-center family \((1,1,\lambda)\).
3. We show that the finite lapse critical structure changes as

   $$
   2S_1
   \longrightarrow
   3S_1+S_2
   \longrightarrow
   2S_1.
   $$

4. We certify representative phase structures using interval arithmetic and Krawczyk-type validated root isolation.
5. We extend the analysis to the asymmetric family \((1+\varepsilon,1-\varepsilon,\lambda)\), obtaining a cusp unfolding and a validated \(2/3\)-power scaling.
6. We compare the gradient skeleton with contour/Reeb-type scalar-field representations.

---

## 2. Lapse Graphs and Gravitational Gradient Skeletons

Consider a static spacetime with metric

$$
ds^2
=
-N^2c^2dt^2
+
h_{ij}dx^idx^j.
$$

For a static observer,

$$
N=\frac{d\tau}{dt}.
$$

Thus \(N:S\rightarrow\mathbb{R}_+\) is naturally interpreted as a spatial redshift field on the static spatial manifold \(S\).

We may lift this scalar field into an auxiliary product space \(S\times\mathbb{R}\) by defining

$$
F_\ell(p)
=
(p,\ell N(p)),
$$

where \(\ell>0\) is an arbitrary visualization scale.

The resulting graph is

$$
\Sigma_N^\ell
=
\left\{
(p,\ell N(p))
\mid
p\in S
\right\}.
$$

The induced metric is

$$
\tilde h
=
h+\ell^2dN\otimes dN.
$$

This construction does **not** introduce a new physical dimension. It is only an auxiliary embedding used to represent the scalar field geometrically.

### Definition 2.1. Gravitational Gradient Skeleton

We define the gravitational gradient skeleton of the lapse field as

$$
\mathcal S_N
=
(C,I,V,\Gamma,\mathcal B),
$$

where

- \(C\): finite critical points of \(N\),
- \(I\): Morse index of each critical point,
- \(V\): critical lapse values \(N(c)\),
- \(\Gamma\): separatrix and heteroclinic connectivity,
- \(\mathcal B\): source/horizon basin structure.

This descriptor contains more information than basin adjacency alone.

---

## 3. Physical Meaning of the Lapse Gradient

For a static observer, the four-acceleration satisfies

$$
a_i
=
c^2D_i\ln N.
$$

Equivalently,

$$
a_i
=
\frac{c^2}{N}D_iN.
$$

Thus the spatial gradient of the lapse is directly related to the proper acceleration required to remain static.

For the Schwarzschild spacetime,

$$
N(r)
=
\sqrt{
1-\frac{2GM}{rc^2}
}.
$$

The spatial gradient norm satisfies

$$
|DN|_h
=
\frac{GM}{c^2r^2}.
$$

Hence the proper acceleration of a static observer is

$$
a
=
\frac{GM}
{r^2\sqrt{1-\frac{2GM}{rc^2}}}.
$$

In the weak-field limit,

$$
N
=
1+\frac{\Phi}{c^2}
+
O(c^{-4}),
$$

so the Newtonian gravitational acceleration is recovered as

$$
\mathbf g_{\rm N}
=
-c^2\nabla N
+
O(c^{-2}).
$$

---

## 4. Majumdar-Papapetrou Spacetimes

The Majumdar-Papapetrou metric can be written as

$$
ds^2
=
-U^{-2}dt^2
+
U^2d\mathbf{x}^2,
$$

where

$$
U(\mathbf{x})
=
1+\sum_{a=1}^{n}\frac{M_a}{r_a},
\qquad
r_a
=
|\mathbf{x}-\mathbf{x}_a|.
$$

The lapse is therefore

$$
N=U^{-1}.
$$

Since

$$
\nabla N
=
-U^{-2}\nabla U,
$$

finite critical points of \(N\) are exactly the finite critical points of \(U\):

$$
\nabla N=0
\iff
\nabla U=0.
$$

Furthermore, because the spatial metric is

$$
h=U^2\delta,
$$

the physical and Euclidean gradient fields differ only by a positive scalar factor. Their integral curves are therefore identical up to reparameterization.

---

## 5. Binary Majumdar-Papapetrou Case

Consider two masses \(M_1\) and \(M_2\) located at \(x=-1\) and \(x=1\).

The finite critical point between the two sources satisfies

$$
\frac{M_1}{(x_*+1)^2}
=
\frac{M_2}{(1-x_*)^2},
$$

which gives

$$
x_*
=
\frac{\sqrt{M_1}-\sqrt{M_2}}
{\sqrt{M_1}+\sqrt{M_2}}.
$$

The Hessian signature of the lapse at this point is

$$
(-,+,+),
$$

so the binary MP system contains exactly one finite Morse-index-1 lapse saddle,

$$
1\times S_1.
$$

---

## 6. Equilateral Three-Center Family

Place the three sources at

$$
P_1=(1,0),
$$

$$
P_2=
\left(
-\frac12,
\frac{\sqrt3}{2}
\right),
$$

and

$$
P_3=
\left(
-\frac12,
-\frac{\sqrt3}{2}
\right).
$$

We first study the one-parameter family

$$
(M_1,M_2,M_3)
=
(1,1,\lambda).
$$

Here the varying mass \(\lambda\) is attached to \(P_1=(1,0)\), while the two
equal unit masses sit at \(P_2,P_3\); the reflection \(P_2\leftrightarrow P_3\)
is therefore the exact \(\mathbb Z_2\) symmetry of the configuration, with
fixed point \(P_1\).

For the fully symmetric case \(\lambda=1\), there are four finite lapse critical points.

The critical structure is

$$
3S_1+S_2.
$$

The central critical point has

$$
N=\frac14
$$

and Morse index 2, while the other three critical points are symmetry-related Morse-index-1 saddles.

---

## 7. Exact Bifurcation Analysis

### 7.1 Symmetry-Adapted Coordinates

Let \(r\) be the coordinate along the symmetry axis through \(P_1\), and let \(t\) be the transverse coordinate in the source plane.

The squared distances are

$$
R_1^2
=
r^2+r+t^2-\sqrt3\,t+1,
$$

$$
R_2^2
=
r^2+r+t^2+\sqrt3\,t+1,
$$

and

$$
R_3^2
=
(1-r)^2+t^2.
$$

Therefore,

$$
U(r,t;\lambda)
=
1+
\frac1{\sqrt{R_1^2}}
+
\frac1{\sqrt{R_2^2}}
+
\frac{\lambda}{\sqrt{R_3^2}}.
$$

Because

$$
U(r,-t;\lambda)
=
U(r,t;\lambda),
$$

the system has an exact \(\mathbb Z_2\) reflection symmetry.

---

### 7.2 First Bifurcation: Symmetry-Breaking Pitchfork

On the symmetry axis, the critical branch satisfies

$$
\lambda(r)
=
\frac{(2r+1)(1-r)^2}
{(r^2+r+1)^{3/2}}.
$$

The first degeneracy occurs at

$$
r_-
=
\frac{-5+\sqrt{33}}{4},
$$

corresponding to

$$
\lambda_-
=
0.67387747447073801137622337017\ldots.
$$

A Lyapunov-Schmidt reduction gives

$$
g(t,\mu)
=
t
\left[
a\mu
+
\beta t^2
+
\cdots
\right],
$$

where

$$
\mu
=
\lambda-\lambda_-,
$$

and

$$
a
=
3.55212027148\ldots
>
0,
$$

$$
\beta
=
-16.0449298058\ldots
<
0.
$$

Thus the off-axis branches satisfy

$$
t_\pm
=
\pm
0.4705165678
\sqrt{\lambda-\lambda_-}
+
O\!\left(
(\lambda-\lambda_-)^{3/2}
\right).
$$

The local critical-point structure changes as

$$
S_1
\longrightarrow
S_2+2S_1.
$$

Including the second pre-existing index-1 saddle, the full structure changes as

$$
2S_1
\longrightarrow
3S_1+S_2.
$$

---

### 7.3 Second Bifurcation: Fold

The second critical parameter is obtained at

$$
r_+
=
\frac{-7+\sqrt{33}}{8},
$$

with

$$
\lambda_+
=
1.13625221066468090151835125249\ldots.
$$

At this point,

$$
U_r=0,
\qquad
U_{rr}=0,
$$

while

$$
U_{r\lambda}\neq0,
$$

and

$$
U_{rrr}\neq0.
$$

Hence the second bifurcation is a generic fold.

An index-1 and an index-2 saddle annihilate,

$$
S_1+S_2
\longrightarrow
\varnothing.
$$

The complete finite critical structure therefore follows the sequence

$$
\boxed{
2S_1
\longrightarrow
3S_1+S_2
\longrightarrow
2S_1
}.
$$

---

## 8. Symmetry Breaking and Cusp Unfolding

We now perturb the masses as

$$
(M_1,M_2,M_3)
=
(1+\varepsilon,1-\varepsilon,\lambda).
$$

The reduced normal form becomes

$$
g(t,\mu,\varepsilon)
=
k\varepsilon
+
a\mu t
+
\beta t^3
+
\cdots,
$$

with

$$
k
=
1.28410256136\ldots.
$$

This is an imperfect pitchfork, or cusp unfolding.

The fold boundary obeys

$$
\lambda_{\rm fold}(\varepsilon)
-
\lambda_-
\sim
C|\varepsilon|^{2/3},
$$

where

$$
C
=
1.58535311676\ldots.
$$

Validated calculations gave

$$
R(10^{-5})=1.58564,
$$

$$
R(10^{-4})=1.58666,
$$

$$
R(10^{-3})=1.59132,
$$

and

$$
R(10^{-2})=1.61074,
$$

for

$$
R(\varepsilon)
=
\frac{
\lambda_{\rm fold}(\varepsilon)-\lambda_-
}{
|\varepsilon|^{2/3}
}.
$$

The values converge toward \(C\) as \(\varepsilon\to0\).

The two-dimensional parameter-space continuation further gives the terminal symmetry-restoration value

$$
\varepsilon_c
=
\frac{1-\lambda_-}{1+\lambda_-}
=
0.194830583781157\ldots,
$$

with

$$
\lambda_c
=
1+\varepsilon_c
=
1.194830583781157\ldots.
$$

---

## 9. Gradient-Skeleton Visualization

The following figure compares the planar lapse-gradient skeleton in three representative phases:

- \(\lambda=0.5\): \(2S_1\),
- \(\lambda=0.8\): \(3S_1+S_2\),
- \(\lambda=1.2\): \(2S_1\).

![Planar lapse-gradient skeleton across representative MP phases](figures/mp_gradient_skeleton_phase_comparison.png)

**Figure 1.** Planar lapse-gradient skeleton across representative MP phases. Solid curves indicate unstable branches of the descending lapse flow, while dashed curves indicate stable manifolds traced backward. The middle phase uniquely contains an additional index-2 saddle.

For the middle phase \(\lambda=0.8\), the three-dimensional stable-manifold separator structure was approximated using trajectory fans emitted from the positive eigenspaces of the index-1 lapse saddles.

![Three-dimensional lapse-gradient skeleton at lambda = 0.8](figures/mp_lambda_08_stable_manifolds.png)

**Figure 2.** Three-dimensional lapse-gradient skeleton and separator-sheet approximation for \(\lambda=0.8\). Stable manifolds of index-1 saddles are visualized as trajectory fans; unstable branches are shown as thicker curves.

The asymmetric phase diagram is shown below.

![Critical-structure phase diagram](figures/mp_three_center_phase_diagram.png)

**Figure 3.** Critical-structure phase diagram for the asymmetric family \((1+\varepsilon,1-\varepsilon,\lambda)\). The interior region corresponds to the finite critical structure \(3S_1+S_2\), while the exterior corresponds to \(2S_1\).

---

## 10. Validated Numerics

### 10.1 Independently Reproduced Results

The following results were independently reproduced using outward-rounded interval arithmetic and arbitrary-precision ball arithmetic:

- Arb enclosures for the principal constants,
- local pitchfork inequalities \(V_{rr}>0\) and \(Q_s<0\),
- off-axis pitchfork branch existence,
- off-axis pitchfork branch uniqueness,
- fold nondegeneracy,
- complete critical-point sets at \(\lambda=0.5\),
- complete critical-point sets at \(\lambda=0.8\),
- complete critical-point sets at \(\lambda=1.2\),
- Brouwer-index and boundary-winding consistency,
- asymmetric cusp scaling,
- the on-axis middle-slab mechanism used in the phase analysis.

The off-axis global exclusion and asymptotic tail slabs remain open.

---

### 10.2 Representative Completeness Certificates

At the three representative values,

$$
\lambda=0.5,
\qquad
\lambda=0.8,
\qquad
\lambda=1.2,
$$

the complete finite critical structures were certified as

$$
\lambda=0.5:
\qquad
2S_1,
$$

$$
\lambda=0.8:
\qquad
3S_1+S_2,
$$

and

$$
\lambda=1.2:
\qquad
2S_1.
$$

The certification procedure used:

1. analytic reduction to the source plane,
2. localization to the convex hull,
3. source-neighborhood dominance estimates,
4. interval gradient exclusion,
5. Krawczyk existence and uniqueness tests,
6. interval Hessian-sign classification.

No unresolved boxes remained at the representative parameter values.

---

### 10.3 Brouwer Index and Boundary Winding

The sum of the Brouwer indices of the finite regular critical points is

$$
-2.
$$

Each of the three \(1/r\)-type source singularities contributes \(+1\) to the winding of an outer boundary enclosing all sources.

Hence

$$
+1
=
(-2)+3.
$$

This resolves the distinction between the index sum of regular finite critical points and the total winding measured on a boundary enclosing the singular sources.

---

### 10.4 Certified Local Pitchfork Monotonicity

To avoid the removable singularity associated with \(U_t/t\) at \(t=0\), define

$$
s=t^2.
$$

Let

$$
a=r^2+r+1+s,
$$

and

$$
D=a^2-3s.
$$

The two symmetry-related source terms can be combined exactly as

$$
\frac1{\sqrt{a-\sqrt{3s}}}
+
\frac1{\sqrt{a+\sqrt{3s}}}
=
\sqrt{
\frac{2a}{D}
+
\frac{2}{\sqrt D}
}.
$$

Thus define

$$
V(r,s,\lambda)
=
1+
\sqrt{
\frac{2a}{D}
+
\frac{2}{\sqrt D}
}
+
\frac{\lambda}
{\sqrt{(1-r)^2+s}}.
$$

On the implicit radial branch \(V_r=0\),

$$
Q
=
\frac{U_t}{t}
=
2V_s,
$$

and

$$
Q_s
=
2
\left(
V_{ss}
-
\frac{V_{sr}^2}{V_{rr}}
\right).
$$

On the explicit interval box

$$
r
\in
[r_- -0.0007,\,
 r_- +0.0014],
$$

$$
s\in[0,0.0004],
$$

and

$$
\lambda
\in
[\lambda_- -0.001,\,
 \lambda_- +0.001],
$$

the independent interval calculation certified

$$
V_{rr}>0
$$

and

$$
Q_s<0.
$$

This gives uniqueness of the positive off-axis branch in \(s\), and hence exactly two symmetry-related off-axis branches in \(t\).

---

## 11. Main Certified Results

### Theorem 11.1. Certified Local Pitchfork Bifurcation

For the equilateral three-center MP family

$$
(M_1,M_2,M_3)
=
(1,1,\lambda),
$$

a symmetry-breaking pitchfork bifurcation occurs at

$$
\lambda
=
\lambda_-
=
0.673877474470738\ldots.
$$

Within an explicit validated local box,

$$
V_{rr}>0
$$

and

$$
Q_s<0.
$$

Therefore:

- for \(\lambda<\lambda_-\), there is no off-axis finite critical point in the certified neighborhood;
- for \(\lambda>\lambda_-\), exactly one positive \(s=t^2\) solution exists, producing exactly two symmetry-related off-axis critical points.

The local Morse structure changes as

$$
S_1
\longrightarrow
S_2+2S_1.
$$

---

### Theorem 11.2. Certified Local Fold Bifurcation

A generic fold bifurcation occurs at

$$
\lambda
=
\lambda_+
=
1.136252210664681\ldots.
$$

At the fold,

$$
U_r=0,
$$

$$
U_{rr}=0,
$$

while

$$
U_{r\lambda}\neq0,
$$

and

$$
U_{rrr}\neq0.
$$

Thus locally,

$$
S_1+S_2
\longrightarrow
\varnothing.
$$

---

### Proposition 11.3. Certified Representative Phases

The complete finite critical structures at the representative parameters are

$$
\lambda=0.5:
\qquad
2S_1,
$$

$$
\lambda=0.8:
\qquad
3S_1+S_2,
$$

and

$$
\lambda=1.2:
\qquad
2S_1.
$$

---

### Corollary 11.4. Partially Certified Phase Pattern

The analytic bifurcation structure, certified local bifurcations, complete representative critical sets, and validated on-axis middle-slab analysis jointly support the phase pattern

$$
2S_1
\longrightarrow
3S_1+S_2
\longrightarrow
2S_1.
$$

However, a complete global theorem for all

$$
0<\lambda<\infty
$$

is not yet claimed because off-axis degeneracy exclusion and the asymptotic tail slabs have not yet been fully interval-certified.

---

## 12. Comparison with Contour and Reeb Representations

Contour trees and Reeb graphs provide powerful representations of scalar-field level-set connectivity.

They can identify topological changes of level sets at critical values and therefore can, in principle, detect many of the same scalar-field bifurcation events.

The present work does **not** claim that lapse bifurcations are invisible to contour/Reeb methods.

The distinction lies in the information represented.

### Contour/Reeb representation

Primarily encodes

$$
\text{level-set component connectivity}.
$$

### Coarse basin adjacency

Primarily encodes

$$
\text{which source/horizon basins touch}.
$$

### Gravitational gradient skeleton

Encodes

$$
\text{critical positions}
+
\text{Morse indices}
+
\text{stable manifolds}
+
\text{unstable manifolds}
+
\text{heteroclinic organization}
+
\text{basin structure}.
$$

For the representative MP configurations, the coarse source-basin adjacency remains the complete graph \(K_3\) across all three tested phases.

Thus basin adjacency alone does not distinguish the phases.

In contrast, the gradient skeleton immediately distinguishes

$$
(n_1,n_2)
=
(2,0)
$$

from

$$
(n_1,n_2)
=
(3,1).
$$

This provides a compact discrete descriptor of the internal organization of the redshift field even when the coarse basin adjacency remains unchanged.

The proposed method should therefore be interpreted as complementary to contour/Reeb topology rather than as a replacement for it.

A future fully quantitative comparison should compare, for the same parameter sweep:

1. number of contour/Reeb critical events,
2. number and indices of lapse critical points,
3. heteroclinic connectivity changes,
4. representation size,
5. stability under perturbation,
6. computational cost.

---

## 13. Discussion

The original intuition behind this study was to reinterpret gravitationally distorted space through the spatial variation of proper-time flow.

The analysis shows that this intuition becomes mathematically useful when reformulated not as a new physical spatial dimension, but as the global topology of the lapse/redshift field.

The principal results are:

1. The lapse field has a physically meaningful gradient because

   $$
   a_i=c^2D_i\ln N.
   $$

2. Exact MP geometries provide a tractable relativistic setting in which lapse critical topology can be studied analytically.
3. The three-center system exhibits nontrivial Morse bifurcations.
4. Symmetry breaking produces a cusp unfolding with a validated \(2/3\)-power scaling.
5. Coarse basin adjacency is insufficient to distinguish all phases.
6. The combined Morse-index and gradient-manifold structure provides a more informative descriptor.
7. Validated numerics proved useful not only for confirmation but also for detecting and correcting conceptual and numerical errors.

During the validation process, several genuine issues were identified and corrected:

- boundary-root pathology in symmetric bisection,
- severe interval dependency in direct function evaluation,
- the sign of \(U_{zz}\),
- precision-floor effects under deep subdivision,
- the distinction between regular-critical-point index sum and total boundary winding in the presence of singular sources.

These corrections strengthen the credibility of the final formulation.

---

## 14. Limitations and Open Problems

The main unresolved issue is the global exclusion of additional off-axis degeneracies over the entire parameter range.

The currently certified scope includes:

- local pitchfork certification,
- local fold certification,
- complete representative critical sets,
- on-axis middle-slab degeneracy analysis,
- asymmetric cusp scaling.

The following remain open:

1. global off-axis degeneracy exclusion,
2. fully validated small-\(\lambda\) tail,
3. fully validated large-\(\lambda\) tail,
4. fully certified global phase theorem,
5. quantitative contour/Reeb versus gradient-skeleton benchmarking.

A natural mathematical question is therefore:

> Can all off-axis degeneracies be globally excluded for the full \((1,1,\lambda)\) MP family?

---

## 15. Conclusion

We introduced a lapse-gradient topological framework for static multi-source spacetimes and applied it to Majumdar-Papapetrou multi-black-hole geometries.

For the equilateral three-center family, we derived two exact bifurcation values,

$$
\lambda_-
=
0.673877474470738\ldots,
$$

and

$$
\lambda_+
=
1.136252210664681\ldots,
$$

and showed that the finite lapse critical structure follows the pattern

$$
2S_1
\longrightarrow
3S_1+S_2
\longrightarrow
2S_1.
$$

The local pitchfork and fold mechanisms, complete representative critical sets, and asymmetric cusp scaling were independently reproduced using validated numerics.

The main conceptual contribution is not the introduction of a new physical dimension, but the use of the physically meaningful lapse/redshift field as a basis for a global gradient-topological classification of static gravitational geometry.

The gradient skeleton supplements traditional level-set representations by explicitly encoding Morse indices, stable and unstable manifolds, and heteroclinic organization.

A complete global phase theorem remains open pending full off-axis and tail-slab certification.

---

## 16. Validated-Numerics Scope Statement

The current validated-numerics scope is summarized by the following statement:

> All computer-assisted inequalities and root-isolation results were independently reproduced using outward-rounded interval arithmetic and arbitrary-precision ball arithmetic, covering the local pitchfork and fold analysis, representative critical-point completeness, Brouwer-index consistency, and asymmetric cusp scaling in full, together with the on-axis mechanism of the middle parameter slab. The global off-axis and asymptotic tail-slab portions remain open.

---

## 17. Figure Files

The following figure files are referenced by this Markdown document:

- `mp_gradient_skeleton_phase_comparison.png`
- `mp_lambda_08_stable_manifolds.png`
- `mp_three_center_phase_diagram.png`

For a repository layout, a convenient structure is:

```text
paper/
├── manuscript.md
└── figures/
    ├── mp_gradient_skeleton_phase_comparison.png
    ├── mp_lambda_08_stable_manifolds.png
    └── mp_three_center_phase_diagram.png
```

If the figures are moved into a `figures/` directory, update the Markdown image references accordingly.
