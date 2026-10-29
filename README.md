# Categorifying the Penrose polynomial using nonorientable surfaces

This repository contains supplementary computational material for joint work with
Daniel W. Collison, titled

> Categorifying the Penrose polynomial using nonorientable surfaces

The arXiv version will be available at

> https://arxiv.org/abs/2607.01632

The repository contains one main file:

1. `mobius_homology_colab_v2.ipynb`, a self-contained Google Colab notebook for
   computing the homology associated with a perfect matching graph using
   the Möbius Frobenius algebra construction.

The goal is transparency: the notebook gives an executable version of the
construction, including the local Frobenius algebra maps, the cube of
smoothings, the differentials, the check that the differential squares to zero,
and the resulting homology dimensions.

## Contact

If you find any errors in the paper **please email me**:

[dtubbenhauer@gmail.com](mailto:dtubbenhauer@gmail.com?subject=[GitHub]%web-reps)

Same goes for any errors related to this page.

## Main file

### `mobius_homology_colab_v2.ipynb`

This notebook is the recommended entry point. It is self-contained and is meant
to run directly in Google Colab.

The input is a finite graph. If the graph is trivalent, a perfect matching may be specified. In all other cases, the canonical perfect matching of the blowup of the graph is taken. The
matching edges are then smoothed, and each smoothing choice gives a single state of a hypercube. Each state then corresponds to a collection of circles.

For a state $\alpha$, the notebook forms the chain-group summand

$$
    V_\alpha = V^{\otimes k_\alpha},
$$

where $k_\alpha$ is the number of circles in the smoothing associated with
$\alpha$. The homological degree is the number of one-smoothings in
$\alpha$.

The two algebras used in the notebook are

$$
    V_{n, \xi}=\mathbb C[x,y]/(x^n-\xi, y^3-xy),
$$

with $\xi =0, 1$. For the default value `n = 1`, this specializes to

$$
    V_{1, 0} =\mathbb C[y]/(y^3) \hspace{1cm} \text{and} \hspace{1cm} V_{1, 1}=\mathbb C[y]/(y^3-y),
$$

for example. In general, the algebras can both be given the basis

$$   
B = \\{ 1, y, y^2, x, xy, xy^2, \ldots, x^{n-1}, x^{n-1}y, x^{n-1}y^2 \\}.
$$

Given $\omega = e^{2\pi i/n}$ and $p_k = \frac{1}{n} \sum_{j=0}^{n-1} \omega^{-kj}x^j$ with $k =0, \ldots, n-1$, there is also the color basis

$$   
E_{k,0} = p_k(1-x^{n-1}y^2), \hspace{1cm}
E_{k, \pm} = \frac{1}{2}p_k(x^{n-1}y^2 \pm \omega^{-k/2}y)
$$

associated with $V_{n,1}$.

The algebras are also graded $V_{n, \xi} = \oplus_{g \in \mathcal{G}_{n, \xi}}W_g$  with 

$$
\mathcal{G}_{n, \xi} =\begin{cases*}
\mathbb{Z}/n\mathbb{Z} \times \mathbb{Z}/2\mathbb{Z} & $n$ even\\
\mathbb{Z}/(n+1)\mathbb{Z} \times \mathbb{Z}/2\mathbb{Z} & $n=1$, or $n$ odd and $\xi=0$
\\
\mathbb{Z}/2\mathbb{Z} & else
\end{cases*}
$$

according to

$$
x^iy^j \in \begin{cases*}
W_{(2i+j, j)} & $n$ even\\
W_{(2i+j, j)} & $n=1$, or $n$ odd and $\xi=0$ \\
W_{j} & else
\end{cases*}.
$$

where $0 \leq i \leq n-1$ and $0 \leq j \leq2$.

For $n$ even, the Frobenius trace is chosen for both algebras so that

$$
    \epsilon(x^{n-1}y^2)=3n,
$$

and the trace is zero on the other basis elements. From this data, the notebook
constructs multiplication, the Frobenius pairing, and the corresponding
comultiplication. For $n$ odd, a shifted comultiplication is used instead.

The differential is assembled edge-by-edge in the cube of smoothings. Along a
cube edge, exactly one smoothing changes. The notebook compares the actual
circle decompositions before and after the change and applies the corresponding
local map:

- multiplication $\mu$, when two circles merge into one;
- comultiplication $\Delta$, when one circle splits into two;
- the Möbius map $m$, when one circle remains one circle but passes through a
  Möbius-type local change.

The notebook tracks the actual circle components, not just the number of
circles. This matters because the local map has to be applied to the correct
tensor factor or tensor factors inside $V^{\otimes k}$.

The signs in the differential are the usual cube signs: when changing the
$j$-th smoothing coordinate, the sign is determined by the number of
one-smoothings before coordinate $j$.

After building the differentials, the notebook checks

$$
    d_{i+1}d_i=0
$$

for all homological degrees, with a numerical tolerance of 1e^-9 for the complex number evaluation. It then computes the rank of

$$
    H^{i, *} = \ker(d_i)/\mathrm{im}(d_{i-1})
$$

numerically for each given homological degree both in total and by fixed graded (internal) degree.

The final output includes:

- the number of circles in every smoothing state;
- whether $d^2=0$;
- the homology dimensions by homological degree;
- the homology dimensions by internal degree;
- the corresponding Poincare polynomial.

Furthermore, there is an option to compute the following quantity:

$$
\text{Tr}\big(\tau_A |_{H^0}\big)
$$

where $A \subseteq \\{0, 1, \ldots, n\\}$ and $\tau_A: V_{n,1} \rightarrow V_{n,1}$ is a map defined on the color basis as follows:

$$
\tau_A:  \hspace{0.5cm} E_{a, +} \mapsto E_{a, -}, \hspace{0.5cm} E_{a, -} \mapsto E_{a, +},
$$

for all $a \in A$ and identity on the rest. The resulting quantity gives the number of $(3n-2|A|)$-face colorings of the input graph, given without a perfect matching. 

## Requirements

The notebook is designed for Google Colab. No Sage installation is required.

It uses:

1. `sympy` for algebra element construction, matrix computations, and polynomial manipulation;
2. `numpy` for fast rank and nullspace computations (SVD), and efficient matrix computations;
3. `itertools`, `collections`, and other standard-library Python modules.

Google Colab already provides `sympy` and `numpy`. If running locally instead, install them
with

```bash
python3 -m pip install sympy numpy
```

## Examples

Open `mobius_homology_colab_v2.ipynb` in Google Colab and run all cells from top to
bottom.

The main editable input cell has the form

```python
# Replace with your ribbon graph data to obtain the blowup
# cycles = [[[1, 2], [1, 1], [1, 0]], [[0, 2], [0, 1], [0, 0]]]

# Alternatively, you can input perfect matching data
cycles = [[1, 1, 1, 0], [0, 0, 0, 0]]
matching = [[0, 1]]

# Choice of parameter n
n = 1

# Boolean input: True for xi = 1 or False for xi = 0
xi = True

try:
    matching
except NameError:
    MH = MobiusHomology(cycles=cycles, n=n, xi=xi)
else:
    MH = MobiusHomology(cycles=cycles, matching=matching, n=n, xi=xi)
print('circle counts by homological degree:', MH.state_circle_counts())
print('d^2 = 0:', MH.check_d_squared())
print('homology dimensions by homological degree:', MH.homology_dimensions())
print('homology dimensions by internal degree:', MH.homology_dimensions_by_internal_degree())
if xi:
  print(f'n = {n}, xi = 1 Poincare polynomial:', MH.poincare_polynomial())
else:
  print(f'n = {n}, xi = 0 Poincare polynomial:', MH.poincare_polynomial())
```

Here:

- `cycles` encodes the clockwise ordering of edges at each vertex. Multiple edges require some more care to input; please see the corresponding paper for details;
- `matching` encodes the perfect-matching edges whose smoothings form the cube;
- `n` is the Frobenius algebra parameter;
- `xi` is either `True` or `False` for the selections $V_{n,1}$ or $V_{n, 0}$ respectively.

For this example, the expected output is

```text
circle counts by homological degree: [[2], [1]]
d^2 = 0: True
homology dimensions by homological degree: [6, 0]
homology dimensions by internal degree: {0: [3, 0], 1: [3, 0]}
n = 1, xi = 1 Poincare polynomial: 3*b*q + 3
```

This means that the cube of smoothings has one state with two circles in
homological degree zero and one state with one circle in homological degree one. The resulting homology
has dimensions

$$
    \dim H^{0, (0, 0)}=3,\qquad \dim H^{0, (1, 1)} = 3,
$$

with the rest zero, so the Poincare polynomial is

$$
    3bq+3.
$$

There is also a shortcut for inputing the data for cycle graphs $C_m$, path graphs $P_m$, complete graphs $K_m$, and complete bipartite graphs $K_{m, p}$. The second editable cell allows for the selection of both the family and specified member, before running the main computation performed with the canonical perfect matching. 

Finally, the third editable cell counts face colorings for input graphs and a given choice of set $A\subseteq \\{0,1, \ldots, n\\}$ as follows:

```python
# Replace the cycles input with your graph cycle data
cycles = [[[1, 2], [1, 1], [1, 0]], [[0, 2], [0, 1], [0, 0]]]

# Choice of parameter n
n = 1

# Choice of set A
A = []

color_num = 3*n - 2*len(A)
trace = get_face_color(cycles, n, A)
print(f"Number of {color_num}-face colorings: {int(trace)}")
```
The output for this example is 


```text
Number of 3-face colorings: 6
```

## What the notebook is doing

The main computation proceeds in five steps.

First, the notebook enumerates all smoothing states. If there are $r$
matching edges, then there are $2^r$ states. A state is a binary vector of
length $r$.

Second, for each state, the notebook constructs the actual smoothed circle
components. These components are stored combinatorially, so that the code knows
which tensor factor corresponds to which circle.

Third, the notebook constructs the cube differential. It only connects states
that differ in exactly one smoothing coordinate. This avoids the common mistake
of connecting all states in adjacent homological degrees.

Fourth, along each genuine cube edge, the notebook identifies the local change
in the circle decomposition. Depending on the change, it inserts $\mu$,
$\Delta$, or $m$, together with identity maps on all unaffected tensor
factors.

Fifth, the notebook assembles the differential matrices, checks that their
successive products are zero, and computes the homology. The differential matrix corresponding to homological degree zero is later used to compute face colorings via the trace of $\tau_A$ using standard linear algebra.

Thus the notebook computes the actual chain complex, not only the graded Euler
characteristic and not only the number of circles in each smoothing.

## Important conventions

The notebook uses the following conventions.

- Homological degree is the number of one-smoothings;
- Cube signs are determined by the number of earlier one-smoothings;
- The local map is chosen from the actual circle change:
  $2\to 1$ gives $\mu$, $1\to 2$ gives $\Delta$, and $1\to 1$ gives
  the Möbius map $m$;
- Linear algebra is over $\mathbb{C}$;
- The displayed Poincare polynomial records homological degree and the graded dimensions of the chain groups.

## Relation to the paper

The computations are intended to support examples and checks in the paper

Categorifying the Penrose polynomial using nonorientable surfaces

joint with Daniel W. Collison.

The arXiv placeholder is

```text
arXiv:2607.01632
```

The suggested citation is:

```bibtex
@misc{CollisonTubbenhauerMobiusHomologyData,
  author = {Collison, D.W. and Tubbenhauer, D.},
  title = {Code, data and more for ``{C}ategorifying the {P}enrose polynomial using nonorientable surfaces''},
  year = {2026},
  url = \url{https://github.com/dtubbenhauer/graphhomology}
}
```
In addition, the file

### `additional_material.pdf`

contains more details for some of the proofs found in the above paper; the proofs for which this occurs are signposted in both documents.

## Repository contents

Current intended contents:

```text
README.md
LICENSE
mobius_homology_colab.ipynb
additional_material.pdf
```

Additional examples, corrected notebooks, or output files may be added later.

## Erratum

Empty so far.
