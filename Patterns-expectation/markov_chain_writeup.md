# Markov Chain Writeup for `pattern_expectation.py`

This note explains the mathematics behind
[`pattern_expectation.py`](/Users/theodore_gerasimov/Documents/Github/LP-linal-DRP26/Patterns-expectation/pattern_expectation.py),
which computes the expected number of Bernoulli trials needed before a binary
pattern first appears.

## Recommendation

For this repository, a Markdown writeup is the better primary format.
The main reason is that the heart of the method is a sequence of definitions,
lemmas, and short proofs, and that material is easier to read, review, diff,
and cite in plain Markdown than in a notebook.

A Jupyter notebook would still be valuable as a companion artifact if the goal
is teaching through experiments. In particular, a notebook is better for:

- visualizing the transition graph,
- checking small examples numerically,
- comparing different patterns and different values of `p`.

So the best split is:

- Markdown for the canonical mathematical writeup,
- notebook later if we want an interactive demo layer.

## 1. Problem Statement

Fix a nonempty binary word

$$
w = w_1 w_2 \cdots w_m \in \{0,1\}^m.
$$

Let $(X_n)_{n \ge 1}$ be i.i.d. Bernoulli random variables with

$$
\Pr(X_n = 1) = p, \qquad \Pr(X_n = 0) = q := 1-p.
$$

Define the stopping time

$$
T_w := \inf \left\{ n \ge m : X_{n-m+1} X_{n-m+2} \cdots X_n = w \right\}.
$$

The function
[`expected_time_for_pattern`](/Users/theodore_gerasimov/Documents/Github/LP-linal-DRP26/Patterns-expectation/pattern_expectation.py#L84)
returns

$$
\mathbb{E}[T_w].
$$

The code does this by converting the pattern-matching problem into a finite
Markov chain with one absorbing state.

## 2. State Space from Prefix Matches

For each $i \in \{0,1,\dots,m\}$, state $i$ means:

> among the symbols observed so far, the longest suffix that is also a prefix
> of $w$ has length $i$.

Thus:

- state `0` means no nontrivial prefix of $w$ is currently matched,
- state `m` means the full pattern has appeared, so it is absorbing.

This is exactly what the transition builder
[`_build_transitions`](/Users/theodore_gerasimov/Documents/Github/LP-linal-DRP26/Patterns-expectation/pattern_expectation.py#L29)
encodes.

### Definition of the transition function

Let `delta(i, b)` be the next state after we are in state $i$ and observe bit
$b \in \{0,1\}$. Formally, `delta(i, b)` is the largest integer $j$ such that

$$
w_1 w_2 \cdots w_j
$$

is a suffix of the word formed by taking a length-$i$ matched prefix and then
appending $b$.

The code computes this efficiently using the KMP prefix function in
[`_prefix_function`](/Users/theodore_gerasimov/Documents/Github/LP-linal-DRP26/Patterns-expectation/pattern_expectation.py#L16).

## 3. Why the KMP Transition Construction Is Correct

### Lemma 1

For each state $i < m$ and bit $b \in \{0,1\}$, the value produced by
[`_build_transitions`](/Users/theodore_gerasimov/Documents/Github/LP-linal-DRP26/Patterns-expectation/pattern_expectation.py#L29)
is exactly `delta(i, b)`.

### Proof

Suppose we currently match the prefix $w_1 \cdots w_i$ and we append $b$.
If $w_{i+1} = b$, then the matched prefix length increases to $i+1$.

If $w_{i+1} \ne b$, then the only possible new matched length is a proper
border of $w_1 \cdots w_i$, meaning a word that is both a prefix and a suffix
of that prefix. The KMP prefix-function value `pi[i-1]` is the length of the
longest such border. If matching fails at length $i$, the next candidate must
therefore have length `pi[i-1]`. Repeating this fallback step enumerates the
possible border lengths from longest to shortest. The first candidate that can
be extended by $b$ is therefore the longest prefix of $w$ that is a suffix of
the updated observed word.

That is exactly the definition of `delta(i, b)`, so the construction is
correct. `QED`

## 4. The Markov Chain

Define a process $(Y_n)_{n \ge 0}$ by letting $Y_n$ be the matched-prefix
length after the first $n$ flips. Then $Y_0 = 0$, and for $i < m$,

$$
\Pr(Y_{n+1} = \delta(i,0) \mid Y_n = i) = q,
$$

$$
\Pr(Y_{n+1} = \delta(i,1) \mid Y_n = i) = p.
$$

Also, state $m$ is absorbing:

$$
\Pr(Y_{n+1} = m \mid Y_n = m) = 1.
$$

This is the transition rule implemented in
[`expected_time_for_pattern`](/Users/theodore_gerasimov/Documents/Github/LP-linal-DRP26/Patterns-expectation/pattern_expectation.py#L124).

### Proposition 2

$(Y_n)$ is a time-homogeneous Markov chain on the state space
$\{0,1,\dots,m\}$.

### Proof

By definition, $Y_n$ records the longest suffix of the observed sequence that
matches a prefix of $w$. Once $Y_n = i$ is known, the distribution of $Y_{n+1}$
depends only on the next independent Bernoulli flip and on `delta(i,0)` and
`delta(i,1)`. It does not depend on earlier history except through $Y_n$.

Therefore

$$
\Pr(Y_{n+1} = j \mid Y_n, Y_{n-1}, \dots, Y_0)
= \Pr(Y_{n+1} = j \mid Y_n),
$$

and the transition probabilities do not depend on $n$. So $(Y_n)$ is a
time-homogeneous Markov chain. `QED`

## 5. Connection Between Hitting Time and Pattern Occurrence

Let

$$
\tau_m := \inf\{n \ge 0 : Y_n = m\}.
$$

### Proposition 3

The hitting time $\tau_m$ is exactly the first time the pattern $w$ has appeared
as a contiguous block in the Bernoulli sequence. Hence

$$
\tau_m = T_w.
$$

### Proof

By the definition of the states, $Y_n = m$ if and only if the longest suffix of
the first $n$ observed bits that is also a prefix of $w$ has length $m$. Since
$m$ is the full length of $w$, this is equivalent to saying that the last $m$
bits observed at time $n$ are exactly $w$.

So the event `{Y_n = m}` is exactly the event that the pattern has just been
completed at time $n$. Taking the first such time shows $\tau_m = T_w$. `QED`

## 6. First-Step Equations for the Expected Waiting Time

For each state $i$, define

$$
E_i := \mathbb{E}[\tau_m \mid Y_0 = i].
$$

Clearly,

$$
E_m = 0.
$$

For each transient state $i < m$, conditioning on the next flip gives

$$
E_i = 1 + q E_{\delta(i,0)} + p E_{\delta(i,1)},
$$

with the convention that $E_m = 0$ if one of the transitions lands directly in
the absorbing state.

This is exactly the linear system assembled by the matrix code in
[`expected_time_for_pattern`](/Users/theodore_gerasimov/Documents/Github/LP-linal-DRP26/Patterns-expectation/pattern_expectation.py#L120).

### Proposition 4

The vector $(E_0,\dots,E_{m-1})$ is the unique solution of the linear system

$$
(I - Q)E = \mathbf{1},
$$

where $Q$ is the transient-to-transient submatrix of the Markov chain.

### Proof

The displayed recursion comes from first-step analysis, so any true expected
hitting-time vector must satisfy it. Written in matrix form over transient
states, this becomes

$$
E = \mathbf{1} + Q E,
$$

hence

$$
(I-Q)E = \mathbf{1}.
$$

To prove uniqueness, it is enough to show that $I-Q$ is invertible.

If the pattern is feasible under the Bernoulli law, then there is some
$\alpha > 0$ such that, from any current state, the probability that the next
$m$ flips equal the pattern exactly is $\alpha$. Indeed:

- if $0 < p < 1$, then $\alpha = p^{\#1(w)} q^{\#0(w)} > 0$,
- if $p = 0$, feasibility means $w = 00\cdots0$, in which case $\alpha = 1$,
- if $p = 1$, feasibility means $w = 11\cdots1$, in which case $\alpha = 1$.

If the next $m$ flips equal $w$, then regardless of the current state, the last
$m$ observed bits at the end of that block are exactly $w$, so the chain must
hit the absorbing state by then.

Now partition future time into disjoint blocks of length $m$. In every block,
conditional on not yet being absorbed at the start of the block, the
probability of absorption during that block is at least $\alpha$. Therefore

$$
\Pr_i(\tau_m > km) \le (1-\alpha)^k \to 0
$$

for every starting state $i < m$.

So $\tau_m < \infty$ almost surely, which means every nonabsorbing state is
transient. For a finite Markov chain, the transient matrix $Q$ then satisfies
$Q^n \to 0$, so $1$ cannot be an eigenvalue of $Q$. Hence $I-Q$ is invertible.

Therefore the linear system has a unique solution, namely the expected hitting
times. `QED`

## 7. Why the Code Returns Infinity in Degenerate Impossible Cases

The code handles the cases

- `p = 0` with at least one `1` in the pattern,
- `p = 1` with at least one `0` in the pattern,

before building the linear system; see
[`expected_time_for_pattern`](/Users/theodore_gerasimov/Documents/Github/LP-linal-DRP26/Patterns-expectation/pattern_expectation.py#L109).

### Proposition 5

In those cases, $\mathbb{E}[T_w] = \infty$.

### Proof

If `p = 0`, then every trial equals `0` almost surely, so any pattern containing
a `1` can never occur. Thus $T_w = \infty$ almost surely, hence
$\mathbb{E}[T_w] = \infty$.

The case `p = 1` and a pattern containing `0` is identical. `QED`

## 8. Worked Example: Pattern `101` for a Fair Coin

Take $w = 101$ and $p = q = 1/2$.

The states are:

- `0`: no prefix currently matched,
- `1`: the current suffix is `1`,
- `2`: the current suffix is `10`,
- `3`: the full pattern `101` has appeared.

The transition table produced by
[`_build_transitions`](/Users/theodore_gerasimov/Documents/Github/LP-linal-DRP26/Patterns-expectation/pattern_expectation.py#L29)
is:

| state | on `0` | on `1` |
| --- | --- | --- |
| `0` | `0` | `1` |
| `1` | `2` | `1` |
| `2` | `0` | `3` |
| `3` | `3` | `3` |

So the first-step equations are

$$
E_0 = 1 + \frac12 E_0 + \frac12 E_1,
$$

$$
E_1 = 1 + \frac12 E_2 + \frac12 E_1,
$$

$$
E_2 = 1 + \frac12 E_0 + \frac12 E_3,
$$

with $E_3 = 0$. Solving:

$$
E_0 - E_1 = 2,
$$

$$
E_1 - E_2 = 2,
$$

$$
E_2 = 1 + \frac12 E_0.
$$

Substituting gives

$$
E_0 = 10.
$$

This matches the implementation:

```python
expected_time_for_pattern(Fraction(1, 2), "101") == 10
```

## 9. Worked Example: Pattern `11` for a Fair Coin

For $w = 11$ and $p = q = 1/2$, the transition table is:

| state | on `0` | on `1` |
| --- | --- | --- |
| `0` | `0` | `1` |
| `1` | `0` | `2` |
| `2` | `2` | `2` |

The equations become

$$
E_0 = 1 + \frac12 E_0 + \frac12 E_1,
$$

$$
E_1 = 1 + \frac12 E_0,
$$

with $E_2 = 0$. Solving gives

$$
E_0 = 6.
$$

Again this matches the code:

```python
expected_time_for_pattern(Fraction(1, 2), "11") == 6
```

## 10. Matrix Interpretation of the Implementation

If the transient states are indexed by `0,1,...,m-1`, then the code constructs
a matrix `A` and vector `b` such that

$$
A E = b,
$$

where:

- `b` is the all-ones vector,
- `A` starts as the identity matrix,
- for each state `i`, the code subtracts `q` from column `delta(i,0)` if that
  transition remains transient,
- for each state `i`, the code subtracts `p` from column `delta(i,1)` if that
  transition remains transient.

So `A = I - Q`. The routine
[`_gaussian_elimination`](/Users/theodore_gerasimov/Documents/Github/LP-linal-DRP26/Patterns-expectation/pattern_expectation.py#L51)
solves this linear system, and the function returns `E_0`, the expected waiting
time from the initial empty-match state.

## 11. Final Summary

The code implements the following mathematical pipeline:

1. Build the KMP prefix automaton for the target pattern.
2. Interpret prefix lengths as states of a finite Markov chain.
3. Make the full-pattern state absorbing.
4. Write first-step equations for the expected hitting time.
5. Solve the linear system `(I-Q)E = 1`.
6. Return `E_0`.

This is mathematically sound because:

- the KMP fallback logic computes the correct suffix/prefix state,
- those states form a Markov chain,
- hitting the absorbing state is exactly the same as first observing the pattern,
- the first-step equations characterize the expected waiting time uniquely,
- impossible degenerate cases are correctly assigned infinite expectation.
