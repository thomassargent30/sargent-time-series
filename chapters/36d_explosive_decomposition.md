# Decomposing an Explosive Autoregression

```{note}
This section is an addition to the book, with no source outside it. It applies the root-sorting
device of {doc}`Chapter IX <ch09_difference_equations>` to the simplest unstable process. It
connects to the martingale and bubble representations of {doc}`36_bubbles` and to Exercise 4 of
{doc}`37_exercises`.
```

Throughout this chapter the device of *solving stable roots backward and unstable roots forward*
(introduced in {doc}`Chapter IX <ch09_difference_equations>` and used in
{doc}`36a_interpreting_vars`, {doc}`36b_exact_linear_re`, and
{doc}`Chapter XIV <ch14_investment_uncertainty>`) has been a way to pick out a *bounded* or
*stationary* solution of an unstable difference equation. This short section applies the same device to
the simplest possible unstable process — a scalar explosive autoregression — and reads the result as a
statement about how to *represent a nonstationary process*.

## An explosive AR(1) and its two solutions

Consider the explosive first-order autoregression

```{math}
:label: eq-expl-ar
x_{t+1} = a\,x_t + \varepsilon_{t+1}, \qquad \varepsilon_t \sim \text{i.i.d.}(0,\sigma^2), \qquad a > 1 .
```

Iterating {eq}`eq-expl-ar` *backward* from an initial condition $x_0$ gives the **causal** solution

```{math}
:label: eq-expl-causal
x_t = a^t x_0 + \sum_{j=0}^{t-1} a^{\,t-1-j}\,\varepsilon_{j+1},
```

which expresses $x_t$ through current and past shocks but explodes as $t\to\infty$. Because $a>1$, the
root sits *outside* the region that would make {eq}`eq-expl-causal` a stationary Wold representation.
Following the chapter's recipe, we instead solve the unstable root **forward**. Rewrite
{eq}`eq-expl-ar` as $x_t = a^{-1}x_{t+1} - a^{-1}\varepsilon_{t+1}$ and define the discounted process
$y_t \equiv a^{-t}x_t$, which obeys the *stable* recursion $y_t = y_{t-1} + a^{-t}\varepsilon_t$. Since
$a>1$, the sum $\sum_{j\ge1} a^{-j}\varepsilon_j$ converges in mean square, so

```{math}
:label: eq-expl-Z
Z \equiv \lim_{t\to\infty} y_t = y_0 + \sum_{j=1}^{\infty} a^{-j}\varepsilon_j
```

is a well-defined, **time-invariant** random variable — the discounted value, from the vantage point of
the infinite future, of the initial condition and the entire shock sequence. Writing
$y_t = Z - \sum_{j>t} a^{-j}\varepsilon_j$ and multiplying by $a^{t}$ yields the decomposition

```{math}
:label: eq-expl-decomp
\boxed{\;x_t = a^{t} Z \;-\; \sum_{k=1}^{\infty} a^{-k}\,\varepsilon_{t+k}\;}
```

## Reading the decomposition

The two terms in {eq}`eq-expl-decomp` separate the process into an explosive part with a *random
amplitude* and a stationary part that looks only *forward*:

- **Explosive component $a^{t}Z$.** All of the nonstationarity is carried by a single deterministic
  explosive path $a^{t}$ scaled by the time-invariant random amplitude $Z$ of {eq}`eq-expl-Z`. Once the
  realization is fixed, $Z$ is a number; the process diverges along $a^{t}Z$ exactly as the homogeneous
  solution of {doc}`Chapter IX <ch09_difference_equations>` diverges along $\lambda^{t}$ times an
  arbitrary constant, except that here the "constant" is pinned down by the whole shock history.
- **Stationary anticipative component $\tilde x_t = -\sum_{k\ge1} a^{-k}\varepsilon_{t+k}$.** This is the
  unstable root solved forward: a geometric sum of *future* shocks, discounted by the stable factor
  $a^{-1}<1$. It is covariance stationary, with autocovariance

  ```{math}
  :label: eq-expl-cov
  \operatorname{Cov}(\tilde x_t,\tilde x_{t+h}) = \frac{\sigma^2}{a^2-1}\,a^{-|h|},
  ```

  which is precisely the covariogram of a *stable* AR(1) with root $a^{-1}$. Unlike the one-sided
  moving averages that dominate this book, $\tilde x_t$ is **anticipative** — it depends on shocks dated
  later than $t$ — so it is exactly the kind of *non-realizable* representation that
  {doc}`Chapter IX, Section 13 <ch09_difference_equations>` and
  {doc}`Chapter XIV <ch14_investment_uncertainty>` contrast with nonanticipative ones.

## Two connections

**To solving unstable roots forward.** The stationary anticipative term is the forward solution of the
unstable root, and it is exactly the object that {ref}`Exercise 4 <ex-4>` of this chapter asks the
reader to construct and then push one step further. There the deviation-from-trend
$s_t = \tilde x_t$ is shown to admit a *fundamental*, one-sided Wold representation
$s_t = (1-a^{-1}L)^{-1}u_t$ in a **new** white noise $u_t \ne \varepsilon_t$ with the *smaller* variance
$\sigma_u^2 = \sigma^2/a^2$. Thus "solving the unstable root forward" trades the original explosive
innovation $\varepsilon_t$ for a fundamental innovation, at the cost of shrinking the innovation
variance by $a^{-2}$ — the scalar analogue of the innovation-relabeling that reappears whenever a
{doc}`vector system is factored into stable and unstable roots <36a_interpreting_vars>`.

**To representing a nonstationary process.** Decomposition {eq}`eq-expl-decomp` is the explosive-root
cousin of the *unit-root* decompositions studied in {ref}`Exercise 25 <ex-25>` and its
{doc}`solution <38_exercise_solutions>`. When $a=1$, {cite:t}`beveridge1981new` split a difference-
stationary series into a random-walk **permanent** component built from *past* shocks plus a stationary
transitory part; when $a>1$, the permanent piece becomes the explosive $a^{t}Z$ and — tellingly — its
amplitude $Z$ depends on the *entire* shock sequence, while the stationary remainder looks purely
*forward*. The direction in which the "trend" accumulates information flips as the root crosses the unit
circle. The same forward-looking, martingale-like object $Z$ underlies the nonstationary bubble
solutions of {doc}`Section 36 <36_bubbles>`, where an arbitrary martingale can be added to a stationary
particular solution.

```{note}
Explosive autoregressions were studied early by Quenouille (1957) and Chow (1983); their
representation theory is developed by McCabe and Tremayne (1989), and their estimation theory by
Phillips and Magdalinos (2007). Textbook treatments of the forward (anticipative) solution appear in
Hamilton (1994). The permanent/transitory reading is due to Beveridge and Nelson (1981).
```

## References

```{bibliography}
:labelprefix: EX
:filter: key in {"beveridge1981new", "chow1983econometrics", "hamilton1994time", "mccabetremayne1989representations", "phillips2007treatment", "quenouille1957analysis"}
```
