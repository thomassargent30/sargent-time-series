# Chapter XIV — Investment Under Uncertainty

```{note}
This chapter is Chapter XIV of Thomas J. Sargent, *Macroeconomic Theory*, 2nd ed. (Academic
Press, 1987). It builds directly on two chapters of this book: the Euler-equation and
stable-roots-backward / unstable-roots-forward machinery of
{doc}`Chapter IX <ch09_difference_equations>`, and the prediction formulas of Chapter XI —
in particular the geometric-lead formula {eq}`eq-90` and the compact vector predictor
{eq}`eq-108`. Some cross-references in the original ("Chapter I", "Chapter III", "Chapter XII",
"Chapter XVII", "the next chapter") point to other chapters of the 1987 textbook that are not
part of this Jupyter book; they are retained as prose.
```

This chapter studies aspects of the capital accumulation process in setups where firms are
uncertain about the future. Our first task is to extend our earlier study of quadratic dynamic
optimization problems to the case in which there is uncertainty about future values of the
exogenous processes facing agents. Then we shall present a simple version of Lucas and Prescott's
model of firms' investment behavior in a competitive industry. In the process we shall be able to
give a precise characterization of the concept of a rational expectations equilibrium.

## 1. Optimum Decision Rules Under a Quadratic Objective

We consider the problem: maximize (at each point in time $t$) the discounted present value

```{math}
:label: eq-14-1
v_t = E_t \sum_{j=0}^{\infty} b^j\, g(n_{t+j-1},\, n_{t+j},\, z_{t+j})
```

over stochastic processes for $\{n_{t+j}\}_{j=0}^{\infty}$ subject to $n_{t-1}=\bar n_{t-1}$ given.
Here $E_t(x)=E x\mid\Omega_t$, where $E$ is the mathematical expectation operator and $\Omega_t$ is
an information set to be specified. We assume that the discount factor obeys $0<b<1$ and that
$g(n_{t+j-1},n_{t+j},z_{t+j})$ is concave in $n_{t+j-1},n_{t+j}$. Here $z_{t+j}$ is a vector of
random variables that are exogenous to the decision maker. At time $t+j$ the decision maker will
have available an information set $\Omega_{t+j}$ on which to base his decision, with
$\Omega_t\supset\Omega_{t-1}$ for all $t$. The decision maker chooses $n_t$ and a *strategy* —
contingency plans $\tilde n_{t+1}(\cdot),\tilde n_{t+2}(\cdot),\ldots$ giving
$n_{t+1}=\tilde n_{t+1}(\Omega_{t+1})$, etc. — that make $n_{t+j}$ a function of the information
$\Omega_{t+j}$ that will be available when it must be set.[^fn-14-1]

To match the notation of {eq}`eq-14-1` with a problem that interests us, let

```{math}
:label: eq-14-2
g(n_{t-1},n_t,z_t) = (f_0+a_t)n_t - \tfrac12 f_1 n_t^2 - w_t n_t - \tfrac12 d(n_t-n_{t-1})^2
```

where $f_0,f_1,d>0$, $w_t$ is the real wage, $n_t$ is employment, and $a_t$ is a random shock to the
productivity of labor. When $w_t$ and $a_t$ are stochastic processes, the solution to {eq}`eq-14-1`
becomes a stochastic version of the demand for labor studied in
{doc}`Chapter IX <ch09_difference_equations>`.

Equating to zero the derivative of $v_t$ in {eq}`eq-14-1` with respect to $n_t$ gives one
first-order necessary condition,

```{math}
:label: eq-14-3
g_2(n_{t-1},n_t,z_t) + b E_t g_1(n_t,n_{t+1},z_{t+1}) = 0 .
```

Here $n_{t-1}$ and $z_t$ are known at $t$ while $n_{t+1}$ and $z_{t+1}$ are still random. At $t+1$
the decision maker faces a problem of the same form as {eq}`eq-14-1`, whose first-order condition
is $g_2(n_t,n_{t+1},z_{t+1}) + b E_{t+1}g_1(n_{t+1},n_{t+2},z_{t+2}) = 0$. Continuing in this way,
the plan must satisfy the system of stochastic difference equations — the **stochastic Euler
equations**

```{math}
:label: eq-14-6
g_2(n_{t+j-1},n_{t+j},z_{t+j}) + b E_{t+j}g_1(n_{t+j},n_{t+j+1},z_{t+j+1}) = 0,
\qquad j=0,1,2,\ldots
```

The transversality condition is obtained, as in {doc}`Chapter IX <ch09_difference_equations>`, by
taking the finite-$T$ first-order condition for $n_{t+T}$ and letting $T\to\infty$:

```{math}
:label: eq-14-7
\lim_{T\to\infty} b^T E_t\, g_2(n_{T-1},n_T,z_T)\, n_T = 0 .
```

**The labor-demand example.** With $g$ given by {eq}`eq-14-2` the firm maximizes

$$
v_t = E_t \sum_{j=0}^{\infty} b^j\Big\{(f_0+a_{t+j}-w_{t+j})n_{t+j} - \tfrac12 f_1 n_{t+j}^2
      - \tfrac12 d(n_{t+j}-n_{t+j-1})^2\Big\}
$$

subject to $n_{t-1}$ given. Assume the exogenous processes $\{a_{t+j}\},\{w_{t+j}\}$ are of
exponential order less than $1/\sqrt b$: for some $K>0$ and $1\le x<1/\sqrt b$,
$|E_t w_{t+j}|<K(x)^{j+t}$ and $|E_t a_{t+j}|<K(x)^{j+t}$. The Euler equations are
$f_0+a_{t+j}-w_{t+j}-f_1 n_{t+j}-d(n_{t+j}-n_{t+j-1})+dbE_{t+j}(n_{t+j+1}-n_{t+j})=0$, or

```{math}
:label: eq-14-8
b E_{t+j}n_{t+j+1} + \phi\, n_{t+j} + n_{t+j-1} = d^{-1}(w_{t+j}-a_{t+j}-f_0),
\qquad \phi = -\big[f_1 d^{-1}+(1+b)\big],\quad j=0,1,2,\ldots
```

with transversality condition

```{math}
:label: eq-14-9
\lim_{T\to\infty} E_t b^T\big\{f_0+a_{t+T}-w_{t+T}-f_1 n_{t+T}-d(n_{t+T}-n_{t+T-1})\big\}n_{t+T}=0 .
```

Equations {eq}`eq-14-8`–{eq}`eq-14-9` generalize the Euler equation and transversality condition
of the nonstochastic problem of {doc}`Chapter IX <ch09_difference_equations>`, since in the
nonstochastic case $E_{t+j}n_{t+j+1}=n_{t+j+1}$.

For convenience let $z_{t+j}=d^{-1}(w_{t+j}-a_{t+j}-f_0)$. In the nonstochastic case the solution
was $n_{t+j}=\lambda_1 n_{t+j-1}-\lambda_1\sum_{i=0}^\infty(1/\lambda_2)^i z_{t+j+i}$, where
$\lambda_1<1<\lambda_2$ solve $1+\frac{\phi}{b}L+\frac1b L^2=(1-\lambda_1 L)(1-\lambda_2 L)$. It is
natural to guess that a solution to {eq}`eq-14-8`–{eq}`eq-14-9` is

```{math}
:label: eq-14-10
n_{t+j} = \lambda_1 n_{t+j-1} - \lambda_1 \sum_{i=0}^{\infty}\Big(\tfrac{1}{\lambda_2}\Big)^i
          E_{t+j}\, z_{t+j+i}, \qquad j=0,1,2,\ldots
```

That {eq}`eq-14-10` solves the Euler equation can be verified directly. Shifting forward one period
and using the law of iterated expectations to condition on $\Omega_{t+j}$ gives

```{math}
:label: eq-14-11
E_{t+j}n_{t+j+1} = \lambda_1 n_{t+j} - \lambda_1 \sum_{i=0}^{\infty}\Big(\tfrac{1}{\lambda_2}\Big)^i
                   E_{t+j}\, z_{t+j+1+i} .
```

Substituting {eq}`eq-14-10` and {eq}`eq-14-11` into {eq}`eq-14-8` and using $-\phi=b(\lambda_1+\lambda_2)$,
$b\lambda_2=1/\lambda_1$, so that $\{b\lambda_1+\phi\}\lambda_1=-1$, collapses the identity to
$-\sum_{i=0}^\infty(1/\lambda_2)^{i+1}E_{t+j}z_{t+j+i+1}+\sum_{i=0}^\infty(1/\lambda_2)^i E_{t+j}z_{t+j+i}=z_{t+j}$,
which holds identically. Under the assumption that $z_t$ is of exponential order less than
$1/\sqrt b$, the transversality condition holds as in
{doc}`Chapter IX <ch09_difference_equations>`.

**A constructive derivation.** There is no need to guess. Write the Euler equation as

```{math}
:label: eq-14-12
b E_s n_{s+1} + \phi E_s n_s + E_s n_{s-1} = E_s z_s ,
```

or $(bB^{-2}+\phi B^{-1}+1)E_s n_{s-1}=E_s z_s$, where the operator $B$ is defined by
$B^{-j}E_{s-1}n_s=E_{s-1}n_{s+j}$ (it advances the *date being forecast* but leaves the information
set fixed).[^fn-14-2] This can be written

```{math}
:label: eq-14-13
b(\lambda_1-B^{-1})(\lambda_2-B^{-1})E_s n_{s-1} = E_s z_s ,
```

with $|\lambda_1|<1$ and $\lambda_2=1/\lambda_1 b$. Operating on both sides with
$[(\lambda_2-B^{-1})]^{-1}$ — the only legitimate forward inverse[^fn-14-3] — and setting the
constant multiplying $\lambda_2^s$ to zero to satisfy transversality, then using
$\lambda_2^{-1}=\lambda_1 b$, gives

```{math}
:label: eq-14-14
n_s - \lambda_1 n_{s-1} = -\lambda_1 \sum_{i=0}^{\infty}(\lambda_1 b)^i E_s z_{s+i} .
```

**Certainty equivalence.** The solution depends only on the conditional *means* $E_t z_{t+j}$, not
on higher moments — the "certainty equivalence" or "separation" principle possessed by quadratic
$g$'s. The problem separates into two stages: first form the forecasts $E_t z_{t+j}$; second, solve
the nonstochastic problem

$$
\max_{n_t,n_{t+1},\ldots} V_t = \sum_{j=0}^{\infty} b^j\, g(n_{t+j-1},n_{t+j},E_t z_{t+j}),
\qquad n_{t-1}\ \text{fixed}.
$$

This separation of forecasting from optimization explains why quadratic objectives are assumed in
much applied work; for general $g$ the two problems do not separate.[^fn-14-4]

## 2. Optimal Linear Policies

The decision rule {eq}`eq-14-10` sets $n_t$ as a linear function of $n_{t-1}$ and the conditional
expectations $E_t w_{t+i}, E_t a_{t+i}$. In general these conditional expectations are *nonlinear*
functions of the information in $\Omega_t$, so the optimal decision rule is nonlinear.

Suppose instead we restrict ourselves to decision rules that express $n_t$ as *linear* functions of
$n_{t-1}$ and the information in $\Omega_t$. When $g$ is quadratic, the optimal linear rule is
obtained by replacing the conditional mathematical expectations in {eq}`eq-14-10` with the
corresponding *linear least squares projections* on $\Omega_t$ (the regressions of Chapter X). In
the special case where $\{z_{t+j}\}$ is a multivariate normal process, conditional expectations are
linear and equal the projections; then linear least squares policies are optimal among *all*
decision rules. With Gaussian $\{w_t\},\{a_t\}$ the rule can, with positive probability, call for
negative employment; it is usual to assume the shock variances are small relative to the constants
so that $n_t<0$ occurs with negligible probability.

## 3. Lucas's Critique

In Chapter XI we derived a formula for the linear least squares prediction of a geometric lead of
the kind in {eq}`eq-14-14`. Suppose $z_t$ has the autoregressive representation

```{math}
:label: eq-14-15
a(L)z_t = \epsilon_t, \qquad a(L)=1-a_1 L-\cdots-a_r L^r,\quad \epsilon_t=z_t-P[z_t\mid z_{t-1},\ldots].
```

Then formula {eq}`eq-90` of Chapter XI asserts that

```{math}
:label: eq-14-16
E_t \sum_{i=0}^{\infty}\lambda_2^{-i} z_{t+i}
   = \frac{1-\lambda_1 b\, a(\lambda_1)^{-1} a(L) L^{-1}}{1-\lambda_1 b L^{-1}}\, z_t
   = a(\lambda_1 b)^{-1}\Big[1 + \sum_{j=1}^{r-1}\Big(\sum_{k=j+1}^{r}(\lambda_1 b)^{k-j}a_k\Big)L^j\Big] z_t ,
```

where $\lambda_1 b=\lambda_2^{-1}$. Substituting {eq}`eq-14-16` into {eq}`eq-14-14` gives the *pair*
of equations

```{math}
:label: eq-14-17
n_t = \lambda_1 n_{t-1} - \lambda_1\, \frac{1-\lambda_1 b\, a(\lambda_1)^{-1}a(L)L^{-1}}{1-\lambda_1 b L^{-1}}\, z_t,
\qquad a(L)z_t=\epsilon_t .
```

Equation {eq}`eq-14-17` is an optimal linear *decision rule* for employment, expressing $n_t$ as a
linear function of $n_{t-1}$ and $z_t,z_{t-1},\ldots,z_{t-r+1}$,

```{math}
:label: eq-14-18
n_t = h_1 n_{t-1} + h_2 z_t + h_3 z_{t-1} + \cdots + h_r z_{t-r+1} .
```

The form of dependence on the $z$'s depends on the parameters of $a(L)$, because the $z$'s appear
only insofar as they help predict the geometric sum $\sum_{j\ge0}(\lambda_1 b)^j z_{t+j}$. Thus the
optimal decision rule **inherits parameters from the stochastic process for $\{z_t\}$**: the
coefficients $h_2,\ldots,h_r$ depend on the law of motion {eq}`eq-14-15`. It is fruitless to search
for a single decision rule $n_t=h(n_{t-1},\Omega_t)$ that is invariant across hypothetical
environments (laws of motion for $z_t$). This principle, stated for labor supply by Gordon and
Hynes (1970), underlies Robert E. Lucas's (1976) critique of the econometric policy-evaluation
procedures that existed in 1973, which treated decision rules like {eq}`eq-14-18` as *structural*
(invariant under interventions in the process for $z_t$). Policy evaluation should instead take
into account the dependence of private decision rules on the government's choice of a policy rule —
a theme pursued in the cross-equation restrictions of
{doc}`rational expectations models <22_rational_expectations>` and
{doc}`exact linear rational expectations models <36b_exact_linear_re>`.[^fn-14-5]

## 4. Investment

We now apply these methods to Lucas and Prescott's (1971) model of investment under uncertainty,
which gives a precise illustration of a **rational expectations equilibrium**. We first study a
"perfect foresight" (nonstochastic) model and then move to the stochastic setting.

Consider an industry of $n$ identical competitive firms using capital $k_t$ to produce output
$f_0 k_t$, with $f_0>1$. The industry demand curve is

```{math}
:label: eq-14-19
p_t = A_0 - A_1 Y_t + u_t, \qquad A_0,A_1>0,
```

where $Y_t=nf_0 k_t$ is industry output and $u_t$ is a demand shock. The representative firm is a
price-taker with respect to output prices $\{p_{t+j}\}$ and prices of capital $\{J_{t+j}\}$. In the
nonstochastic case, taking known sequences $\{p_{t+j}\},\{J_{t+j}\}$ of exponential order less than
$1/\sqrt b$, the firm chooses $\{k_{t+j}\}$ to maximize

```{math}
:label: eq-14-20
v_t = \sum_{j=0}^{\infty} b^j\Big\{ p_{t+j}(f_0 k_{t+j}) - J_{t+j}(k_{t+j}-k_{t-1+j})
      - \tfrac{d}{2}(k_{t+j}-k_{t-1+j})^2 \Big\}
```

subject to $k_{t-1}$ given, where $d>0$ is an adjustment-cost coefficient. The Euler equation is

```{math}
:label: eq-14-21
p_{t+j}f_0 - J_{t+j} + bJ_{t+j+1} - d(1+b)k_{t+j} + b\,d\,k_{t+j+1} + d\,k_{t+j-1} = 0,
\qquad j=0,1,\ldots
```

which we rewrite as $b\big(1-\tfrac1b L\big)(1-L)k_{t+j+1}=\tfrac1d(J_{t+j}-bJ_{t+j+1}-p_{t+j}f_0)$.
The solution satisfying the transversality condition[^fn-14-6] is

```{math}
:label: eq-14-22
(1-L)k_{t+j+1} = \frac{-d^{-1}}{1-bL^{-1}}\big(J_{t+j+1}-bJ_{t+j+2}-p_{t+j+1}f_0\big),
```

giving the firm's rate of investment as a function of *future* values of the output price and the
price of capital.

While each firm perceives $p_t$ as independent of its own decisions, the price is determined by all
firms together through $p_t=A_0-A_1 nf_0 k_t+u_t$. We seek an **equilibrium** pair of sequences
$\{\bar p_{t+j}\},\{\bar k_{t+j}\}$ satisfying:

- (i) Given $\{\bar k_{t+j}\}$, prices clear the market: $\bar p_{t+j}=A_0-A_1 nf_0\bar k_{t+j}+u_{t+j}$.
- (ii) Facing $\{\bar p_{t+j}\}$ as a price-taker, $\{\bar k_{t+j}\}$ maximizes present value {eq}`eq-14-20`.

Substituting $A_0-A_1 f_0 nk_{t+j}+u_{t+j}$ for $p_{t+j}$ in the Euler equation {eq}`eq-14-21` —
*after* differentiating, so that the firm acts as a price taker[^fn-14-7] — gives

```{math}
:label: eq-14-23
A_0 f_0 + f_0 u_{t+j} - J_{t+j} + bJ_{t+j+1} - \{d(1+b)+A_1 f_0^2 n\}k_{t+j}
   + d\,k_{t+j-1} + d\,b\,k_{t+j+1} = 0 .
```

Writing this as $b\,k_{t+j+1}+\phi\,k_{t+j}+k_{t+j-1}=d^{-1}\{J_{t+j}-bJ_{t+j+1}-f_0 u_{t+j}-A_0 f_0\}$
with $\phi=-\big[(1+b)+\tfrac{A_1 f_0^2 n}{d}\big]$, and factoring
$1+\tfrac\phi b L+\tfrac1b L^2=(1-\lambda_1 L)(1-\lambda_2 L)$ with $\lambda_1<\tfrac1b<\lambda_2$,
the solution satisfying the firm's transversality condition is

```{math}
:label: eq-14-24
k_{t+j+1} = \lambda_1 k_{t+j} - \frac{\lambda_1 d^{-1}}{1-\lambda_2^{-1}L^{-1}}
   \big\{J_{t+j+1}-bJ_{t+j+2}-f_0 u_{t+j+1}-A_0 f_0\big\}, \qquad j=-1,0,1,\ldots
```

and the equilibrium output price is

```{math}
:label: eq-14-25
p_{t+j} = A_0 - A_1 nf_0 k_{t+j} + u_{t+j} .
```

## 5. A Digression on the Relation Between Equilibrium and Optimality

Is the equilibrium difference equation {eq}`eq-14-23` the Euler equation of an interesting maximum
problem? Consider maximizing

```{math}
:label: eq-14-26
W_t = \sum_{j=0}^{\infty} b^j\Big\{\big[A_0 f_0 nk_{t+j}-\tfrac12 A_1(f_0^2 n^2 k_{t+j}^2)+f_0 nu_{t+j}k_{t+j}\big]
      - nJ_{t+j}(k_{t+j}-k_{t+j-1}) - \tfrac12 dn(k_{t+j}-k_{t+j-1})^2\Big\}
```

subject to $k_{t-1}$ given. Its Euler equation is exactly {eq}`eq-14-23`. The bracketed term is the
area under the demand curve,

$$
\int_0^{Y_t}(A_0-A_1 x+u_t)\,dx = A_0 Y_t - \tfrac12 A_1 Y_t^2 + Y_t u_t
   = A_0 f_0 nk_t - \tfrac12 A_1 f_0^2 n^2 k_t^2 + f_0 nk_t u_t .
$$

Thus the equilibrium implicitly maximizes the social welfare criterion {eq}`eq-14-26`, which equals
discounted consumer surplus minus producer surplus. Lucas and Prescott used the observation that an
equilibrium implicitly solves a social-welfare problem to characterize equilibria in settings where
the direct method is unavailable: their device replaces a "fixed point" problem with a maximization
problem.

## 6. Investment Under Uncertainty

Now take $\{u_{t+j}\},\{J_{t+j}\}$ to be exogenous stochastic processes of exponential order less
than $1/\sqrt b$. The firm, a price-taker with respect to $\{J_{t+j}\}$ and the equilibrium price
process $\{p_{t+j}\}$, chooses $\{k_{t+j}\}$ to maximize

```{math}
:label: eq-14-27
v_t = E_t \sum_{j=0}^{\infty} b^j\Big\{p_{t+j}f_0 k_{t+j}-J_{t+j}(k_{t+j}-k_{t+j-1})
      - \tfrac{d}{2}(k_{t+j}-k_{t+j-1})^2\Big\} .
```

Paralleling Section 1, the firm's optimum plan is

```{math}
:label: eq-14-28
k_{t+j} = k_{t+j-1} - \frac1d \sum_{i=0}^{\infty} b^i E_{t+j}\big\{J_{t+j+i}-bJ_{t+j+1+i}-f_0 p_{t+j+i}\big\},
```

which agrees with {eq}`eq-14-22` in the nonstochastic case. We seek an equilibrium pair of
stochastic processes $\{\bar p_{t+j}\},\{\bar k_{t+j}\}$ such that (i) given the firm's plan,
$\bar p_{t+j}=A_0-A_1 nf_0\bar k_{t+j}+u_{t+j}$ clears the market; and (ii) facing $\{\bar p_{t+j}\}$
as a price-taker, $\{\bar k_{t+j}\}$ maximizes {eq}`eq-14-27`. Such an equilibrium is a **rational
expectations equilibrium**: firms form the forecasts of future prices appearing in {eq}`eq-14-27`
by taking conditional expectations with respect to the stochastic process that actually governs
prices. Proceeding in complete analogy with the nonstochastic analysis,

```{math}
:label: eq-14-29
k_{t+j+1} = \lambda_1 k_{t+j} - \frac{\lambda_1}{d}\sum_{i=0}^{\infty}\Big(\tfrac1{\lambda_2}\Big)^i
   E_{t+j+1}\big\{J_{t+j+1+i}-bJ_{t+j+2+i}-f_0 u_{t+j+1+i}-A_0 f_0\big\},
```

```{math}
:label: eq-14-30
p_{t+j} = A_0 - A_1 nf_0 k_{t+j} + u_{t+j} .
```

Now suppose $J_t,u_t$ are governed by

```{math}
:label: eq-14-31
a(L)J_t = \epsilon_{Jt}, \qquad g(L)u_t = \epsilon_{ut},
```

with $a(L)=1-a_1 L-\cdots-a_r L^r$, $g(L)=1-g_1 L-\cdots-g_r L^r$, and innovations
$\epsilon_{Jt}=J_t-P[J_t\mid J_{t-1},u_{t-1},\ldots]$, $\epsilon_{ut}=u_t-P[u_t\mid\ldots]$, where the
zeros of $a(z),g(z)$ exceed unity in modulus. Using formula {eq}`eq-90` of Chapter XI, {eq}`eq-14-29`
can be represented as the equilibrium motion of capital

```{math}
:label: eq-14-32
\begin{aligned}
k_{t+j+1} = \lambda_1 k_{t+j} - \lambda_1 d^{-1}\Big\{
  &\Big[\tfrac{1-\lambda_1 b\, a(\lambda_1)^{-1}a(L)L^{-1}}{1-\lambda_1 b L^{-1}}\Big]J_{t+j+1}
  - b\Big[\tfrac{L^{-1}I-L^{-1}a(b\lambda_1)^{-1}a(L)}{1-\lambda_1 b L^{-1}}\Big]J_{t+j+1} \\
  &- f_0\Big[\tfrac{1-\lambda_1 b\, g(\lambda_1)^{-1}g(L)L^{-1}}{1-\lambda_1 b L^{-1}}\Big]u_{t+j+1}
  - A_0 f_0/(1-\lambda_1 b)\Big\} .
\end{aligned}
```

## 7. Supply, Demand and Identification

Equations {eq}`eq-14-29` and {eq}`eq-14-32` can be used to derive a dynamic "supply curve." For
convenience specialize {eq}`eq-14-31` to $g(L)=I$ (so $u_t$ is white noise) and $a(L)=(1-\rho L)$,
$|\rho|<1$ (so $J_t$ is first-order autoregressive). The firm's stock obeys {eq}`eq-14-28`:

```{math}
:label: eq-14-28b
k_t = k_{t-1} - d^{-1}\sum_{j=0}^{\infty} b^j E_t\{J_{t+j}-bJ_{t+j+1}\} + d^{-1}f_0\sum_{j=0}^{\infty} b^j E_t p_{t+j} .
```

Let $K_t=nk_t$ be aggregate capital, $\bar f_0=f_0\lambda_1 d^{-1}n$, $\bar A_0=A_0/(1-\lambda_1 b)$.
The aggregate capital stock follows

```{math}
:label: eq-14-33
K_t = \lambda_1 K_{t-1} - \lambda_1 d^{-1}n\Big\{\tfrac{1-b\rho}{1-\lambda_1 b\rho}\Big\}J_t
      + \bar f_0(u_t+\bar A_0) .
```

Using {eq}`eq-14-33` and $J_t=\rho J_{t-1}+\epsilon_{Jt}$, the law of motion of the vector
$(K_t,J_t,1)$ is

```{math}
:label: eq-14-34
\begin{pmatrix} K_{t+1}\\ J_{t+1}\\ 1\end{pmatrix}
=\begin{pmatrix}
\lambda_1 & \dfrac{-\lambda_1\rho d^{-1}n(1-b\rho)}{1-\lambda_1 b\rho} & \bar f_0\bar A_0\\[2ex]
0 & \rho & 0\\[1ex]
0 & 0 & 1
\end{pmatrix}
\begin{pmatrix} K_t\\ J_t\\ 1\end{pmatrix}
+\begin{pmatrix}
\dfrac{-\lambda_1 d^{-1}n(1-b\rho)}{1-\lambda_1 b\rho}\,\epsilon_{Jt+1}+\bar f_0 u_{t+1}\\[2ex]
\epsilon_{Jt+1}\\[1ex]
0
\end{pmatrix}
```

or $x_{t+1}=A x_t+\epsilon_{t+1}$, where $x_t=(K_t,J_t,1)^T$. From Chapter XI we have the compact
predictor {eq}`eq-108`,

```{math}
:label: eq-14-36
E_t \sum_{j=0}^{\infty}\mu^j x_{t+j} = (I-\mu A)^{-1}x_t ,
```

valid for any scalar $\mu$ whose product with the modulus of the largest eigenvalue of $A$ is less
than one. Applying {eq}`eq-14-36` to {eq}`eq-14-34` gives

```{math}
:label: eq-14-37
E_t \sum_{j=0}^{\infty}\mu^j \begin{pmatrix}K_{t+j}\\ J_{t+j}\\ 1\end{pmatrix}
=\begin{pmatrix}
\dfrac{1}{1-\mu\lambda_1}K_t
   - \dfrac{\mu n\lambda_1\rho d^{-1}(1-b\rho)}{(1-\lambda_1 b\rho)(1-\mu\lambda_1)(1-\rho\mu)}J_t
   + \dfrac{\bar f_0\bar A_0}{(1-\mu\lambda_1)(1-\mu)}\\[2.5ex]
\dfrac{1}{1-\rho\mu}J_t\\[2ex]
\dfrac{1}{1-\mu}
\end{pmatrix} .
```

Returning to {eq}`eq-14-28b` with $a(L)=(1-\rho L)$ and substituting the demand curve
$p_{t+j}=A_0-A_1 f_0 K_{t+j}+u_{t+j}$ for $p_{t+j}$,

```{math}
:label: eq-14-38
k_t = k_{t-1} - d^{-1}J_t + \frac{d^{-1}f_0 A_0}{1-b} + d^{-1}f_0 u_t
      - d^{-1}f_0^2 A_1 \sum_{j=0}^{\infty} b^j E_t K_{t+j} .
```

Using {eq}`eq-14-37` to evaluate the geometric sum, eliminating $K_t$ via the demand curve, and
using $y_t=f_0 k_t$, gives the **dynamic supply curve**

```{math}
:label: eq-14-39
\begin{aligned}
y_t = y_{t-1} &- d^{-1}f_0\Big\{1-\tfrac{f_0^2 A_1 b\lambda_1\rho d^{-1}n}{(1-\lambda_1 b\rho)(1-b\lambda_1)}\Big\}J_t
      + \Big[\tfrac{d^{-1}f_0^2 A_0}{1-b}-\tfrac{d^{-1}f_0^3 A_1 b\bar f_0\bar A_0 n}{(1-b\lambda_1)(1-b)}
        -\tfrac{d^{-1}f_0^2 A_0}{1-b\lambda_1}\Big]\\
      &- \tfrac{d^{-1}f_0^2 b\lambda_1}{1-b\lambda_1}u_t + \tfrac{d^{-1}f_0^2}{1-b\lambda_1}p_t ,
\end{aligned}
```

expressing the firm's output as a function of lagged output and current factor and output prices.
**All** of the demand-curve parameters appear in the supply curve — directly and through their
influence on $\lambda_1$ — and the demand disturbance $u_t$ appears on the right side. Any variable
that helps predict future prices appears in firms' supply curve as an "information variable." This
subverts the exclusion restrictions ordinarily relied upon to identify a supply curve. The
identifying restrictions that *are* available come from the *cross-equation restrictions* linking
the equilibrium law of motion {eq}`eq-14-32` to the laws of motion {eq}`eq-14-31` — the common
appearance of $a(L),g(L)$ in both — exploited by Taylor (1979, 1980) and Hansen and Sargent (1980).

## 8. Investment Under Uncertainty and an Externality

Now let the technology exhibit an externality,

```{math}
:label: eq-14-40a
y_t = (f_0 k_t + f_1 K_t), \qquad f_1\neq 0,
```

where $K_t=nk_t$ is the aggregate capital stock (a linear version of Romer 1983). Define the
one-period rental on capital $w_t=J_t-bE_t J_{t+1}$, equivalently
$J_t=\sum_{j=0}^\infty b^j E_t w_{t+j}$.[^fn-14-def] The representative firm maximizes
$E_0\sum_{t=0}^\infty b^t\{p_t(f_0 k_t+f_1 K_t)-w_t k_t-\tfrac d2(k_t-k_{t-1})^2\}$ subject to
$p_t=A_0-A_1(f_0+nf_1)K_t+u_t$ and the laws of motion for $w_t,u_t$. Its Euler equation is
$E_t\{p_t f_0-w_t-d(k_t-k_{t-1})+bd(k_{t+1}-k_t)\}=0$. Substituting the demand curve and multiplying
by $n$ gives the equilibrium difference equation

```{math}
:label: eq-14-42
-bE_t K_{t+1} + \big[(1+b)+A_1 nd^{-1}(f_0^2+nf_1 f_0)\big]K_t - K_{t-1}
   = d^{-1}\big[A_0 f_0 n + f_0 nu_t - nw_t\big] .
```

The social planning problem maximizes expected discounted consumer plus producer surplus,

```{math}
:label: eq-14-43
E_0 \sum_{t=0}^{\infty} b^t\Big\{A_0(f_0+nf_1)K_t - \tfrac{A_1}{2}[(f_0+nf_1)K_t]^2
   + u_t(f_0+nf_1)K_t - w_t K_t - \tfrac{d}{2n}(K_t-K_{t-1})^2\Big\},
```

whose Euler equation is

```{math}
:label: eq-14-44
-bE_t K_{t+1} + \Big[(1+b)+\tfrac{A_1 n(f_0+nf_1)^2}{d}\Big]K_t - K_{t-1}
   = d^{-1}\big[nA_0(f_0+nf_1)+n(f_0+nf_1)u_t-nw_t\big] .
```

Comparing {eq}`eq-14-42` and {eq}`eq-14-44`, when $f_1\neq 0$ the competitive equilibrium is *not*
optimal relative to the welfare criterion {eq}`eq-14-43`: the externality drives a wedge.

**Correcting the externality with a tax.** Let the government levy $\tau_t k_t$ with

```{math}
:label: eq-14-45
\tau_t = \alpha_0 + \alpha_1 u_t + \alpha_2 w_t + \alpha_3 k_t .
```

Under this schedule the firm's Euler equation, after substituting the demand curve and multiplying
by $n$, becomes

```{math}
:label: eq-14-46
-bE_t K_{t+1} + \big[(1+b)+A_1 nd^{-1}(f_0^2+nf_1 f_0)+2\alpha_3 d^{-1}\big]K_t - K_{t-1}
   = d^{-1}\big[(nA_0 f_0-n\alpha_0)+(nf_0-n\alpha_1)u_t-(n+\alpha_2 n)w_t\big] .
```

Equating {eq}`eq-14-46` to the planner's {eq}`eq-14-44` requires the government to set

```{math}
:label: eq-14-47
\alpha_0 = -A_0 nf_1, \quad \alpha_1 = -nf_1, \quad \alpha_2 = 0,
\quad \alpha_3 = \tfrac{A_1 n}{2}\big[(f_0+nf_1)^2-(f_0^2+nf_1 f_0)\big] .
```

The optimal tax {eq}`eq-14-47` requires the government to know the parameters of preferences
$(A_0,A_1)$, technology $(f_0,f_1)$, and industry structure $(n)$. The reason rational expectations
econometrics aims to estimate the "deep" parameters of preferences and technologies is precisely
that these are the parameters needed to derive optimal policy interventions.

## 9. Conclusions

The rational expectations competitive equilibrium has the attractive property that firms in the
industry forecast the output price optimally: firms forecast prices as well as the economist
modelling them. The output price is endogenous, influenced by firms' behavior in light of their
forecasts. There is a mapping from firms' *perceived* law of motion for the output price to the
*actual* law of motion; a rational expectations equilibrium is a **fixed point** of this mapping.
Linear Lucas–Prescott models have served as the basis for econometric implementation (Sargent 1981;
Hansen and Sargent 1980), multi-factor extensions (Hansen and Sargent 1981; Eichenbaum 1983), and
studies of government interventions and exhaustible resources (Eckstein and Eichenbaum 1985a,b;
Hansen, Epple, and Roberds 1985; Townsend 1983).

## Exercises

**1.** Assume that

$$
J_{t+1}=\alpha J_t+\eta_{t+1},\ \ |\alpha|<1/b, \qquad
u_{t+1}=\beta u_t+\varepsilon_{t+1},\ \ |\beta|<1/b,
$$

where $E_t\eta_{t+1}=\bar\eta$ and $E_t\varepsilon_{t+1}=\bar\varepsilon$. Use this information and
{eq}`eq-14-29` to calculate a "reduced form" for investment of the form

$$
k_{t+j+1}-\lambda_1 k_{t+j}=\gamma_0+\gamma_1 J_{t+1+j}+\gamma_2 u_{t+1+j},
$$

giving explicit formulas for $\gamma_0,\gamma_1,\gamma_2$.

```{admonition} Solution to Exercise 1
:class: dropdown

For an AR(1) with drift, $E_s J_{s+i}=\alpha^i J_s+\bar\eta\frac{1-\alpha^i}{1-\alpha}$ and likewise for
$u$. Set $s=t+j+1$ and evaluate the geometric sums in {eq}`eq-14-29` using $\lambda_2^{-1}=\lambda_1 b$.
The $J$-terms combine as $\sum_i\lambda_2^{-i}E_s(J_{s+i}-bJ_{s+1+i})=\frac{1-b\alpha}{1-\alpha\lambda_1 b}J_s+\text{const}$,
and $-f_0\sum_i\lambda_2^{-i}E_s u_{s+i}=-\frac{f_0}{1-\beta\lambda_1 b}u_s+\text{const}$. Since
$k_{t+j+1}-\lambda_1 k_{t+j}=-\tfrac{\lambda_1}{d}\sum_i\lambda_2^{-i}E_s\{J_{s+i}-bJ_{s+1+i}-f_0 u_{s+i}-A_0 f_0\}$,

$$
\gamma_1=-\frac{\lambda_1}{d}\,\frac{1-b\alpha}{1-\alpha\lambda_1 b},
\qquad
\gamma_2=\frac{\lambda_1 f_0}{d\,(1-\beta\lambda_1 b)},
$$

$$
\gamma_0=-\frac{\lambda_1}{d}\Big\{\frac{\bar\eta(1-b)}{(1-\alpha)}\Big[\tfrac{1}{1-\lambda_1 b}-\tfrac{1}{1-\alpha\lambda_1 b}\Big]
   -\frac{f_0\bar\varepsilon}{1-\beta}\Big[\tfrac{1}{1-\lambda_1 b}-\tfrac{1}{1-\beta\lambda_1 b}\Big]
   -\frac{A_0 f_0}{1-\lambda_1 b}\Big\}.
$$

Investment responds to the current factor price and demand shock with coefficients that depend on the
*persistence* parameters $\alpha,\beta$ (through the geometric-lead weights) — the reduced-form
coefficients are not invariant to the $J,u$ processes, an instance of Lucas's critique.
```

**2.** *(Certainty equivalence principle.)* Let $x$ be a random variable with density $g(x)$, and let
$\alpha$ be a parameter set by a decision-maker. Let $f(x,\alpha)$ be concave and twice continuously
differentiable. Consider *Problem 1*: choose $\alpha$ to maximize $E f(x,\alpha)=\int f(x,\alpha)g(x)\,dx$.

&nbsp;&nbsp;**A.** Find the first-order condition for choosing $\alpha$.

&nbsp;&nbsp;**B.** Suppose $f(x,\alpha)=(x,\alpha)A(x,\alpha)'+(x,\alpha)B$, where $B$ is $2\times1$,
$(x,\alpha)$ is $1\times2$, and $A$ is a $2\times2$ negative definite matrix. Prove that in this case
choosing $\alpha$ to solve Problem 1 gives the same $\alpha$ as *Problem 2*: choose $\alpha$ to
maximize $f(Ex,\alpha)$.

```{admonition} Solution to Exercise 2
:class: dropdown

**A.** Differentiating under the integral, $\frac{d}{d\alpha}\int f(x,\alpha)g(x)\,dx=\int
f_\alpha(x,\alpha)g(x)\,dx=E f_\alpha(x,\alpha)=0$. Concavity gives the second-order condition.

**B.** Write $A=\begin{pmatrix}a_{11}&a_{12}\\a_{21}&a_{22}\end{pmatrix}$, $B=(b_1,b_2)'$, so

$$
f(x,\alpha)=a_{11}x^2+(a_{12}+a_{21})x\alpha+a_{22}\alpha^2+b_1 x+b_2\alpha,
\qquad
f_\alpha=(a_{12}+a_{21})\,x+2a_{22}\,\alpha+b_2 .
$$

Because $f_\alpha$ is **linear in $x$**, $E f_\alpha(x,\alpha)=(a_{12}+a_{21})(Ex)+2a_{22}\alpha+b_2=
f_\alpha(Ex,\alpha)$. Hence the Problem-1 condition $E f_\alpha=0$ and the Problem-2 condition
$f_\alpha(Ex,\alpha)=0$ are identical, giving

$$
\alpha^*=-\frac{(a_{12}+a_{21})Ex+b_2}{2a_{22}} ,
$$

a maximizer since $A$ negative definite implies $a_{22}<0$. Only the mean $Ex$ enters — the
certainty-equivalence (separation) principle. It fails for non-quadratic $f$, where $f_\alpha$ is
nonlinear and $Ef_\alpha(x,\alpha)\neq f_\alpha(Ex,\alpha)$.
```

**3.** *(Eckstein 1983.)* A large number $n$ of identical farms produce corn; each maximizes

$$
E_0\sum_{t=0}^{\infty} b^t\Big\{p_{t+1}y_{t+1}-c_0 a_t-\tfrac{c_1}{2}a_t^2-c_2 a_t a_{t-1}-w_t a_t\Big\}
$$

subject to $a_{-1}$ given, with $c_0,c_1,c_2>0$, $0<b<1$, output $y_{t+1}=fa_t$ ($f>0$), and demand
$p_{t+1}=\beta_0-\beta_1 Y_{t+1}+u_{t+1}$ where $Y_{t+1}=ny_{t+1}$ and $u_t=\rho u_{t-1}+\varepsilon_t$,
$|\rho|<1$. Fertilizer cost $w_t=d_0\eta_t+d_1\eta_{t-1}$ ($\eta$ white noise). At $t$ the farm sees
$w_t,\ldots,u_t,\ldots,a_{t-1}$ and $A_{t-1}=na_{t-1}$.

&nbsp;&nbsp;**A.** Carefully *define* a rational expectations equilibrium.

&nbsp;&nbsp;**B.** Describe how to compute it; get as far as you can for $\beta_1=10^{-6}$, $n=10^6$,
$f=1$, $b=1$, $c_1=16$, $c_2=4$.

&nbsp;&nbsp;**C.** Describe the effect on the equilibrium law of motion for $A_t$ of changing the
fertilizer-cost process to $w_t=g_0\varepsilon_t+g_1\varepsilon_{t-1}+g_2\varepsilon_{t-2}$.

```{admonition} Solution to Exercise 3
:class: dropdown

**A.** A rational expectations equilibrium is a stochastic process $\{a_t\}$ (equivalently $\{p_t\}$)
such that (i) taking the price process $\{p_{t+1}\}$ and the exogenous $w_t,u_t$ as given, the
representative farm's $\{a_t\}$ maximizes its objective; and (ii) the price process satisfies the
demand curve $p_{t+1}=\beta_0-\beta_1 nf a_t+u_{t+1}$ when $a_t$ is the representative farm's choice —
so farms' forecasts $E_t p_{t+1}$ are conditional expectations under the *actual* equilibrium price
process.

**B.** The farm's Euler equation (differentiate w.r.t. $a_t$; it enters the date-$t$ revenue
$p_{t+1}fa_t$ and costs and the date-$(t{+}1)$ cost $c_2 a_{t+1}a_t$) is

$$
bc_2\,E_t a_{t+1}+c_1 a_t+c_2 a_{t-1}=f\,E_t p_{t+1}-c_0-w_t .
$$

Imposing equilibrium $E_t p_{t+1}=\beta_0-\beta_1 nf a_t+\rho u_t$ gives

$$
bc_2\,E_t a_{t+1}+(c_1+\beta_1 nf^2)a_t+c_2 a_{t-1}=f\beta_0+f\rho u_t-c_0-w_t .
$$

Factor the characteristic polynomial, solve the stable root backward / unstable root forward (as in
the text), and use the AR(1) forecast of $u$ and the MA(1) forecast of $w$. With the given numbers
$\beta_1 nf^2=10^{-6}\cdot10^{6}\cdot1=1$, so $c_1+\beta_1 nf^2=17$, and (dividing by $bc_2=4$) the
homogeneous characteristic equation is $z^2+\tfrac{17}{4}z+1=(z+\tfrac14)(z+4)=0$, with roots
$-\tfrac14$ and $-4$ (a reciprocal pair since $b=1$). The stable root is $\lambda_1=-\tfrac14$, so

$$
a_t=-\tfrac14\,a_{t-1}+\text{(geometric-lead terms in }E_t u_{t+i},E_t w_{t+i}\text{)} .
$$

The negative root makes acreage oscillate period-to-period (a cobweb-like alternation), damped since
$|\lambda_1|<1$.

**C.** Changing $w$ from MA(1) to MA(2) changes the term $\sum_{i\ge0}(\lambda_1 b)^i E_t w_{t+i}$ in
the decision rule: for MA(1) only $E_t w_{t+1}\neq0$, whereas for MA(2) both $E_t w_{t+1}$ and
$E_t w_{t+2}$ are nonzero, adding another lag of the innovation to the rule. Hence the equilibrium law
of motion for $A_t$ *inherits the new cost-process parameters* $g_0,g_1,g_2$ — the decision rule is
not invariant to the $w$-process (Lucas's critique).
```

**4.** A small country produces bananas competitively with free entry. The world price
$p_t$ is exogenous with Wold representation $p_t=c(L)\varepsilon_t$ ($c(L)^{-1}$ one-sided,
square-summable). Output $y_t=f(L)n_t$, $f(L)=\sum_{j\ge0}f_j L^j$, $n_t$ employment. Each firm takes
$w_t$ parametrically, but the country wage obeys $w_t=\beta_0+\beta_1 N_t$ ($N_t$ total employment).
The representative firm faces $p_t$ and $w_t$ parametrically and solves
$\max E_0\sum_{t=0}^{\infty}\{p_t f(L)n_t-w_t n_t\}$ (constant returns to scale).

&nbsp;&nbsp;**A.** Find the "marginal expected present value" of employing an additional worker at $t$.

&nbsp;&nbsp;**B.** Impose free entry (zero expected present value) and derive the equilibrium condition
$w_t=E_t\sum_{j=0}^{\infty}h_j p_{t+j}$, giving the $h_j$.

&nbsp;&nbsp;**C.** With $f(L)=1/(1-\lambda L)$ and $c(L)=1/(1-\rho_1 L-\rho_2 L^2)$, derive
$w_t=\sum_{j=0}^{\infty}g_j p_{t-j}$, giving $g_j$ in terms of $\lambda,\rho_1,\rho_2$.

&nbsp;&nbsp;**D.** How would this change if a banana cartel changed $c(L)$ to $c(L)=1+0.99L$? How does
this illustrate Lucas's critique?

```{admonition} Solution to Exercise 4
:class: dropdown

**A.** A worker hired at $t$ raises output at $t+j$ by $f_j$ (through $f(L)$), sold at $p_{t+j}$. The
marginal expected present value is therefore $E_t\sum_{j=0}^{\infty}f_j\,p_{t+j}$ (with a discount
factor $b$: $E_t\sum_j b^j f_j p_{t+j}$).

**B.** Free entry with constant returns drives the marginal worker's value to the wage:

$$
w_t=E_t\sum_{j=0}^{\infty}h_j\,p_{t+j}, \qquad h_j=f_j\ \ (\text{or }b^j f_j\text{ with discounting}).
$$

**C.** With $f(L)=1/(1-\lambda L)$, $f_j=\lambda^j$, so $w_t=E_t\sum_j(b\lambda)^j p_{t+j}$ (take
$b=1$: $\sum_j\lambda^j$). Since $p_t$ is AR(2) with $a(L)=1-\rho_1 L-\rho_2 L^2$, the geometric-lead
formula {eq}`eq-90` gives (with $b=1$)

$$
w_t=\frac{1+\lambda\rho_2 L}{1-\rho_1\lambda-\rho_2\lambda^2}\,p_t
   =g_0 p_t+g_1 p_{t-1},\qquad
g_0=\frac{1}{1-\rho_1\lambda-\rho_2\lambda^2},\quad
g_1=\frac{\lambda\rho_2}{1-\rho_1\lambda-\rho_2\lambda^2},
$$

and $g_j=0$ for $j\ge2$: the equilibrium wage is a two-term distributed lag on the price.

**D.** If the cartel makes $p_t=(1+0.99L)\varepsilon_t$ (MA(1)), then $E_t p_{t+1}=0.99\varepsilon_t$
and $E_t p_{t+j}=0$ for $j\ge2$, so $w_t=p_t+0.99\lambda\,\varepsilon_t$ with
$\varepsilon_t=c(L)^{-1}p_t=\sum_j(-0.99)^j p_{t-j}$, giving an *entirely different*, sign-alternating
distributed lag $\{g_j\}$. The wage–price relation (5) is not invariant to the price process $c(L)$:
the "structural" $g_j$ change with the exogenous regime — Lucas's critique.
```

**5.** *(Stabilizing prices vs. quantities.)* An industry of $n$ identical firms has $y_t=fk_t$,
$f>0$. The representative firm maximizes
$E_0\sum_{t=0}^{\infty}\beta^t\{p_t y_t-(w_t+\tau_t)k_t-\tfrac{d}{2}(k_t-k_{t-1})^2\}$, $0<\beta<1$,
with demand $p_t=A_0-A_1 Y_t+u_t$, $Y_t=fK_t$, $\tau_t$ a tax-subsidy on capital,
$w_t=\bar\omega+\varepsilon_{wt}$ ($\varepsilon_{wt}$ serially independent, mean zero) and $u_t$
serially independent mean zero.

&nbsp;&nbsp;**A.** With $\tau_t\equiv0$, define a rational expectations competitive equilibrium and
display the equilibrium law of motion for $K_t$.

&nbsp;&nbsp;**B.** Let the government set $\tau_t$ as a linear function of
$\{u_{t-1},\ldots,\varepsilon_{w,t-1},\ldots\}$ (zero constant). Find the feedback rule that minimizes
the stationary variance of $Y_t$; then the rule that minimizes the stationary variance of $p_t$.

&nbsp;&nbsp;**C.** Under what more general assumption about $\{u_t\}$ would the two objectives give
*different* rules for $\tau_t$?

&nbsp;&nbsp;**D.** Is minimizing the stationary variance of price a good policy goal? Of quantity?

```{admonition} Solution to Exercise 5
:class: dropdown

**A.** A rational expectations competitive equilibrium is a process $\{K_t\}$ (and price) such that
each firm, price-taking on $\{p_t\}$ and using the equilibrium law of motion for $K_t$, maximizes its
value, and the market clears $p_t=A_0-A_1 fK_t+u_t$. The firm's Euler equation
$E_t\{p_t f-w_t-d(k_t-k_{t-1})+bd(k_{t+1}-k_t)\}=0$, with the demand curve substituted and multiplied
by $n$, gives (as in {eq}`eq-14-42` with $f_1=0$)

$$
-bE_tK_{t+1}+\big[(1+b)+A_1 nf^2 d^{-1}\big]K_t-K_{t-1}=d^{-1}[A_0 fn+fnu_t-nw_t].
$$

Factoring $1+\tfrac\phi b L+\tfrac1b L^2=(1-\lambda_1 L)(1-\lambda_2 L)$, $\lambda_1<1$, and solving the
stable root backward, and because $u_t,w_t$ are **white noise** (so $E_t u_{t+i}=E_t w_{t+i}=0$ for
$i\ge1$ and the geometric-lead sums collapse to their $i=0$ terms),

$$
K_t=\lambda_1 K_{t-1}+\pi_0+\pi_u u_t+\pi_w w_t,\qquad
\pi_u=\tfrac{\lambda_1 fn}{d},\ \ \pi_w=-\tfrac{\lambda_1 n}{d},
$$

(with $\pi_0$ the constant from $A_0 f$).

**B.** Since $\tau_t$ depends only on *past* shocks, it is predetermined at $t$ and enters like a
known forcing term, shifting the constant/feedback part of $K_t$ but not offsetting the
contemporaneous white-noise innovations $u_t,w_t$ that hit $K_t$. Write the closed-loop law of motion
$K_t=g_0+g_1(L)u_t+g_2(L)\varepsilon_{wt}$ (with the feedback folded in). *Minimizing
$\operatorname{var}(Y_t)=f^2\operatorname{var}(K_t)$* selects the feedback that flattens the response
of $K$ to the shocks; *minimizing $\operatorname{var}(p_t)$* with $p_t=A_0-A_1 fK_t+u_t$ has the
*extra direct term $u_t$*, so it trades off the $K$-response against the direct demand shock. With
white-noise $u$, the predictable part $\tau$ can act on is only the *lagged* shocks, so the two rules
coincide except for how they treat the contemporaneous $u_t$ that appears directly in $p$ but only
through $K$ in $Y$.

**C.** If $\{u_t\}$ is **serially correlated**, then part of $u_t$ is predictable from $\Omega_{t-1}$,
so $\tau_t$ (a function of past shocks) *can* offset it. That predictable component enters price
*directly* (coefficient $1$) but quantity only through $K$ (coefficient $-A_1 f\cdot$feedback). Hence
the price- and quantity-stabilizing feedback rules differ whenever $u$ is serially correlated (they
coincide only in the white-noise case).

**D.** Neither is obviously desirable. Welfare in this model is the surplus criterion (cf.
{eq}`eq-14-26`), not the variance of $p$ or of $Y$. Stabilizing price can *increase* quantity
variance and vice versa; minimizing either variance is a proxy that need not track the surplus
objective. (This is the discrete-time rational-expectations analogue of Weitzman's "prices vs.
quantities.")
```

**6.** An industry of $n$ firms faces demand $p_t=A_0-A_1 Y_t+u_t$ ($u_t$ serially uncorrelated, mean
zero) with representative-firm cost
$c_t=c_0+c_1 y_t+\tfrac{c_2}{2}y_t^2+\tfrac{c_3}{2}(y_t-y_{t-1})^2+\tfrac{c_4}{2}(Y_t-Y_{t-1})^2+J_t y_t+\tau_t y_t$,
where $Y_t=ny_t$, $J_t=\lambda J_{t-1}+\varepsilon_t$, and the tax rule is
$\tau_t=\delta_0+\delta_1 y_t+\delta_2 y_{t-1}$ (each firm recognizes the dependence of $\tau_t$ on its
own $y_t$). The term $\tfrac{c_4}{2}(Y_t-Y_{t-1})^2$ is an industry-wide adjustment cost — an
externality. Firms maximize $E\sum_{t=0}^{\infty}\beta^t\{p_t y_t-c_t\}$.

&nbsp;&nbsp;**A.** Define a rational expectations competitive equilibrium.
&nbsp;&nbsp;**B.** Compute it.
&nbsp;&nbsp;**C.** For the "social planning" problem (maximize expected consumer minus producer
surplus, omitting the transfer $\tau_t y_t$), does the competitive equilibrium solve it for arbitrary
$(\delta_0,\delta_1,\delta_2)$?
&nbsp;&nbsp;**D.** Find $(\delta_0,\delta_1,\delta_2)$ that make it do so.

```{admonition} Solution to Exercise 6
:class: dropdown

**A.** A process $\{y_t\}$ (and price) such that each firm — taking the price process $\{p_t\}$ and the
industry aggregate $\{Y_t\}$ in the externality term as given (via a perceived law of motion), while
recognizing $\tau_t$'s dependence on its own $y_t$ — maximizes its value, and the market clears
$p_t=A_0-A_1 ny_t+u_t$ with rational forecasts.

**B.** The firm's Euler equation (differentiate w.r.t. $y_t$; note $Y_t$ is external, $\tau_t y_t$ is
internal) is

$$
E_t\big[p_t-c_1-(c_2+2\delta_1)y_t-c_3(y_t-y_{t-1})+\beta c_3(y_{t+1}-y_t)-J_t-\delta_0-\delta_2 y_{t-1}-\beta\delta_2 y_{t+1}\big]=0.
$$

Substituting $p_t=A_0-A_1 ny_t+u_t$ (after differentiating) gives a second-order stochastic Euler
equation; factor and solve the stable root backward, using the AR(1) forecast of $J$ and white-noise
$u$, to get $y_t=\lambda_1 y_{t-1}+(\text{terms in }J_t,u_t)$.

**C.** No. The planner internalizes the externality: differentiating the true surplus criterion (which
contains $\tfrac{c_4}{2}(Y_t-Y_{t-1})^2$ with $Y_t=ny_t$) produces extra terms
$-c_4 n(y_t-y_{t-1})+\beta c_4 n(y_{t+1}-y_t)$ that the competitive firm omits. So for arbitrary
$\delta$ the two Euler equations differ — the $c_4$ externality is a wedge.

**D.** Choose the tax so its marginal contribution reproduces the marginal external cost. Matching the
tax terms $-\delta_0-2\delta_1 y_t-\delta_2 y_{t-1}-\beta\delta_2 E_t y_{t+1}$ to the missing external
terms $-c_4 n(y_t-y_{t-1})+\beta c_4 n(y_{t+1}-y_t)$ gives

$$
\delta_1=\tfrac{c_4 n}{2}\ (\text{so }2\delta_1=c_4 n),\qquad \delta_2=-c_4 n,\qquad \delta_0=0,
$$

i.e. a Pigouvian tax equal to the marginal external adjustment cost $c_4(Y_t-Y_{t-1})$. (Any constant
$\delta_0$ merely shifts the level; $\delta_0=0$ keeps the tax revenue-neutral on average.) With these
$\delta$'s the competitive Euler equation coincides with the planner's, so the competitive equilibrium
solves the planning problem.
```

**7.** *(Duck decoys.)* A fixed number $n$ of identical firms produce duck decoys, $y_t=fk_t$ ($f>0$),
aggregate $Y_t=fK_t$, $K_t=nk_t$. Demand $p_t=A_0-A_1 Y_t+u_t$ ($A_0,A_1>0$) with
$u_t=\frac{1}{1-\rho L}\varepsilon_{ut}$, $|\rho|<1$. The rental on capital obeys the upward-sloping
industry supply $w_t=B_0+B_1 K_t$ ($B_0,B_1>0$). Net cash flow
$\pi_t=p_t y_t-w_t k_t-\tfrac{d}{2}(k_t-k_{t-1})^2$, present value $E\sum_{t=0}^{\infty}b^t\pi_t$.

&nbsp;&nbsp;**A.** Define a rational expectations competitive equilibrium.
&nbsp;&nbsp;**B.** Compute it.
&nbsp;&nbsp;**C.** What social planning problem does it solve?
&nbsp;&nbsp;**D.** Define a monopolistic (collusive) rational expectations equilibrium.
&nbsp;&nbsp;**E.** Compute it.
&nbsp;&nbsp;**F.** Show that the feedback coefficient of $K_t$ on $K_{t-1}$ is *smaller* under monopoly
than under competition.
&nbsp;&nbsp;**G.** Justice Department economists know $f,n,b,d$ but not $(B_0,B_1,A_0,A_1)$; they observe
$\{(Y_s,K_s)\}$. Can they tell whether the industry is competitive or monopolistic?
&nbsp;&nbsp;**H.** Same, but they also observe $\{(p_s,w_s)\}$. Can they now tell?

```{admonition} Solution to Exercise 7
:class: dropdown

**A.** A process $\{K_t\}$ (and $\{p_t\}$) such that each firm, price-taking on $p_t$ and on the
rental $w_t=B_0+B_1 K_t$ (using a perceived law of motion for $K_t$), maximizes present value, and the
markets clear: $p_t=A_0-A_1 fK_t+u_t$, $w_t=B_0+B_1 K_t$.

**B.** The competitive Euler equation $E_t\{p_t f-w_t-d(k_t-k_{t-1})+bd(k_{t+1}-k_t)\}=0$, with the
demand and rental curves substituted *after* differentiating and multiplied by $n$, is

$$
-bE_tK_{t+1}+\big[(1+b)+d^{-1}n(A_1 f^2+B_1)\big]K_t-K_{t-1}=d^{-1}[n(A_0 f-B_0)+nf u_t].
$$

Factor $(1-\lambda_1^c L)(1-\lambda_2^c L)$ with $\lambda_1^c<1<\lambda_2^c$ and solve the stable root
backward using the AR(1) forecast of $u$: $K_t=\lambda_1^c K_{t-1}+(\text{terms in }u_t)$.

**C.** With no externality, the competitive equilibrium solves the planner's surplus problem: maximize
$E\sum b^t\{\int_0^{Y_t}(A_0-A_1x+u_t)dx-\int_0^{K_t}(B_0+B_1K)dK-\tfrac{d}{2}\sum_{\text{firms}}(k_t-k_{t-1})^2\}$
(area under demand minus area under the capital-supply curve minus adjustment costs).

**D.** A monopolistic (collusive) equilibrium: the $n$ firms jointly choose $\{K_t\}$ to maximize
industry value, *internalizing* the effect of $K_t$ on price — i.e. substitute $p_t=A_0-A_1 fK_t+u_t$
into the objective **before** differentiating.

**E.** The monopoly marginal revenue is $A_0 f-2A_1 f^2 K_t+fu_t$ (the extra $-A_1 f^2 K_t$), giving

$$
-bE_tK_{t+1}+\big[(1+b)+d^{-1}n(2A_1 f^2+B_1)\big]K_t-K_{t-1}=d^{-1}[n(A_0 f-B_0)+nf u_t].
$$

**F.** In both cases the middle coefficient is $-\phi=(1+b)+\theta d^{-1}$ with
$\theta_{\text{comp}}=n(A_1 f^2+B_1)$ and $\theta_{\text{mon}}=n(2A_1 f^2+B_1)>\theta_{\text{comp}}$.
The stable root solves $-\phi=b\lambda_1+\lambda_1^{-1}$; on the branch $\lambda_1\in(0,1)$ the function
$b\lambda+\lambda^{-1}$ is decreasing, so a *larger* $-\phi$ gives a *smaller* $\lambda_1$ (the
Figure-4 argument of {doc}`Chapter IX <ch09_difference_equations>`). Hence
$\lambda_1^{\text{mon}}<\lambda_1^{\text{comp}}$: monopoly adjusts capital with less persistence.

**G.** **No.** From $\{(Y_s,K_s)\}$ one recovers only the reduced-form feedback $\lambda_1$, which
depends on the *composite* $\theta=n(A_1 f^2+B_1)$ (competition) or $n(2A_1 f^2+B_1)$ (monopoly). Since
$(A_1,B_1)$ are unknown, any observed $\lambda_1$ can be rationalized by a competitive model with one
$(A_1,B_1)$ or a monopoly model with another — the market structure is not identified.

**H.** **Yes.** Observing $\{(p_s,w_s)\}$ lets one estimate the demand slope $A_1$ (regress $p$ on $Y$)
and the capital-supply slope $B_1$ (regress $w$ on $K$). With $A_1,B_1$ known and $f,n,b,d$ given, one
computes the predicted $\lambda_1$ under competition versus monopoly ($\theta$ versus $\theta+nA_1 f^2$)
and compares with the observed feedback in the $K_t$ law of motion. Because the two predictions differ
(monopoly's $\lambda_1$ is smaller), the cross-equation restriction linking the demand/supply slopes to
the equilibrium $\lambda_1$ **identifies** whether the industry is competitive or monopolistic.
```

## References

- Eckstein, Z. (1984). A rational expectations model of agricultural supply. *Journal of Political Economy* 92(1), 1–19.
- Eichenbaum, M. S. (1983). A rational expectations equilibrium model of inventories of finished goods and employment. *Journal of Monetary Economics* 12(2), 259–277.
- Gordon, D. F., and A. Hynes (1970). On the theory of price dynamics. In E. S. Phelps et al. (eds.), *Microeconomic Foundations of Employment and Inflation Theory*. Norton.
- Hansen, L. P., and T. J. Sargent (1980). Formulating and estimating dynamic linear rational expectations models. *Journal of Economic Dynamics and Control* 2(1), 7–46.
- Hansen, L. P., and T. J. Sargent (1981). Linear rational expectations models for dynamically interrelated variables. In R. E. Lucas Jr. and T. J. Sargent (eds.), *Rational Expectations and Econometric Practice*. University of Minnesota Press.
- Holt, C. C., et al. (1960). *Planning Production, Inventories, and Work Force*. Prentice-Hall.
- Kwakernaak, H., and R. Sivan (1972). *Linear Optimal Control Systems*. Wiley.
- Lucas, R. E. Jr. (1976). Econometric policy evaluation: a critique. In K. Brunner and A. Meltzer (eds.), *The Phillips Curve and Labor Markets*, Carnegie–Rochester Conference Series 1. North-Holland.
- Lucas, R. E. Jr., and E. C. Prescott (1971). Investment under uncertainty. *Econometrica* 39(5), 659–681.
- Romer, P. M. (1983). Externalities and increasing returns in dynamic competitive analysis. Working Paper, University of Rochester.
- Sargent, T. J. (1981). Interpreting economic time series. *Journal of Political Economy* 89(2), 213–248.
- Taylor, J. B. (1979). Estimation and control of a macroeconomic model with rational expectations. *Econometrica* 47(5), 1267–1286.
- Telser, L. G., and R. L. Graves (1972). *Functional Analysis in Mathematical Economics*. University of Chicago Press.
- Townsend, R. (1983). Forecasting the forecasts of others. *Journal of Political Economy* 91(4), 546–588.

[^fn-14-1]: The approach originally stems from Holt et al. (1960). For a rigorous treatment see
Telser and Graves (1972). Engineers solve these systems as matrix Riccati equations; see Kwakernaak
and Sivan (1972). In their jargon the systems here are *not* "controllable" but are "stabilizable"
and "detectable," so convergence of the Riccati equations in the infinite-horizon problem is
assured.

[^fn-14-2]: It is necessary to distinguish two operators $B$ and $L$. The operator $B$ is defined
by $B^{-1}[Ex_{t+j}\mid\Omega_{t-1}]=Ex_{t+j+1}\mid\Omega_{t-1}$ — application of $B^{-1}$ shifts
forward by one period the date on the variable whose conditional forecast is being computed, but
leaves the information set unaltered. The lag operator $L$ is defined by $L^j x_t=x_{t-j}$; in
particular $L^{-1}(Ex_{t+j}\mid\Omega_{t-1})=Ex_{t+j+1}\mid\Omega_t$, so $L^{-1}$ shifts *both* the
random variable and the information set forward one period.

[^fn-14-3]: The properties of $B$ make the forward inverse of $1-\lambda_2 B$ the only legitimate
one (apart from convergence). Operating with polynomials in nonpositive powers of $B$ is
legitimate, but operating with polynomials in *positive* powers of $B$ is not: for example
$E_t x_{t+1}=E_t y_{t+1}$ does not imply $BE_t x_{t+1}=BE_t y_{t+1}$, i.e. $x_t=y_t$. The operation
in the text involves only nonpositive powers of $B$.

[^fn-14-4]: Another stochastic process for $\{n_{t+j}\}$ that satisfies the Euler equations and the
transversality condition is $n_{t+j}=\lambda_1 n_{t+j-1}-\lambda_1\sum_{i=0}^\infty\lambda_2^{-i}z_{t+j+i}$
(with the *actual* future $z$'s, not their forecasts). But this depends on future values not
observable at $t+j$, so it is not a solution of our problem, which requires $n_{t+j}$ to be a
function of information available at $t+j$. The solution that is also a function of that information
is said to be *realizable* or *nonanticipative*.

[^fn-14-5]: Because of the quadratic objective and linear law of motion, the decision rule
{eq}`eq-14-17` is independent of the variance of the innovation $\epsilon_t$; the points here hold
even if $E\epsilon_t^2=0$, in which case the firm forecasts future $z$'s perfectly from current and
past $z$'s.

[^fn-14-6]: The transversality condition is derived by a procedure analogous to the one used in
Section 1 and in {doc}`Chapter IX <ch09_difference_equations>`.

[^fn-14-7]: To obtain a rational expectations *competitive* equilibrium, $A_0-A_1 f_0 nk_{t+j}+u_{t+j}$
is substituted for $p_{t+j}$ only *after* the Euler equation has been obtained — i.e. after
{eq}`eq-14-20` has been differentiated with respect to $k_{t+j}$ — which ensures the firm acts as a
price taker. Substituting *before* differentiating yields a rational expectations *monopoly*
equilibrium.

[^fn-14-def]: Definition $w_t=J_t-bE_t J_{t+1}$ is chosen so that replacing $J_{t+j}(k_{t+j}-k_{t+j-1})$
with $w_{t+j}k_{t+j}$ in {eq}`eq-14-27` leads to the same marginal condition {eq}`eq-14-29` (a
summation-by-parts argument).
