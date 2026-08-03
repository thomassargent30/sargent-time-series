# Lucas's Two Illustrations of the Quantity Theory and Whiteman's Critique

```{note}
This section describes Robert E. Lucas, Jr., "Two Illustrations of the Quantity Theory of Money,"
*American Economic Review* 70(5), 1005–1014 (1980); Charles H. Whiteman, "Lucas on the Quantity Theory:
Hypothesis Testing without Theory," *American Economic Review* 74(4), 742–749 (1984); and Thomas J.
Sargent and Paolo Surico, "Two Illustrations of the Quantity Theory of Money: Breakdowns and Revivals,"
*American Economic Review* 101(1), 109–128 (2011). It is the third application in this book of
**Sims's approximation-error formula** — after {doc}`33a_seasonality_approximation` and
{doc}`36c_cagan_hyperinflation` — and it is perhaps the cleanest illustration of what that formula is
good for. A computational treatment, with Python code for the filters, the scatter plots, the VAR, and
the Bayesian estimation of the DSGE model, is the QuantEcon lecture
[*Two Illustrations of the Quantity Theory of Money*](https://python.quantecon.org/sargent_surico.html).
```

## Lucas's method

{cite:t}`lucas1980two` set out to illustrate two central quantity-theoretic propositions: a change in the growth
rate of money induces (i) an equal change in the rate of inflation, and (ii) an equal change in the
nominal interest rate. Rather than estimate a structural model — a course he judged to require
"nesting the two hypotheses in question within a complex maintained hypothesis, which must be accepted
as valid in order to carry out the test" — he applied a deliberately *atheoretical* two-step procedure.

**Step 1: filter.** Replace each series by a two-sided, exponentially weighted moving average. For money
growth $\mu_t$,

```{math}
:label: eq-lw-filter
\bar\mu_t(\beta) \;=\; \frac{1-\beta}{1+\beta}\sum_{k=-\infty}^{\infty}\beta^{|k|}\mu_{t-k},
\qquad 0<\beta<1,
```

the weights being normalized to sum to one. **Step 2: scatter.** Plot the filtered inflation–money
pairs $(\bar\mu_t(\beta),\bar\pi_t(\beta))$ and the filtered interest–money pairs
$(\bar\mu_t(\beta),\bar\rho_t(\beta))$ for a range of $\beta$, and look at the slope. For US data over
1955–1975 and $\beta$ near $1$, both scatters hugged the $45^\circ$ line: Lucas read
$a_1\approx a_2\approx 1$ in

```{math}
:label: eq-lw-scatter
\bar\pi_t(\beta) = a_1\,\bar\mu_t(\beta)+\varepsilon_{1t},
\qquad
\bar\rho_t(\beta) = a_2\,\bar\mu_t(\beta)+\varepsilon_{2t},
```

as evidence for the two quantity-theory propositions, and hence *against* the **Mundell–Tobin effect**,
by which higher inflation lowers the return on money, shifts portfolios toward real capital, and makes
the nominal interest rate rise *less* than one-for-one with money growth.

Using the {doc}`filter kit <10_filter_kit>` of this chapter, the transfer function of
{eq}`eq-lw-filter` is

```{math}
:label: eq-lw-gain
B(e^{-i\omega}) \;=\; \frac{(1-\beta)^2}{1+\beta^2-2\beta\cos\omega}.
```

Note that $B(1)=1$ for every $\beta$ — the filter has unit gain at frequency zero — while for any fixed
$\omega\ne0$, $B(e^{-i\omega})\to0$ as $\beta\to1$. So {eq}`eq-lw-filter` is a **low-pass** filter that,
as $\beta\uparrow1$, concentrates all of its power at $\omega=0$. Filtered money growth has spectral
density $\lvert B(e^{-i\omega})\rvert^2 S_{\mu\mu}(e^{-i\omega})$, which for $\beta$ near one is
massed almost entirely at the origin.

## Whiteman's first point: what the method measures

Whiteman's (1984) first contribution was to say precisely what Lucas's graphical procedure *estimates*.
Consider the **two-sided** population projection of the nominal interest rate on the entire money-growth
process,

```{math}
:label: eq-lw-proj
\rho_t = \sum_{k=-\infty}^{\infty}\gamma_k\,\mu_{t-k} + \eta_t,
\qquad E\mu_{t-k}\eta_t = 0 \ \ \forall k,
```

whose transfer function is the familiar ratio of the cross spectrum to the spectrum
({doc}`07_cross_spectrum`, {doc}`30_filtering_projections`),

```{math}
:label: eq-lw-transfer
\gamma(e^{-i\omega}) = \frac{S_{\rho\mu}(e^{-i\omega})}{S_{\mu\mu}(e^{-i\omega})} .
```

Equation {eq}`eq-lw-proj` cannot be estimated — it has infinitely many parameters — so any practitioner
fits some *restricted* lag distribution $\gamma'$. By **Sims's approximation-error formula**
({ref}`Exercise 1 <ex-1>` of this chapter), least squares in population chooses $\gamma'$ to minimize
the spectral-density-weighted distance

```{math}
:label: eq-lw-sims
\int_{-\pi}^{\pi}\big|\gamma(e^{-i\omega})-\gamma'(e^{-i\omega})\big|^{2}\,
S_{\mu\mu}(e^{-i\omega})\,d\omega .
```

Now apply {eq}`eq-lw-sims` to Lucas's procedure. Regressing one filtered series on another with a
*single* coefficient $a_2$ is the extreme restriction $\gamma'(e^{-i\omega})\equiv a_2$, a **constant**
lag-distribution transfer function; and the relevant weighting density is that of *filtered* money
growth. So Lucas's method chooses $a_2$ to minimize

```{math}
:label: eq-lw-lucascrit
\int_{-\pi}^{\pi}\big|\gamma(e^{-i\omega})-a_2\big|^{2}\,
S_{\bar\mu(\beta)}(e^{-i\omega})\,d\omega .
```

Because the weighting density is massed at $\omega=0$ when $\beta$ is near one, essentially all of the
weight falls on getting the approximation right at the origin. Hence

```{math}
:label: eq-lw-h0
a_2 \;\longrightarrow\; \gamma(e^{-i\cdot 0}) \;=\; \gamma(1) \;=\; \sum_{k=-\infty}^{\infty}\gamma_k
\;=\; \frac{S_{\rho\mu}(0)}{S_{\mu\mu}(0)} \qquad\text{as }\beta\to1 .
```

**Lucas's scatter-plot slope is an estimate of the sum of the coefficients in a two-sided distributed
lag regression** — equivalently, the ratio of the cross spectrum to the spectrum at frequency zero.
Sargent and Surico write this object $h(0)$; it is a *population* magnitude that any time series model
delivers. Whiteman emphasized that this holds "regardless of the time-series properties" of the two
series: the method is a robust, if extremely roundabout, estimator of a sum of lag coefficients.

Two features of the book's machinery deserve notice here. First, it is legitimate to filter *both*
series with the same filter precisely because, as {doc}`30_filtering_projections` establishes, applying
a common filter to $y$ and $x$ **leaves the two-sided projection unaltered** (it is the *one-sided*
projection that a common filter disturbs, except when Granger non-causality makes the two coincide).
Lucas's low-pass filter therefore reweights *which frequencies the fit attends to* without changing the
$\gamma_k$'s being estimated. Second, the exercise is a benign cousin of the
{doc}`Slutsky–Kuznets effects <09_slutsky_kuznets>`: filtering can manufacture apparent regularities,
and one must know exactly what a filter does before reading economics off a smoothed plot.

```{admonition} Why not just add up estimated lag coefficients?
:class: dropdown

Whiteman's footnote on this point is a small gem of projection theory, and it explains why Lucas's
roundabout method is actually the *safer* one.

Suppose that instead of Lucas's procedure you fit a parsimonious rational lag distribution
$\gamma'(e^{-i\omega})$ by least squares and then report $\gamma'(1)$. The trouble is that the sum of
coefficients is **not a continuous function of the least-squares criterion** {eq}`eq-lw-sims`. The
criterion is an integral, and the single frequency $\omega=0$ is a set of measure zero: you can change
$\gamma'(1)$ *arbitrarily* — even send it to infinity — while barely moving the value of the integral,
i.e. without materially changing the fit.

Whiteman illustrated this with Mills's (1982) estimated transfer function. Perturbing one estimated
coefficient from $-0.56$ to $-0.58$ — a change amounting to less than a tenth of one standard error,
and invisible in the fit — drives the denominator of the rational lag distribution toward zero at
$z=1$, so that $\gamma'(1)$ swings from unity to $+\infty$. The reported standard error on the sum is
therefore badly misleading.

Lucas's method escapes this trap because it estimates the sum by fitting a *constant*, and a constant
transfer function **is** continuous with respect to the least-squares metric. The price is that
Lucas's procedure deliberately sacrifices accuracy about individual lag coefficients in order to buy
accuracy about their sum — and, as Whiteman notes, the same "controlled" approximation error that makes
the sum well estimated destroys any hope of attaching a standard error to it.
```

## Whiteman's second point: what the measurement means

Knowing that Lucas measured $\gamma(1)$ still leaves the substantive question: does $\gamma(1)\approx1$
constitute evidence against the Mundell–Tobin effect? To answer, Whiteman did what Lucas had declined to
do — he adopted an explicit structural model and computed the *population* version of Lucas's statistic
inside it. Tellingly, he chose **Lucas's own (1975) equilibrium model of the business cycle**, which is
convenient because it isolates the Mundell–Tobin effect in a single parameter.

In that model, log output and the log real return to capital are $y_t=\delta_0'+\delta_1'k_t$ and
$r_t=\delta_0-\delta_1k_t$; the demands for capital and for real balances are

```{math}
:label: eq-lw-model
\begin{aligned}
k_{t+1} &= \alpha_0' + \alpha_1 E_t r_{t+1} + \alpha_2\,(E_t p_{t+1}-p_t) + \alpha_3 k_t,\\
m_t - p_t &= \beta_0' - \beta_1 E_t r_{t+1} - \beta_2\,(E_t p_{t+1}-p_t) + \beta_3 k_t,
\end{aligned}
```

with money supply $m_t = A(L)e_t$ for a fundamental $\{e_t\}$. **The parameter $\alpha_2$ is the
Mundell–Tobin effect**: when $\alpha_2=0$ the demand for capital is independent of expected inflation,
monetary disturbances leave the steady-state capital stock alone, and the second equation of
{eq}`eq-lw-model` reduces to a version of *Cagan's portfolio-balance schedule* — the very schedule
studied in {doc}`36c_cagan_hyperinflation` and {doc}`22_rational_expectations`, whose forward solution
$p_t=\frac{1}{1+\beta_2}\sum_{k\ge0}\big(\frac{\beta_2}{1+\beta_2}\big)^k E_t m_{t+k}$ has lead weights
summing to one. In that case the model does exhibit the two quantity-theory propositions.

Two pieces of the book's apparatus organize the solution. Whiteman assumed that $(k_{t+1},p_t)$ fails to
Granger-cause $m_t$; by {doc}`Sims's theorem <27_granger_causality>` this makes the moving average
representation **triangular**,

```{math}
:label: eq-lw-triangular
\begin{bmatrix} m_t \\ k_{t+1} \\ p_t\end{bmatrix}
=
\begin{bmatrix} A(L) & 0 & 0\\ B(L) & D(L) & G(L)\\ C(L) & F(L) & H(L)\end{bmatrix}
\begin{bmatrix} e_t \\ w_t \\ v_t\end{bmatrix},
```

exactly the device used in {doc}`36c_cagan_hyperinflation` and {doc}`36a_interpreting_vars`. Rational
expectations then determine $B(L)$ and $C(L)$ uniquely from the money-supply rule $A(L)$, by the
familiar factorization into stable and unstable roots and the
{doc}`geometric-lead <20_geometric_leads>` formulas of this chapter. The transfer function of the
inflation-on-money-growth projection turns out to be simply $\gamma^{\pi}(z)=C(z)/A(z)$.

Whiteman's conclusions are these.

- **The first illustration is uninformative about $\alpha_2$.** When $\alpha_2=0$,
  $\gamma^\pi(1)=1-\frac{\beta_2}{1+\beta_2}A\big(\frac{\beta_2}{1+\beta_2}\big)A(1)^{-1}$, so
  $\gamma^\pi(1)=1$ requires $A(1)^{-1}=0$, i.e. a *random walk component in the money supply* — the
  condition Lucas had himself employed in his 1972 and 1975 papers. But when that condition holds,
  $\gamma^\pi(1)=1$ *whatever* the value of $\alpha_2$. Lucas's first result is essentially independent
  of the magnitude of the Mundell–Tobin effect.
- **The second illustration points the other way.** Conditioning on $\gamma^\pi(1)=1$, the requirement
  $\gamma^\rho(1)=1$ becomes a restriction involving *all* the model's parameters together with those
  of the money-supply rule; it implies $\alpha_2=0$ only if $A^*(\lambda_2^{-1})/A^*(1)=0$, which needs
  either a random walk in *money growth* or the absence of an autoregressive representation. Neither is
  plausible. Moreover Lucas's own figures argue against the first: were the money spectrum already
  infinite at $\omega=0$, raising $\beta$ (which enhances low-frequency power) would leave $a_1$ and
  $a_2$ unchanged — yet Whiteman's replication found $a_1$ moving from $0.02$ to $0.08$ to $0.87$ to
  $0.99$ as $\beta$ ran through $0,\,0.5,\,0.9,\,0.95$. The most plausible reading of the evidence taken
  as a whole is that the restriction holds *with $\alpha_2\ne0$* — that is, **Lucas's numbers can be
  read as evidence in favor of the very Mundell–Tobin effect he took them to refute.**

The moral is a **Lucas-critique-based admonition** turned on Lucas's own paper. Reduced-form
low-frequency correlations are complicated functions of the structure of the economy *and* of the laws
of motion of the processes agents care about; they are unlikely to reveal much about the true model
"unless one already has it very much in mind." This is the same non-invariance that drives
{doc}`22_rational_expectations`, {doc}`36a_interpreting_vars`, {doc}`36b_exact_linear_re`,
{doc}`Chapter XIV <ch14_investment_uncertainty>` and — in the guise of a misspecified regression whose
weighting spectrum is generated by the policy rule — {doc}`36c_cagan_hyperinflation`. Whiteman's title
says it: *beware of hypothesis testing without theory.*

## Sargent and Surico: breakdowns and revivals

{cite:t}`sargentsurico2011two` take up Whiteman's reformulation and push it in a quantitative direction,
asking when Lucas's unit slopes "obtain in reality" and when they "break down" — a caveat Lucas himself
had attached to his propositions. Their organizing object is exactly Whiteman's,

```{math}
:label: eq-lw-hzero
h(0) \;=\; \frac{S_{yz}(0)}{S_{z}(0)},
```

which any state-space model delivers. If $X_{t+1}=AX_t+BW_{t+1}$, $Y_t=CX_t$ with $A$ stable, then
writing the long-run multiplier $G=C(I-A)^{-1}B$ gives $2\pi S_Y(0)=GG'$ and hence
$h_{y,z}(0)=[GG']_{yz}/[GG']_{zz}$ — a one-line computation from the
{doc}`compact state-space notation <24_compact_notation>` of this chapter.

Their findings:

- **The slopes are unstable.** Extending the data to 1900–2005 and computing Lucas's regressions on
  filtered data, the inflation-on-money-growth slope is $1.13$ for 1900–1928, $0.39$ for 1929–1954,
  $0.86$ for Lucas's own 1955–1975, and $-0.03$ for 1984–2005; the interest-rate slope is $0.62$ for
  1955–1975 but $0.06$ or negative in several other subperiods. Lucas's sample turns out to be
  *atypical*. A VAR with time-varying coefficients and stochastic volatility tells the same story: the
  posterior median of $h(0)$ drifts substantially, peaking in the 1970s and reaching its lowest values
  in the 1940s and in the most recent twenty years.
- **Monetary policy accounts for the instability.** Estimating a New Keynesian DSGE model over
  1960–1983 reproduces Lucas's unit slopes (implied $h_{\pi,\Delta m}(0)\approx1.01$,
  $h_{R,\Delta m}(0)\approx0.82$). Then, *freezing every nonpolicy parameter* and varying only the
  policy-rule coefficients, the implied $h(0)$'s sweep across essentially the whole range from $0$ to
  $1$ — covering the entire span of the empirical estimates. **The low-frequency slopes are not policy
  invariant.**
- **The economics of the breakdown.** A policy rule that responds aggressively to incipient inflation
  prevents persistent movements in money growth from emerging and thereby *eradicates* the unit slopes;
  a rule that responds too weakly, acceding to persistent money growth, *revives* them. As the authors
  put it, Lucas's illustrations "come back" precisely when a monetary authority allows persistent
  movements in money growth — so, as citizens, they prefer the times when the propositions break down.

This is Whiteman's argument made operational and quantitative: the same statistic, computed inside a
fully articulated model, moves with the money-supply rule. The frequency-zero slope measures a genuine
feature of the data, but *what it means* depends on the policy regime that generated the money process
— which is, once more, the weighting-spectrum lesson of Sims's approximation formula and the
substantive lesson of the Lucas critique.

```{admonition} Exercise
:class: tip

Let $\mu_t$ have spectral density $S_{\mu\mu}$, and let $\bar\mu_t(\beta)$ be the Lucas filter
{eq}`eq-lw-filter` with transfer function {eq}`eq-lw-gain`.

**A.** Verify that $B(1)=1$ for all $\beta\in(0,1)$, and that $B(e^{-i\omega})\to0$ as $\beta\to1$ for
each fixed $\omega\ne0$.

**B.** Suppose $y_t=\gamma(L)\mu_t+\eta_t$ as in {eq}`eq-lw-proj`. Show that the population regression
coefficient of $\bar y_t(\beta)$ on $\bar\mu_t(\beta)$ is a weighted average of
$\gamma(e^{-i\omega})$, and identify the weights. Conclude that the coefficient $\to\gamma(1)$ as
$\beta\to1$.
```

```{admonition} Solution
:class: dropdown

**A.** At $\omega=0$, $\cos\omega=1$, so the denominator of {eq}`eq-lw-gain` is
$1+\beta^2-2\beta=(1-\beta)^2$, which equals the numerator; hence $B(1)=1$ for every $\beta$. For fixed
$\omega\ne0$ the denominator tends to $2-2\cos\omega>0$ as $\beta\to1$ while the numerator $(1-\beta)^2$
tends to $0$, so $B(e^{-i\omega})\to0$. The filter's power is thus squeezed onto an ever-shrinking
neighborhood of the origin while its gain there stays pinned at one — a low-pass filter converging to a
"delta" at $\omega=0$.

**B.** Filtering both series with the common filter $B$ gives $\bar y_t=\gamma(L)\bar\mu_t+\bar\eta_t$
with $\bar\eta_t$ still orthogonal to the whole $\bar\mu$ process (this is the invariance of the
two-sided projection under a common filter, {doc}`30_filtering_projections`). The population regression
coefficient of $\bar y$ on $\bar\mu$ is therefore

$$
b(\beta)=\frac{E\,\bar y_t\bar\mu_t}{E\,\bar\mu_t^2}
=\frac{\int_{-\pi}^{\pi}\gamma(e^{-i\omega})\,|B(e^{-i\omega})|^2 S_{\mu\mu}(e^{-i\omega})\,d\omega}
       {\int_{-\pi}^{\pi}|B(e^{-i\omega})|^2 S_{\mu\mu}(e^{-i\omega})\,d\omega},
$$

a weighted average of the transfer function $\gamma(e^{-i\omega})$ with weights proportional to the
spectral density of *filtered* money growth, $|B(e^{-i\omega})|^2S_{\mu\mu}(e^{-i\omega})$. (Equivalently,
$b(\beta)$ minimizes {eq}`eq-lw-lucascrit` over constants — differentiate and solve.) By part A those
weights concentrate at $\omega=0$ as $\beta\to1$, so provided $S_{\mu\mu}$ is continuous and positive at
the origin, $b(\beta)\to\gamma(1)=\sum_k\gamma_k=S_{y\mu}(0)/S_{\mu\mu}(0)$. $\blacksquare$
```

## References

```{bibliography}
:labelprefix: LW
:filter: key in {"boschenotrok1994longrun", "cagan1956monetary", "ireland2004technology", "lucas1972expectations", "lucas1975equilibrium", "lucas1976econometric", "lucas1980two", "mccallum1984low", "mills1982signal", "mundell1963inflation", "sargentsurico2011two", "sims1972approx", "sims1972money", "tobin1965money", "vogel1974dynamics", "whiteman1983linear", "whiteman1984lucas"}
```
