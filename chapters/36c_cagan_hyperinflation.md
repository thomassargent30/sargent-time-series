# Money Demand in Hyperinflations: A Misspecified Regression and Sims's Approximation Formula

```{note}
This section is based on Thomas J. Sargent, "The Demand for Money During Hyperinflations under
Rational Expectations: I," *International Economic Review* 18(1), 59–82 (1977). It is a capstone
application that ties together several threads of this chapter: the
{doc}`rational-expectations Cagan model <22_rational_expectations>`, the
{doc}`Granger-causality / econometric-exogeneity theorem <27_granger_causality>` and its
{doc}`money–income application <28_sims_money_income>`, the cross-spectral notions of *phase* and
*leading* of {doc}`07_cross_spectrum` and {doc}`08_leading_indicators`, and — the organizing idea —
**Sims's frequency-domain approximation-error formula**, stated as
[Exercise 1](37_exercises.md) of this chapter and deployed in {doc}`33a_seasonality_approximation`.
The full computational treatment, with Python code for the simulations, the maximum-likelihood
estimator, and the empirical tables, is the QuantEcon lecture
[*Demand for Money during Hyperinflations under Rational Expectations*](https://python-advanced.quantecon.org/cagan_rational_expectations.html);
here we summarize the quantitative findings and dwell on the projection theory.
```

## Cagan's paradox

Cagan's (1956) classic study of seven hyperinflations fit a demand schedule for real balances,

```{math}
:label: eq-cg-demand
m_t - p_t = \alpha\,\pi_t + u_t, \qquad \alpha < 0,
```

where $m_t$ is log money, $p_t$ log price level, $\pi_t$ the public's expected inflation, and $u_t$ a
mean-zero disturbance. Cagan estimated $\alpha$ — the semi-elasticity of real balances with respect to
expected inflation — and used it to compute the inflation rate $-1/\alpha$ that would *maximize* the
inflation-tax revenue a money-printing government can extract. For every one of the seven
hyperinflations, the reciprocal of Cagan's estimate of $-\alpha$ came out **far below** the average
actual inflation rate: money creators appeared to be inflating at rates wildly in excess of the
revenue-maximizing rate. The natural suspicion is that this paradox is a *statistical artifact* — a
consequence of a biased estimate of $\alpha$. This section shows that it is, and that the bias is an
instance of Sims's approximation-error formula applied to a misspecified regression.

## Cagan's model under rational expectations

Cagan assumed adaptive expectations, $\pi_t=\dfrac{1-\lambda}{1-\lambda L}\,x_t$, where $x_t\equiv p_t-p_{t-1}$
is inflation and $L$ the lag operator. Sargent and Wallace (1973) observed that
under **rational expectations**, $\pi_t = E_t x_{t+1}$, and solving the resulting forward difference
equation with $\lvert{-\alpha}/(1-\alpha)\rvert<1$ gives

```{math}
:label: eq-cg-re
\pi_t = E_t x_{t+1}
= \frac{1}{1-\alpha}\sum_{j=1}^{\infty}\Big(\tfrac{-\alpha}{1-\alpha}\Big)^{j-1}E_t\mu_{t+j}
  - \frac{1}{1-\alpha}\sum_{j=1}^{\infty}\Big(\tfrac{-\alpha}{1-\alpha}\Big)^{j-1}\big(E_t u_{t+j}-E_t u_{t+j-1}\big),
```

with $\mu_t\equiv m_t-m_{t-1}$ the money-growth rate. Equation {eq}`eq-cg-re` is a
{doc}`geometric distributed lead <20_geometric_leads>` of the same kind studied throughout this
chapter: expected inflation is a discounted sum of *expected future money growth*, so — with the
stochastic process for money creation held fixed — "money causes inflation." Adaptive expectations
coincide with rational expectations only under restrictions on $u$ and $\mu$; Sargent studies two
sufficient conditions: $u_t$ is a random walk, $u_t=u_{t-1}+\eta_t$ (so $E_t(u_{t+j}-u_{t+j-1})=0$),
and the money-creation rule

```{math}
:label: eq-cg-rule
\mu_t = \frac{1-\lambda}{1-\lambda L}\,x_t + \varepsilon_t \;\;(=E_t x_{t+1}+\varepsilon_t),
```

with $\varepsilon_t,\eta_t$ serially uncorrelated, mean zero. Under {eq}`eq-cg-rule` the geometric sum
in {eq}`eq-cg-re` collapses to $\pi_t=\frac{1-\lambda}{1-\lambda L}x_t$, so adaptive expectations *are*
rational. Rule {eq}`eq-cg-rule` describes a government that prints money in response to inflation — a
"real-bills" regime of the sort German officials repeatedly invoked to argue that money was responding
to inflation rather than causing it.

Under these conditions inflation and money growth form a bivariate system that reduces, after
first-differencing, to a first-order moving average driven by $(\varepsilon_t,\eta_t)$. Writing
$\phi\equiv(\lambda+\alpha(1-\lambda))^{-1}$,

```{math}
:label: eq-cg-ma
\begin{aligned}
(1-L)\,x_t &= \phi\,(1-\lambda L)\,(\varepsilon_t-\eta_t),\\
(1-L)\,\mu_t &= \phi(1-\lambda)(\varepsilon_t-\eta_t) - \varepsilon_{t-1}+\varepsilon_t .
\end{aligned}
```

## The bivariate Wold representation and Granger causality

To find what Cagan's regression converges to we need the *fundamental* (Wold) representation of
$(\Delta x_t,\Delta\mu_t)$. Project the money-supply shock on the inflation shock: write
$\varepsilon_t=\rho(\varepsilon_t-\eta_t)+v_t$ with

```{math}
:label: eq-cg-rho
\rho = \frac{E[\varepsilon_t(\varepsilon_t-\eta_t)]}{E[(\varepsilon_t-\eta_t)^2]}
     = \frac{\sigma_\varepsilon^2-\sigma_{\varepsilon\eta}}{\sigma_\varepsilon^2-2\sigma_{\varepsilon\eta}+\sigma_\eta^2},
```

so that $v_t\perp(\varepsilon_t-\eta_t)$. Substituting into {eq}`eq-cg-ma` gives the **triangular**
bivariate Wold representation

```{math}
:label: eq-cg-wold
\begin{aligned}
(1-L)\,x_t &= \phi(1-\lambda L)\,(\varepsilon_t-\eta_t),\\
(1-L)\,\mu_t &= \big[\phi(1-\lambda)+\rho(1-L)\big]\,(\varepsilon_t-\eta_t) + (1-L)\,v_t,
\end{aligned}
```

with jointly fundamental white noises $(\varepsilon_t-\eta_t,\,v_t)$. The **lower-triangular** structure
— $\Delta x$ loads only on the first shock, $\Delta\mu$ on both — is exactly the condition of
{doc}`Sims's theorem <27_granger_causality>`: $\Delta x$ is *econometrically exogenous* with respect to
$\Delta\mu$, so **inflation Granger-causes money creation but not conversely**. This is the same
one-sidedness that organizes Sims's {doc}`money–income test <28_sims_money_income>`, here running from
prices to money rather than money to income.

A caution from {doc}`08_leading_indicators` applies with full force. That $x$ *causes* $\mu$ in
Granger's sense does **not** mean $x$ *leads* $\mu$ in any National-Bureau sense. Indeed, the
equilibrium implies $\mu_t = x_t + a_{2t}-a_{1t}$ (with $a_{1t},a_{2t}$ the innovations in $x,\mu$), so
$x_t$ and $\mu_t$ are *in phase* — their cross-spectrum has **zero phase at every frequency**
({doc}`07_cross_spectrum`). Evidence that inflation *leads* money would be evidence *against* the
model. Granger causality and phase leads are different things, exactly as {doc}`08_leading_indicators`
warned.

Finally, the Granger-causal pattern here is a property of the *particular money rule* {eq}`eq-cg-rule`,
not an invariant feature of the economy: change the money-supply regime and the causal ordering can
change. This is the distinction — central to {doc}`Chapter XIV <ch14_investment_uncertainty>` and to
{doc}`36b_exact_linear_re` — between Granger causality and *invariance under intervention*.

## Cagan's regression as a misspecified distributed lag

Under the two sufficient conditions the structural money-demand schedule {eq}`eq-cg-demand` becomes a
one-parameter **distributed-lag regression** of real balances on inflation,

```{math}
:label: eq-cg-struct
m_t - p_t = \alpha\,\pi_t + u_t = \underbrace{\frac{\alpha(1-\lambda)}{1-\lambda L}}_{b^1(L;\,\alpha)}\,x_t + u_t .
```

Cagan read the least-squares fit of {eq}`eq-cg-struct` as delivering the structural $\alpha$. But this
regression is **misspecified**: the disturbance $u_t$ is *not* orthogonal to the inflation process. Under
the equilibrium money rule {eq}`eq-cg-rule`, money growth responds to inflation, so $\pi_t$ (the
regressor) is correlated with $u_t$ (the error) — a simultaneity that biases ordinary least squares.
Computing the actual population projection of $m_t-p_t$ on the whole inflation process — substitute the
Wold representation {eq}`eq-cg-wold` and apply the summation operator $(1-L)^{-1}$ — gives

```{math}
:label: eq-cg-popreg
(m_t-p_t) = \underbrace{\Big[\rho\alpha - (1-\rho)\tfrac{\lambda}{1-\lambda}\Big]}_{\text{plim }\hat\alpha}
            \,\frac{1-\lambda}{1-\lambda L}\,x_t + \bar u_t,
\qquad \bar u_t=\bar u_{t-1}+v_t,
```

with $\bar u_t$ a random walk orthogonal to the $x$ process. Comparing {eq}`eq-cg-popreg` with the
structural form {eq}`eq-cg-struct`:

- Cagan's estimator of the **shape parameter** $\lambda$ is *consistent*;
- Cagan's estimator of the **slope** $\alpha$ is *not*, converging instead to

```{math}
:label: eq-cg-plim
\operatorname{plim}\hat\alpha = \rho\,\alpha - (1-\rho)\,\frac{\lambda}{1-\lambda}.
```

The probability limit is a *weighted average* of the structural $\alpha$ and the Wallace–Sargent value
$-\lambda/(1-\lambda)$, with weight $\rho$. Two polar cases pin down the intuition: if there is no
noise in portfolio balance ($\eta_t\equiv0$, so $\rho=1$), OLS is consistent, $\operatorname{plim}\hat\alpha=\alpha$;
if $\rho=0$, the estimator collapses to the pure Wallace–Sargent value $-\lambda/(1-\lambda)$,
independent of the true $\alpha$. When the true $\alpha$ is more negative than $-\lambda/(1-\lambda)$ —
the empirically relevant case — the bias pulls $\hat\alpha$ *toward zero*, shrinking $-1/\hat\alpha$
below the true revenue-maximizing rate. That is Cagan's paradox, manufactured by misspecification.
(Consistent with the random-walk residual $\bar u_t$ in {eq}`eq-cg-popreg`, both Cagan and Barro
reported highly serially correlated residuals and very low Durbin–Watson statistics.)

## Sims's approximation-error formula

Equation {eq}`eq-cg-plim` is an instance of **Sims's frequency-domain approximation-error formula**.
Recall ([Exercise 1](37_exercises.md) of this chapter, and the seasonal-adjustment application of
{doc}`33a_seasonality_approximation`): if the true projection of $y_t$ on a covariance-stationary
process $x_t$ has transfer function $b^0(L)$, and an econometrician fits, by least squares, a
*constrained* distributed lag $b^1(L)$ drawn from a restricted class, then in population least squares
chooses $b^1$ to minimize the **spectral-density-weighted mean-squared error**

```{math}
:label: eq-cg-sims
\int_{-\pi}^{\pi}\big|\,b^0(e^{-i\omega}) - b^1(e^{-i\omega})\,\big|^{2}\,g_x(e^{-i\omega})\,d\omega ,
```

where $g_x$ is the spectral density of the regressor process. Cagan's regression is precisely such a
constrained fit: with $y_t=m_t-p_t$, the fitted class is the one-parameter geometric-lag family
$b^1(L;\alpha)=\alpha\,\dfrac{1-\lambda}{1-\lambda L}$, and the true projection {eq}`eq-cg-popreg` is the
member of that family with coefficient $\operatorname{plim}\hat\alpha$. Least squares therefore drives
$\hat\alpha$ to the value minimizing {eq}`eq-cg-sims`, which is exactly {eq}`eq-cg-plim`. The
"approximation error" $\operatorname{plim}\hat\alpha-\alpha$ is the gap between the structural
coefficient and the spectral-density-weighted best fit — non-zero because, under the equilibrium money
rule, the regressor is correlated with the error.

Reading the bias through {eq}`eq-cg-sims` makes two things vivid. First, the weighting density $g_x$ is
the **spectrum of inflation**, which is *generated by the money-supply rule* {eq}`eq-cg-rule`: change
the monetary regime and $g_x$ changes, so the "structural" coefficient Cagan estimates is not invariant
across regimes. Sims's approximation formula here wears the clothes of **Lucas's critique** — the same
regime-dependence of a fitted decision/demand rule that runs through {doc}`22_rational_expectations`,
{doc}`Chapter XIV <ch14_investment_uncertainty>`, and {doc}`36b_exact_linear_re`. Second, it locates the
inconsistency where it belongs: not in the *shape* $\dfrac{1-\lambda}{1-\lambda L}$ (which Cagan gets
right) but in the *interpretation* of a projection coefficient as a structural elasticity.

## Consistent estimation, identification, and findings

Because inflation and money growth are determined *simultaneously*, recovering $\alpha$ requires a
system method. The bivariate model {eq}`eq-cg-ma` can be written as a vector ARMA(1,1) whose innovations
$a_t=(a_{1t},a_{2t})'$ are the one-step-ahead forecast errors of $(x_t,\mu_t)$; crucially, $\alpha$ does
**not** enter the innovation recursions, so a Gaussian full-information maximum-likelihood estimator
(Wilson 1973) identifies $\lambda$ and the innovation covariance
$D_a=(\sigma_{11},\sigma_{12},\sigma_{22})$ by minimizing $\lvert\hat D_a(\lambda)\rvert$ over the single
parameter $\lambda$. But the mapping from the four identified quantities
$(\lambda,\sigma_{11},\sigma_{12},\sigma_{22})$ to the five structural parameters
$(\alpha,\lambda,\sigma_\varepsilon^2,\sigma_\eta^2,\sigma_{\varepsilon\eta})$ is not invertible:
$\lambda$ and $\sigma_\varepsilon^2$ are identified, but $\alpha$ and $\sigma_{\varepsilon\eta}$ are
not *separately* identified — offsetting changes in the two leave the likelihood unchanged. To estimate
$\alpha$ at all one must impose a restriction; Sargent sets $\sigma_{\varepsilon\eta}=0$ (money-supply
and portfolio shocks uncorrelated), which yields an estimator of $\alpha$ that depends sensitively on
the estimated covariance matrix of the forecast errors and must be regarded as *delicate*. (The
QuantEcon lecture implements the estimator and, when $\sigma_{\varepsilon\eta}=0$, an equivalent
instrumental-variables procedure that uses the fitted inflation innovations as an instrument for
expected inflation.)

The quantitative findings, summarized from the paper's tables, bear out the "statistical artifact"
reading of the paradox:

- **The estimates of $\alpha$ are very loose.** For most hyperinflations the maximum-likelihood standard
  error on $\hat\alpha$ is of the same order as the point estimate; only **Hungary I** is estimated with
  much precision. For Hungary I, $\hat\alpha\approx-1.84$ implies a revenue-maximizing inflation of
  $-1/\hat\alpha\approx54\%$ per month against an observed $46\%$ — the one case where the estimate
  substantially *weakens* Cagan's paradox. For the others the point estimates do not eliminate the
  paradox, but two-standard-error bands comfortably include values of $\alpha$ that would.
- **The one-parameter model survives overfitting for several countries.** Testing the restricted
  representation — whose systematic part contains the *single* free parameter $\lambda$ — against six
  vector-ARMA overparameterizations by likelihood-ratio $\chi^2$ statistics, the model is *not* rejected
  at the 5% level for **Germany, Greece, and Poland**; Hungary I and Austria are rejected by several
  parameterizations, Russia by one. That a representation with one systematic free parameter survives at
  all is striking.

The upshot: the demand for money during hyperinflations may not have been as sharply isolated as
Cagan's tight estimates suggested, and the slope of the portfolio-balance schedule is difficult or
impossible to estimate precisely under the money-supply regimes that actually prevailed. The apparent
paradox is largely what Sims's approximation-error formula predicts a misspecified least-squares
regression will produce.

## References

- Anderson, T. W. (1971). *The Statistical Analysis of Time Series.* Wiley. [pp. 159–161].
- Barro, R. J. (1970). Inflation, the payments period, and the demand for money. *Journal of Political Economy* 78(6), 1228–1263.
- Cagan, P. (1956). The monetary dynamics of hyperinflation. In M. Friedman (ed.), *Studies in the Quantity Theory of Money*, 25–117. University of Chicago Press.
- Granger, C. W. J. (1969). Investigating causal relations by econometric models and cross-spectral methods. *Econometrica* 37(3), 424–438.
- Muth, J. F. (1960). Optimal properties of exponentially weighted forecasts. *Journal of the American Statistical Association* 55(290), 299–306.
- Sargent, T. J. (1977). The demand for money during hyperinflations under rational expectations: I. *International Economic Review* 18(1), 59–82.
- Sargent, T. J., and N. Wallace (1973). Rational expectations and the dynamics of hyperinflation. *International Economic Review* 14(2), 328–350.
- Sims, C. A. (1972a). Money, income, and causality. *American Economic Review* 62(4), 540–552.
- Sims, C. A. (1972b). Approximate prior restrictions in distributed lag estimation. *Journal of the American Statistical Association* 67(337), 169–175.
- Wilson, G. T. (1973). The estimation of parameters in multivariate time series models. *Journal of the Royal Statistical Society, Series B* 35(1), 76–85.
```
