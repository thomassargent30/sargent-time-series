# Introduction

````{admonition} About this part of the book — and its debt to Chapter XI
:class: note

The nucleus of what follows is **Chapter XI, "Time Series," of Thomas J. Sargent,
*Macroeconomic Theory*, 2nd edition (Academic Press, 1987)**. Sections 1–38 below follow that
chapter's development and, in most places, its text — corrected, re-typeset, and re-illustrated with
modern data and code.

But this is no longer a reprint of that chapter, and it would mislead the reader to call it one.
Around the 1987 nucleus have grown a dozen new sections that take up questions the original chapter
raised but did not pursue, or that had not yet been written when it appeared: the
{doc}`uncertainty principle <05a_uncertainty_principle>` for transform pairs,
{doc}`FFT-based estimation <07a_fft_estimation>` of spectra and bispectra, the
{doc}`residue theorem <18a_partial_fractions>` behind the partial-fraction calculus,
{doc}`seasonality and approximation error <33a_seasonality_approximation>`,
{doc}`Sims's aggregation formula in the time domain <34a_sims_expository_note>`, two chapters of
Hansen and Sargent's *Rational Expectations Econometrics* on
{doc}`interpreting vector autoregressions <36a_interpreting_vars>` and
{doc}`exact linear rational expectations models <36b_exact_linear_re>`,
{doc}`money demand in hyperinflations <36c_cagan_hyperinflation>`,
{doc}`explosive autoregressions <36d_explosive_decomposition>`,
{doc}`Lucas's quantity-theory illustrations and Whiteman's critique <36e_lucas_whiteman_quantity_theory>`,
and two postscripts on {doc}`nonlinear representations <39_nonlinear_representation>` and
{doc}`complex demodulation <41_comp_demod>`. Together the added sections now account for roughly a
third of this part, and — counting the two preparatory chapters and the Chapter XIV capstone that
frame it — the 1987 chapter supplies well under half of the book.

So read what follows as material **in the spirit of** Chapter XI, and **applying its methods** to a
wider range of problems than that chapter could take up. The organizing ideas are the ones Chapter XI
was built on — covariance stationarity, the spectrum, Wold's theorem, Wiener–Kolmogorov prediction,
and the one-sided projection — and the added sections are there because those ideas kept paying off.

```{admonition} A note on the numbering
:class: tip, dropdown

The section numbers preserve the 1987 chapter's order so that readers of *Macroeconomic Theory* can
navigate by them. Sections **1–38** are the original sequence. Sections carrying a **letter suffix**
(5a, 7a, 18a, 33a, 34a, 36a–36e) are new, and are placed immediately after the original section whose
argument they extend. Sections **39–41** are postscripts that leave the linear, stationary theory
behind.
```
````

In {doc}`Chapter IX <ch09_difference_equations>`, we studied linear difference equations of the form

$$
(1-a_1 L - \cdots - a_n L^n)y_t = x_t,
$$ (eq-1)

where $\{x_t\}_{t=-\infty}^{\infty}$ was taken to be a known sequence. We studied how to find the class of sequences $\{y_t\}$ that satisfy the difference equation and a set of prescribed boundary conditions on the $\{y_t\}$ sequence. Such a $\{y_t\}$ sequence was said to solve the difference equation.

The present chapter studies linear difference equations of the form {eq}`eq-1` in which, rather than being a sequence of unknown numbers, $\{x_t\}$ is a sequence of independently and identically distributed random variables with known variance and mean. With this choice of mechanism for generating $\{x_t\}$, equation {eq}`eq-1` is called a **linear stochastic difference equation**. A solution of such a difference equation is a sequence of random variables $\{y_t\}$. A sequence of random variables is called a *stochastic process*. While the $x_t$ sequence is by assumption a stochastic process consisting of random variables that are independently and identically distributed over time, the $y_t$ process that solves {eq}`eq-1` will in general be correlated over time. That is, while the $\{x_t\}$ process by assumption satisfies

$$
E(x_t - Ex_t)(x_{t+s}-Ex_{t+s}) = 0
$$

for $s \neq 0$, for the $y_t$ process in general $E(y_t - Ey_t)(y_{t+s}-Ey_{t+s}) \neq 0$ for $s \neq 0$. One way to characterize the solution of the difference equation is to summarize the second moments of the $\{y_t\}$ process and to describe how they depend on the $a_j$'s of {eq}`eq-1`.

Stochastic difference equations provide a natural tool for interpreting and modeling economic time series. Macroeconomists spend much of their time interpreting sample first and second moments of observed time series. For example, for an observed sample on two variables $(y_t, z_t,\ t=1,\ldots,T)$ we often calculate various sample moments

$$
T^{-1}\sum_{t=1}^{T}y_t,\quad T^{-1}\sum_{t=1}^{T}z_t,\quad
(T-k)^{-1}\sum_{t=k+1}^{T}y_t y_{t-k},\quad
(T-k)^{-1}\sum_{t=k+1}^{T}y_t z_{t-k},
$$

and

$$
(T-k)^{-1}\sum_{t=k+1}^{T}z_t y_{t-k}
$$

for various values of $k$. It is convenient to adopt a mathematical context in which these sample moments can be regarded as estimators of the population moments $Ey_t,\ Ez_t,\ Ey_t y_{t-k},\ Ey_t z_{t-k}$, and $Ez_t y_{t-k}$, respectively, estimators which converge to these population moments as $T\to\infty$. Linear stochastic difference equations provide such a mathematical context. In studying how to solve stochastic difference equations, one of our intermediate goals is to learn how the coefficients $a_j$ of {eq}`eq-1` can be chosen in order to make the implied pattern of population moments $Ey_t y_{t-k}$ resemble $(T-k)^{-1}\sum_{t=k+1}^{T}y_t y_{t-k}$ as the lag $k$ is varied.

Stochastic processes provide a natural context in which to formulate the problem of prediction. At time $t$, suppose that observations on a stochastic process $(y_{t+1},y_{t+2},\ldots)$ have not yet been revealed, but that observations on $(y_t,y_{t-1},\ldots)$ are available. Suppose that the moments $Ey_t$ and $Ey_t y_{t-k}$ are known for all $t$ and $k$. Then what is the best way to predict $(y_{t+1},y_{t+2},\ldots)$ as a linear function of $(y_t, y_{t-1}, \ldots)$? This linear prediction problem was solved by Wiener and Kolmogorov in the late 1930s; its solution is built on the linear least squares projection developed in {doc}`Chapter X <ch10_regressions>`.

The linear prediction problem is of interest to macroeconomists for at least two reasons. First, macroeconomists are interested in modeling the behavior of agents who are operating in dynamic and uncertain contexts. Typically, the hypothesis of utility or profit maximization ends up confronting those agents with some version of a prediction problem that they must solve in order best to achieve their objective. As we shall see, by using prediction theory, it is possible to extend the solutions of the quadratic dynamic optimization problems that were encountered in {doc}`Chapter IX <ch09_difference_equations>` to the case in which the forcing functions are stochastic processes whose future values are not known at the time when decisions must be made. Thus, prediction theory is an important tool in determining optimization behavior under uncertainty.

Second, macroeconomists are interested in using their own models of economic time series (often a collection of estimated $a_j$'s in {eq}`eq-1` or estimated moments $Ey_t y_{t-k}$) in order to predict the future conditional on the past. When the econometric model occurs in the form of a vector version of {eq}`eq-1`, it is said to be a *vector autoregression*. Linear prediction theory applies directly to such a model.

One of the goals of much recent work in rational expectations economics has been to create models whose equilibria are vector stochastic difference equations. In these models, the outcome of the interaction of a collection of purposeful agents is a stochastic process for, say, prices and quantities whose evolution can be described by a (vector) stochastic difference equation. We shall study versions of such models in which the equilibria are described by linear stochastic difference equations, i.e., vector versions of {eq}`eq-1`. In such models, some of the $a_j$'s become interpretable in terms of purposeful behavior of the agents in the model; that is, they are functions of the parameters of people's objective functions and constraints. One goal of this line of research is to acquire the ability to predict how the equilibrium stochastic process (or difference equation) would change in response to hypothetical changes in particular aspects of the environment confronting the agents in the model.

The idea that low-order linear stochastic difference equations could provide a useful model for business cycles can be traced back at least as far as {cite:t}`slutsky1937summation` and {cite:t}`frisch1933propagation`. We have seen in Chapter IX that low-order nonstochastic linear difference equations with no forcing functions present (i.e., $x_t=0$ for all $t$ in {eq}`eq-1`) result in solutions for $y_t$ that are "smooth", being the weighted sum of a small number of geometric sequences. Such smooth sequences do not resemble observed economic time series. However, if a sufficiently erratic forcing sequence $\{x_t\}$ occurs in {eq}`eq-1`, the resulting $\{y_t\}$ sequence can be sufficiently erratic that it resembles observed economic time series. The idea of Slutsky was to make the $\{x_t\}$ sequence sufficiently erratic by choosing it as the realization of a sequence of independently and identically distributed random variables. The resulting realizations of the $\{y_t\}$ sequence that solves[^fn-intro-1] {eq}`eq-1` would be erratic enough to resemble observed time series. As we shall see, even first-order stochastic linear difference equations ($n=1$ in {eq}`eq-1`) can generate realizations that look like observed economic time series. Furthermore, the hypothesis that $\{x_t\}$ is a sequence of independently and identically distributed random variables in general implies that the future values $(y_{t+1},y_{t+2},\ldots)$ are at best imperfectly predictable from past values $(y_t,y_{t-1},\ldots)$. It is desirable to have models in which both economic agents and econometricians confront uncertainty in this sense. This is one major reason that Slutsky's idea was adopted early in dynamic econometrics, and why it has been retained and expanded upon in work on rational expectations.

[^fn-intro-1]: Here "solve" refers to the ordinary sense used in Chapter IX of finding a $\{y_t\}$ sequence that satisfies {eq}`eq-1` given the realization of the $\{x_t\}$ sequence.
