# Linear Time Series Analysis

This book is a modernized and extended treatment of the linear time series methods that run
through macroeconomics and dynamic econometrics. Its core is **Chapter XI: Time Series** from

> Thomas J. Sargent, *Macroeconomic Theory*, 2nd edition (1987), Academic Press,

now preceded by two chapters — **Chapter IX** and **Chapter X** — that develop the two elementary
tools the time series theory rests on, extended by several new sections, and capped by a full
application chapter, **Chapter XIV — Investment Under Uncertainty**, that puts the whole apparatus to
work in building and interpreting a rational expectations equilibrium.

## A book in three movements, and a capstone

The theoretical material is a single arc that runs from elementary algebra and geometry to a full
theory of linear prediction and its uses in economics. Its logic is a synthesis: **Chapter XI is
what one gets by combining Chapter IX with Chapter X.**

- {doc}`Chapter IX — Difference Equations and Lag Operators <chapters/ch09_difference_equations>`
  is the **algebra**. It introduces the lag operator $L$, the calculus of polynomials in $L$, and
  the use of that calculus to solve linear difference equations — to factor a characteristic
  polynomial into stable and unstable roots, to invert an operator by partial fractions, and to
  solve the linear-quadratic (Euler-equation) optimization problems of dynamic economics by the
  rule *solve stable roots backward and unstable roots forward*.

- {doc}`Chapter X — Linear Least Squares Projections (Regressions) <chapters/ch10_regressions>`
  is the **geometry**. It introduces the linear least squares projection — the regression viewed
  through the orthogonality principle — together with recursive projection (the Kalman filter),
  the law of iterated projections, and the static signal-extraction problem.

- **Chapter XI — Linear Time Series** (the bulk of the book) puts the algebra and the geometry
  together. It studies covariance stationary stochastic processes built from white noise by linear
  difference equations, and it answers, for such processes, the two questions the first two
  chapters were sharpened to ask: *what is the process's structure* (its
  {doc}`spectrum <chapters/06_spectrum>`, its
  {doc}`Wold representation <chapters/13_representation_theory>`), and *how does one predict it*
  (the {doc}`Wiener–Kolmogorov theory <chapters/14_linear_prediction>` of linear prediction,
  {doc}`filtering <chapters/26_optimal_filtering>`, and
  {doc}`signal extraction <chapters/19_signal_extraction>`)?

A fourth element then closes the book: {doc}`Chapter XIV — Investment Under Uncertainty
<chapters/ch14_investment_uncertainty>` is the **capstone**, where the algebra, the geometry, and the
prediction theory are deployed together to construct, compute, and interpret a rational expectations
equilibrium.

## How Chapter XI builds on Chapters IX and X

The two elementary chapters are not prerequisites to be gotten out of the way; their ideas are the
working parts of Chapter XI, reused at every turn. Five threads make the dependence concrete.

**1. Slutsky's reinterpretation: a difference equation driven by chance.** Chapter IX solved
$(1 - a_1 L - \cdots - a_n L^n)\,y_t = x_t$ for a *known* forcing sequence $\{x_t\}$. Chapter XI
opens ({doc}`Section 1 <chapters/01_introduction>`) by making $\{x_t\}$ a sequence of *random*
shocks — a **linear stochastic difference equation**. Slutsky's (1937) insight was that even a
low-order difference equation, if driven by erratic shocks, produces realizations that look like
observed business cycles. Everything Chapter IX taught about roots, stability, and oscillation now
describes the second moments — the {doc}`covariogram <chapters/02_preliminary_concepts>` and
{doc}`spectrum <chapters/06_spectrum>` — of a random process.

**2. Wold's theorem: every stationary process *is* such a difference equation.** The pivot of the
book is {doc}`Representation Theory <chapters/13_representation_theory>`. There, a sequence of the
linear least squares projections of Chapter X — regressing $x_t$ on ever-longer stretches of its
own past — produces the process's **innovations** $\epsilon_t = x_t - P[x_t\mid x_{t-1},x_{t-2},\ldots]$,
and Wold's decomposition theorem shows that *any* covariance-stationary process is a
one-sided moving average $x_t = d(L)\epsilon_t$ of that white noise. This is exactly the fusion of
the two earlier chapters: the object is a stochastic difference equation (Chapter IX), and the
white noise driving it is manufactured by projection (Chapter X). Whittle's spectral
factorization, added to that section, makes the construction *computable* from the spectrum in a
few Fourier transforms.

**3. Prediction is projection organized by lag-operator algebra.** The
{doc}`Wiener–Kolmogorov prediction formula <chapters/14_linear_prediction>` — the central result of
the book — projects the future of a process onto its past. The projection is Chapter X; the
bookkeeping is Chapter IX. Its **annihilation operator** $[\,\cdot\,]_+$ ("discard negative powers
of $L$") is computed pole by pole in the very partial-fraction / first-order
$1/(1-\lambda L)$ basis of Chapter IX — a connection made explicit in
{doc}`The Residue Theorem Behind Partial Fractions <chapters/18a_partial_fractions>`. Selecting the
*fundamental* representation in {doc}`Sections 16–18 <chapters/16_deriving_ma>` is Chapter IX's
root-sorting (roots inside vs. outside the unit circle, in reciprocal pairs) applied to a
covariance generating function.

**4. "Stable roots backward, unstable roots forward" returns as forecasting the future.** The
device that solved Chapter IX's Euler equations reappears throughout Chapter XI whenever agents
must act on forecasts. The {doc}`Hansen–Sargent formula <chapters/20_geometric_leads>` for a
geometric distributed *lead* $P_t\sum_j \lambda^j x_{t+j}$ is the stochastic counterpart of the
forward-looking operator of Chapter IX; it is applied to
{doc}`rational expectations models <chapters/22_rational_expectations>` (Cagan's hyperinflation,
supply and demand with inventories), where factoring $(1-\lambda^{-1}B)(1-\lambda B)$ and inverting
the unstable root forward *inside a projection operator* delivers the equilibrium and its
**cross-equation restrictions**. The {doc}`chain rule of forecasting <chapters/21_chain_rule>` (pure
law of iterated projections) and the {doc}`compact state-space predictor <chapters/25_optimal_prediction>`
$P[x_{t+\tau}\mid x_t]=A^\tau x_t$ are the same idea in recursive and matrix form.

**5. Causality, filtering, and exogeneity are statements about one-sided projections.** Sims's
theorem that {doc}`Wiener–Granger causality equals econometric exogeneity <chapters/27_granger_causality>`
says a projection is *one-sided* — again Chapter X's projection meeting Chapter IX's factorization.
This single idea organizes the applications: Sims's
{doc}`money–income test <chapters/28_sims_money_income>`, the invariance of two-sided (but not
one-sided) projections under {doc}`common filtering <chapters/30_filtering_projections>`, the
{doc}`forward-versus-backward filtering <chapters/31_orthogonality_filtering>` that preserves
orthogonality conditions, and the way {doc}`errors in variables <chapters/35_errors_variables>`
manufacture spurious causality by turning a one-sided projection two-sided.

## The arc of Chapter XI

Read in sequence, Chapter XI moves through five stages:

- **Foundations and the frequency domain** ({doc}`2 <chapters/02_preliminary_concepts>`–{doc}`7 <chapters/07_cross_spectrum>`):
  covariance stationarity, the covariogram and its cross-series analogue, the
  {doc}`Fourier / z-transform <chapters/04_fourier_z_transforms>` machinery (with the
  {doc}`inverse-transform residue calculus <chapters/05_inverse_z_transform>`), and the
  {doc}`spectrum <chapters/06_spectrum>` and {doc}`cross spectrum <chapters/07_cross_spectrum>` —
  gain, phase, and coherence — together with their
  {doc}`FFT-based estimation <chapters/07a_fft_estimation>`.
- **Filters and the business cycle** ({doc}`8 <chapters/08_leading_indicators>`–{doc}`12 <chapters/12_index_models>`):
  what filters do to a spectrum (the {doc}`Slutsky and Kuznets <chapters/09_slutsky_kuznets>`
  spurious cycles), competing {doc}`definitions of the business cycle <chapters/11_business_cycle_definitions>`,
  and the {doc}`index (dynamic-factor) model <chapters/12_index_models>` of comovement.
- **Representation and prediction** ({doc}`13 <chapters/13_representation_theory>`–{doc}`19 <chapters/19_signal_extraction>`):
  Wold's theorem, the Wiener–Kolmogorov formulas, how to build a
  {doc}`Wold representation <chapters/17_wold_ma>` for MA, AR, and ARMA processes, and dynamic
  {doc}`signal extraction <chapters/19_signal_extraction>`.
- **Optimization, expectations, and the vector case** ({doc}`20 <chapters/20_geometric_leads>`–{doc}`26 <chapters/26_optimal_filtering>`):
  geometric leads, the chain rule, rational expectations, vector stochastic difference equations
  and their {doc}`compact notation <chapters/24_compact_notation>`, and optimal prediction and
  filtering.
- **Causality and its econometric pitfalls** ({doc}`27 <chapters/27_granger_causality>`–{doc}`36 <chapters/36_bubbles>`):
  Granger causality, filtering and projections, one-sided-projection theories, seasonal adjustment
  and {doc}`temporal aggregation <chapters/34_aggregation>`, errors in variables,
  {doc}`rational bubbles <chapters/36_bubbles>`, and — the section that Chapter XIV speaks to
  directly — {doc}`the difficulty of interpreting vector autoregressions <chapters/36a_interpreting_vars>`,
  where the equilibrium of a market with forward-looking supply and demand is read as a vector
  autoregression.

## Chapter XIV: a rational expectations equilibrium in action

{doc}`Chapter XIV — Investment Under Uncertainty <chapters/ch14_investment_uncertainty>` is where the
whole apparatus is put to work. It extends the linear-quadratic Euler-equation problems of
{doc}`Chapter IX <chapters/ch09_difference_equations>` to *stochastic* forcing processes — the
certainty-equivalence principle lets forecasting and optimization separate — solves them with the
Wiener–Kolmogorov / Hansen–Sargent geometric-lead formula {eq}`eq-90` of Chapter XI, and assembles
the pieces into Lucas and Prescott's (1971) model of a competitive industry. Its centerpiece is a
precise, constructive account of a **rational expectations equilibrium** and *two* complementary ways
to compute and interpret one:

- **As a fixed point.** Each firm forecasts the endogenous output price by projecting on the very
  price process that all firms' investment decisions jointly generate. The equilibrium is a *fixed
  point* of the mapping from the price law of motion firms **perceive** to the one their behavior
  **actually** produces (Chapter XIV, §§6, 9). This is the operational content of "rational
  expectations": the agents in the model forecast prices as well as the economist who models them.

- **As a social planning problem.** Because the competitive industry equilibrium implicitly maximizes
  a welfare criterion — discounted consumer surplus minus producer surplus — one can replace the
  awkward fixed-point calculation with a straightforward maximization (Chapter XIV, §5). Lucas and
  Prescott's device turns "find the equilibrium" into "solve a planning problem," a trick that
  pervades modern macroeconomics.

**The link back to Chapter XI's dynamic supply and demand curves.** From the equilibrium Chapter XIV
(§7) reads off a *dynamic supply curve*: current output depends on lagged output and on current and
**expected future** prices, so — because forecasting future prices requires the parameters of the
demand process — *the demand curve's parameters appear inside the supply curve*, subverting the
exclusion restrictions that ordinarily identify a supply schedule. The identifying information that
remains lives entirely in the **cross-equation restrictions** that a rational expectations
equilibrium stamps onto the data. This is the very same object studied from the opposite side in
{doc}`A Difficulty in Interpreting Vector Autoregressions <chapters/36a_interpreting_vars>` (Chapter
XI's Hansen–Sargent section): there the equilibrium of a market with **forward-looking supply and
demand** is a covariance-stationary vector process — a vector autoregression — and the lesson is that
its *Wold innovations are generally not the structural supply and demand shocks* that hit agents, so
the innovation accounting of a fitted VAR misattributes its shocks. The two sections are two views of
one equilibrium: Chapter XIV **constructs** the dynamic supply and demand curves and shows why their
parameters are entangled; {doc}`§36a <chapters/36a_interpreting_vars>` takes such an equilibrium as
given and shows why a vector autoregression cannot recover the agents' surprises from it. Both derive
the forward-looking decision rules by the same *stable-roots-backward / unstable-roots-forward*
factorization of {doc}`Chapter IX <chapters/ch09_difference_equations>`, and both turn on the
cross-equation restrictions of {doc}`rational expectations <chapters/22_rational_expectations>` — the
thread that also runs through {doc}`exact linear rational expectations models <chapters/36b_exact_linear_re>`.
Chapter XIV closes with a precise statement of **Lucas's critique**: because decision rules like the
dynamic supply curve inherit the parameters of the processes agents forecast, they are not invariant
to policy interventions in those processes — the reason estimation must reach for the "deep"
parameters of preferences and technology.

## What is new in this edition

Beyond modernizing the 1987 text, this edition adds several sections that extend the theory or
connect it to recent work:

- {doc}`The uncertainty principle for Fourier transform pairs <chapters/05a_uncertainty_principle>` —
  the time–frequency trade-off that limits what any filter or spectral window can resolve.
- {doc}`Estimating spectra, cross spectra, and bispectra with the FFT <chapters/07a_fft_estimation>` —
  following Hinich and Clay (1968), the practical estimation counterpart to the spectrum and
  cross-spectrum theory: the periodogram's inconsistency, the resolution–variance trade-off, and the
  **bispectrum**, whose non-vanishing is a fingerprint of nonlinearity a flat spectrum cannot see.
- {doc}`The residue theorem behind partial fractions <chapters/18a_partial_fractions>` — the
  complex-analysis foundation of the partial-fraction calculus used throughout the prediction and
  signal-extraction chapters, showing why the annihilation operator acts one pole at a time.
- **Whittle's spectral factorization** (in
  {doc}`Representation Theory <chapters/13_representation_theory>`) — a constructive, FFT-based way
  to recover the Wold moving-average kernel and the innovation variance from a spectral density;
  the computational complement to the existence theorems of
  {doc}`Sections 16–18 <chapters/17_wold_ma>`.
- {doc}`Seasonality and approximation errors <chapters/33a_seasonality_approximation>` — following
  Hansen and Sargent (1993), a frequency-domain criterion for when a misspecified model is better
  estimated with seasonally adjusted data.
- {doc}`Sims's formula, derived in the time domain <chapters/34a_sims_expository_note>` — a
  self-contained pedagogical companion to the {doc}`aggregation-over-time <chapters/34_aggregation>`
  section, tying Sims's discrete/continuous-time result to Theil's specification theorem.
- {doc}`A difficulty in interpreting vector autoregressions <chapters/36a_interpreting_vars>` and
  {doc}`exact linear rational expectations models <chapters/36b_exact_linear_re>` — following Hansen
  and Sargent (1991), examples in which the innovations a vector autoregression recovers are *not*
  the shocks that hit agents, with dynamic supply and demand curves derived by the
  stable-roots-backward / unstable-roots-forward method of
  {doc}`Chapter IX <chapters/ch09_difference_equations>` (the counterpart, from the econometrician's
  side, to the equilibrium **constructed** in Chapter XIV).
- {doc}`Money demand in hyperinflations <chapters/36c_cagan_hyperinflation>` — a capstone application,
  following Sargent (1977), that reads Cagan's money-demand regression as a *misspecified distributed
  lag* under rational expectations. It shows that the notorious inconsistency of Cagan's estimator is an
  instance of **Sims's frequency-domain approximation-error formula** (Exercise 1 of Chapter XI, and
  {doc}`§33a <chapters/33a_seasonality_approximation>`), and it ties together the
  {doc}`rational-expectations Cagan model <chapters/22_rational_expectations>`,
  {doc}`Granger causality <chapters/27_granger_causality>` versus mere *leading*
  ({doc}`§08 <chapters/08_leading_indicators>`), and the regime-dependence at the heart of Lucas's
  critique.
- {doc}`Decomposing an explosive autoregression <chapters/36d_explosive_decomposition>` — a short note
  that applies the chapter's "unstable roots forward" device to a scalar explosive AR(1), splitting it
  into an explosive trend with a random amplitude plus a *stationary anticipative* component. It is the
  explosive-root counterpart of the Beveridge–Nelson permanent/transitory decomposition and sets up
  {ref}`Exercise 4 <ex-4>`.
- Two postscripts that move beyond the linear, stationary theory:
  {doc}`nonlinear (Volterra / Wiener–Itô) moving-average representations <chapters/39_nonlinear_representation>`,
  where a nonzero bispectrum is the fingerprint of nonlinearity a flat spectrum cannot see; and
  {doc}`complex demodulation <chapters/41_comp_demod>`, a tool for estimating a *moving* spectrum and
  cross spectrum, applied to the changing seasonal in U.S. interest rates.

## What this book adds throughout

Relative to the 1987 original, this version:

- **Corrects typographical and mathematical errors** present in the LaTeX source.
- **Adds Python code** that generates modern versions of every figure, using current U.S. and
  international data; each figure links to the script that produced it.
- **Adds worked solutions** to the exercises of Chapters IX, X, XI, and XIV as collapsible dropdowns.
- **Extends several sections** with updated empirical examples and commentary, and
  **cross-references QuantEcon lectures** where related code already exists
  (see [ARMA](https://python-advanced.quantecon.org/arma.html) and
  [Spectral Estimation](https://python-advanced.quantecon.org/estspec.html)).

## How to read this book

A reader new to the material can proceed linearly: Chapters IX and X first, then Chapter XI, then the
Chapter XIV capstone. A reader already comfortable with lag operators and regressions can begin
directly at {doc}`chapters/01_introduction` and refer back to
{doc}`Chapter IX <chapters/ch09_difference_equations>` and
{doc}`Chapter X <chapters/ch10_regressions>` as needed;
{doc}`Chapter XIV <chapters/ch14_investment_uncertainty>` can then be read as a self-contained
application once the prediction theory of {doc}`Sections 14 <chapters/14_linear_prediction>` and
{doc}`20 <chapters/20_geometric_leads>` is in hand. Three mathematical or estimation-oriented
digressions — {doc}`Fourier and z-transforms <chapters/04_fourier_z_transforms>`,
{doc}`the uncertainty principle <chapters/05a_uncertainty_principle>`, and
{doc}`FFT estimation <chapters/07a_fft_estimation>` — are marked *optional on first reading* and may
be skipped without loss of continuity. The {doc}`Index <chapters/42_index>` at the end collects the
book's concepts, named results, and people with links to the sections where they appear.

## Notation

Throughout, $L$ denotes the **lag operator**, $Lx_t = x_{t-1}$, introduced in
{doc}`Section 1 of Chapter IX <chapters/ch09_difference_equations>`. $E$ denotes the mathematical
expectation operator, and $\hat E$ or $P[\,\cdot \mid \cdot\,]$ the **linear least squares
projection** operator of {doc}`Chapter X <chapters/ch10_regressions>`. The operator $[\,\cdot\,]_+$
is the **annihilation operator** (discard negative powers of $L$) of the Wiener–Kolmogorov formula.
A polynomial $d(z)$ is **fundamental** when its zeros lie outside the unit circle, so that
$d(L)^{-1}$ is one-sided in nonnegative powers of $L$. All stochastic processes are discrete time
unless otherwise stated.

## References

- Sargent, T.J. (1987). *Macroeconomic Theory*, 2nd ed. Academic Press.
- Lucas, R.E. Jr., and E.C. Prescott (1971). Investment under uncertainty. *Econometrica* 39(5),
  659–681.
- Hansen, L.P., and T.J. Sargent (1991). *Rational Expectations Econometrics*. Westview Press.
- Hansen, L.P., and T.J. Sargent (1993). Seasonality and approximation errors in rational
  expectations models. *Journal of Econometrics* 55, 21–55.
- Hinich, M.J., and C.S. Clay (1968). The application of the discrete Fourier transform in the
  estimation of power spectra, coherence, and bispectra of geophysical data. *Reviews of Geophysics*
  6(3), 347–363.
- Whittle, P. (1983). *Prediction and Regulation*, 2nd ed. University of Minnesota Press.
- Wold, H. (1938). *A Study in the Analysis of Stationary Time Series*. Almqvist & Wiksell.
- Anderson, T.W. (1971). *The Statistical Analysis of Time Series*. Wiley.
- Slutsky, E. (1937). The summation of random causes as the source of cyclic processes.
  *Econometrica* 5, 105–146.
