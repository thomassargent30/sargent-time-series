# Estimating Spectra, Cross Spectra, and Bispectra with the FFT

```{note}
This section recasts the technical review of Hinich and Clay, "The application of the discrete
Fourier transform in the estimation of power spectra, coherence, and bispectra," *Reviews of
Geophysics* 6(3), 347–363 (1968), for economic time series. Their examples were geophysical
(atmospheric pressure, ocean waves); we keep the statistics and drop the geophysics. Throughout we
set the sampling interval to $\Delta=1$ (one observation per period), so frequency $f$ runs over
$[-\tfrac12,\tfrac12]$ and the angular frequency of {doc}`06_spectrum` is $\omega=2\pi f$; reinstating
$\Delta$ merely rescales frequencies to cycles per unit time.
```

The previous sections defined the {doc}`spectrum <06_spectrum>` and {doc}`cross spectrum <07_cross_spectrum>`
as Fourier transforms of population second moments. This section is about *estimating* them from a
finite sample — and about the object that second moments cannot see, the **bispectrum**. The unifying
tool is the **fast Fourier transform (FFT)**: once the data have been transformed, the power spectrum,
the cross spectrum, the coherence, and the bispectrum are all simple products and averages of the same
Fourier coefficients, computable in $O(n\log n)$ operations rather than the $O(n^2)$ of the older
autocovariance-then-transform route.

## 1. The discrete Fourier transform and the periodogram

Given a finite record $X_0,\ldots,X_{n-1}$ of a covariance-stationary process, the discrete Fourier
transform computes the $n$ complex coefficients

```{math}
:label: eq-fft-dft
A_k = \frac{1}{\sqrt{n}}\sum_{t=0}^{n-1} X_t\, e^{\,2\pi i k t/n}, \qquad k=0,1,\ldots,n-1,
```

with inverse $X_t = \frac{1}{\sqrt n}\sum_{k=0}^{n-1}A_k\,e^{-2\pi i t k/n}$. The FFT (Cooley and Tukey,
1965) evaluates {eq}`eq-fft-dft` for all $k$ at once in $O(n\log n)$ operations. Since the $X_t$ are
real, $A_k=A_{n-k}^{*}$, so $|A_k|^2=|A_{n-k}|^2$; the $n/2$ numbers

$$
I_k \equiv |A_k|^2, \qquad k=0,1,\ldots,\tfrac{n}{2}-1,
$$

are the **periodogram ordinates** at the *Fourier frequencies* $f_k=k/n$.

The periodogram is the Fourier transform of the *sample* autocovariance. Straightforward algebra gives

```{math}
:label: eq-fft-perio
I_k = \sum_{r=-n+1}^{\,n-1}\Big(1-\tfrac{|r|}{n}\Big)\,C_r\, e^{-2\pi i r k/n},
\qquad
C_r = \frac{1}{n-r}\sum_{t=0}^{n-r-1}X_t X_{t+r} \ \ (r\ge0),\quad C_{-r}=C_r .
```

Writing $\mu=EX_t$ and $\rho_r=E(X_tX_{t+r})-\mu^2$ for the population autocovariances, and recalling that
the spectrum is $S(f)=\sum_{r}\rho_r e^{-2\pi i r f}$ (the {doc}`06_spectrum` of the process, with
$\omega=2\pi f$), taking the expectation of {eq}`eq-fft-perio` gives, for $k\ge1$,

```{math}
:label: eq-fft-Eperio
E(I_k) = \int_{-1/2}^{1/2} S(f_k-f)\,\frac{\sin^2 n\pi f}{n\sin^2 \pi f}\,df
       \;\approx\; S(f_k)\quad(\text{large }n),
```

while $E(I_0)=S(0)+\mu^2 n$ picks up the mean. Two lessons follow.

**Leakage (bias).** The periodogram's expectation is the true spectrum *convolved with the Fejér kernel*
$W_n(f)=\sin^2(n\pi f)/(n\sin^2\pi f)$. The kernel integrates to one and concentrates as $n\to\infty$, so
$I_k$ is asymptotically unbiased; but for finite $n$ its side lobes **leak** power from frequencies where
$S$ is large into neighboring bands. In economic data the leakage culprits are a **trend** and a
**seasonal**: both put enormous power at $f=0$ and the seasonal frequencies, which then bleeds across the
whole band (see Section 2). The zero-frequency mean is removed by subtracting the sample mean
$\bar X=\frac1n\sum_t X_t$ — equivalently by setting $A_0=0$.

**Inconsistency (variance).** For a linear (not necessarily Gaussian) process the periodogram ordinates
are, for large $n$, approximately *independent* across $k$, with

```{math}
:label: eq-fft-varperio
\operatorname{var}(I_k) \approx S(f_k)^2 \quad (k\ge1),
\qquad
\frac{2\,I_k}{S(f_k)} \ \dot\sim\ \chi^2_2 .
```

The variance **does not shrink with $n$**: the periodogram is *not* a consistent estimator of the
spectrum. Its standard deviation equals 100% of its mean, so $I_k$ can land anywhere from $0$ to
$2S(f_k)$, and a plot of $I_k$ against $k$ looks wildly erratic, riddled with spurious peaks. Adding data
buys *resolution* (more, finer-spaced Fourier frequencies) but never *precision* at any one frequency.
This is the fact that forces every practical spectral estimator to trade resolution for variance.

## 2. Trading resolution for variance

The remedy is to average. Because the periodogram ordinates are asymptotically independent, averaging
many of them cuts the variance in proportion to the number averaged — at the cost of coarser frequency
resolution. Three classic devices, all equivalent in large samples, implement the trade-off.

**Averaged periodogram (Bartlett; Welch).** Split the record into $r$ non-overlapping segments of length
$m=n/r$ and average the segment periodograms:

```{math}
:label: eq-fft-bartlett
S_k = \frac{1}{r}\sum_{p=1}^{r}\big|A_k^{(p)}\big|^2,
\qquad
A_k^{(p)} = \frac{1}{\sqrt m}\sum_{t=0}^{m-1} X_{t+(p-1)m}\, e^{\,2\pi i t k/m},
\qquad k=0,\ldots,\tfrac m2-1 .
```

Then $E(S_k)\approx S(f_k)$ with $f_k=k/m$, and, since the $r$ segments are approximately independent,

```{math}
:label: eq-fft-bartvar
\operatorname{var}(S_k)\approx \frac{1}{r}\,S(f_k)^2 = \frac{m}{n}\,S(f_k)^2 .
```

Resolution is now $1/m$ (coarser than the periodogram's $1/n$), and the proportional variance is $1/r=m/n$
(smaller). The choice of the segment length $m$ *is* the resolution–variance dial.

**Blackman–Tukey.** Transform only the first $M=m/2$ sample autocovariances,
$U_k=\sum_{r=-M}^{M}C_r\,e^{-2\pi i r k/m}$. Then $E(U_k)\approx S(f_k)$ and
$\operatorname{var}(U_k)\approx (m/n)\,S(f_k)^2$ — the same large-sample bias and variance as the averaged
periodogram. Here the implied window is the **Dirichlet kernel** $\sin[(m{+}1)\pi f]/\sin(\pi f)$, whose
side lobes decay only as $O(f^{-1})$ (versus the Fejér kernel's $O(f^{-2})$), so Blackman–Tukey leaks
*more* from sharp peaks than the averaged periodogram does.

**Smoothing (Hanning) and spectral windows.** One can also smooth adjacent ordinates,
$V_k=\sum_j c_j\,U_{k+j}$ with $\sum_j c_j=1$; *Hanning* uses $c_0=\tfrac12,\ c_{\pm1}=\tfrac14$. Then

```{math}
:label: eq-fft-hanning
\operatorname{var}(V_k)\approx \big(c_{-1}^2+c_0^2+c_1^2\big)\frac{m}{n}\,S(f_k)^2
   = \tfrac{3}{8}\,\frac{m}{n}\,S(f_k)^2 ,
```

so Hanning cuts the variance to $3/8$ of the raw estimate, at the price of a wider effective bandwidth.
The three kernels trade off the same way: the Fejér kernel has the smallest bandwidth and lowest side
lobes (least bias) but the averaged periodogram $S_k$ and the Hanned $V_k$ differ in how much variance
reduction they buy for a given bandwidth. **This bias–variance / bandwidth–variance conflict is the
uncertainty principle of {doc}`05a_uncertainty_principle`**: a window narrow in frequency (fine
resolution) is wide in time (few independent averages, high variance), and conversely. There is no free
lunch; one chooses where to sit on the trade-off given the sample length.

**Trends and seasonals.** Suppose the observed series is
$Y_t = m_t + \sum_{\ell} a_\ell e^{2\pi i f_\ell t} + X_t$, a low-order polynomial trend $m_t$ plus a few
sinusoids (a seasonal) plus the stationary part $X_t$ we want. The trend and the sinusoids put a Dirac
spike at $f=0$ and at the seasonal frequencies, and the kernel side lobes spread that power across the
band, **biasing the whole spectrum estimate**. So before estimating $S_X$ one must remove the
low-frequency and seasonal power — by first-differencing (a high-pass filter, then correcting the
estimate by the filter's squared gain, cf. {doc}`09_slutsky_kuznets` and {doc}`10_filter_kit`) or by a
least-squares fit of the trend and seasonal dummies. This is exactly the concern of
{doc}`seasonal adjustment <33_seasonal_adjustment>`: mis-handling the trend and seasonal manufactures
spurious structure in the estimated spectrum, just as an ill-chosen filter manufactures the spurious
cycles of {doc}`09_slutsky_kuznets`.

## 3. Cross spectra and coherence by the FFT

For two jointly stationary series $\{X_t\},\{Y_t\}$ the {doc}`cross spectrum <07_cross_spectrum>` is
$S_{xy}(f)=\sum_r \rho_r^{xy}e^{-2\pi i r f}$ with $\rho_r^{xy}=E(X_t-\mu_x)(Y_{t+r}-\mu_y)$. Its polar
decomposition defines the **coherence** and **phase**,

```{math}
:label: eq-fft-coh
\gamma_{xy}(f) = \frac{|S_{xy}(f)|}{\sqrt{S_x(f)\,S_y(f)}}\in[0,1],
\qquad
\phi(f) = \tan^{-1}\!\frac{\operatorname{Im}S_{xy}(f)}{\operatorname{Re}S_{xy}(f)},
\qquad
S_{xy}=\gamma_{xy}\sqrt{S_xS_y}\,e^{i\phi}.
```

The coherence $\gamma_{xy}(f)^2$ is the **frequency-by-frequency $R^2$** between the two series — the
fraction of the variance of the $Y$-component at frequency $f$ that is linearly predictable from the
$X$-component at that frequency — and the phase measures their lead/lag there (the comovement notion
behind the business-cycle definition of {doc}`11_business_cycle_definitions`). Estimation mirrors the
univariate case. Form the segment transforms $A_k^{(p)}$ of $X$ and $B_k^{(p)}$ of $Y$ as in
{eq}`eq-fft-bartlett`, and average the **cross-periodogram**,

```{math}
:label: eq-fft-crossest
S_k^{xy} = \frac{1}{r}\sum_{p=1}^{r} A_k^{(p)}\,B_k^{*(p)},
\qquad E(S_k^{xy})\approx S_{xy}(f_k),
```

then read off the coherence and phase estimates $\hat\gamma_{xy}(f_k)=|S_k^{xy}|/\sqrt{S_k^x S_k^y}$ and
$\hat\phi(f_k)=\tan^{-1}[\operatorname{Im}S_k^{xy}/\operatorname{Re}S_k^{xy}]$. Their large-sample
variances are

```{math}
:label: eq-fft-cohvar
\operatorname{var}(\hat\gamma_{xy}) \approx \frac{m}{2n}\big[1-\gamma_{xy}^2(f_k)\big]^2,
\qquad
\operatorname{var}(\hat\phi) \approx \frac{m}{2n}\big[\gamma_{xy}^{-2}(f_k)-1\big].
```

The phase is well determined only where the coherence is high: as $\gamma_{xy}\to0$ the phase variance
blows up, so a measured lead/lag between two nearly-incoherent series carries no information — a warning
worth heeding before reading economic significance into an estimated phase lead (compare the
leading-indicator caution of {doc}`08_leading_indicators`).

## 4. Bispectra by the FFT

Everything above lives at the level of *second* moments, and second moments are blind to nonlinearity. If
$\{X_t\}$ is Gaussian — or, more generally, linear in Gaussian innovations — its spectrum and cross
spectrum exhaust its probabilistic structure. To detect a **quadratic** interaction among frequency
components one needs a *third*-moment object, the **bispectrum**.

The motivating phenomenon is **quadratic phase coupling**. Pass a process through a nonlinear filter
(say a squarer). If it contains strong components at $f_1$ and $f_2$ with phases $\phi(f_1),\phi(f_2)$,
the nonlinearity creates a component at $f_1+f_2$ whose phase is

```{math}
:label: eq-fft-coupling
\phi(f_1+f_2) = \phi(f_1)+\phi(f_2)-\theta ,
```

with $\theta$ nearly constant. The bispectrum detects this phase coherence among the triple
$(f_1,f_2,f_1{+}f_2)$. It is the double Fourier transform of the **third-order cumulant**
$C(\sigma,\tau)=E(X_t X_{t+\sigma}X_{t+\tau})$ (for a mean-zero process),

```{math}
:label: eq-fft-bispec
B(f_1,f_2) = \sum_{\sigma}\sum_{\tau} C(\sigma,\tau)\,e^{-2\pi i(f_1\sigma+f_2\tau)},
\qquad B(f_1,f_2)=B(f_2,f_1)=B^{*}(-f_1,-f_2).
```

Two facts make it a diagnostic. **First**, if $\{X_t\}$ is Gaussian then $C(\sigma,\tau)\equiv0$, so
$B\equiv0$; *a non-zero bispectrum therefore certifies non-Gaussianity* — and, since a linear process
driven by Gaussian noise also has $B\equiv0$, a non-zero bispectrum is evidence of **nonlinearity**. This
is the basis of Hinich's linearity and Gaussianity tests and the direct empirical counterpart of the
{doc}`nonlinear (Volterra/Wiener–Itô) theory <39_nonlinear_representation>`, where the bispectrum is the
first polyspectrum beyond the ordinary spectrum. **Second**, its polar form separates magnitude from
phase,

```{math}
:label: eq-fft-bicoh
B(f_1,f_2) = \big[S(f_1)\,S(f_2)\,S(f_1+f_2)\big]^{1/2}\,\rho(f_1,f_2)\,e^{\,i\theta(f_1,f_2)},
```

where $\rho\in[0,1]$ is the **skewness** or **bicoherence** (the degree of quadratic coupling, normalized
like a coherence) and $\theta$ is the **biphase**.

**FFT estimate.** Exactly as the cross-periodogram averages products of *two* Fourier coefficients, the
bispectrum averages *triple* products,

```{math}
:label: eq-fft-bispecest
\hat B(f_j,f_k) = \frac{1}{r}\sum_{p=1}^{r} A_j^{(p)}\,A_k^{(p)}\,A_{j+k}^{*(p)},
\qquad 0<k<j<\tfrac m2 ,
```

with the skewness estimate $\hat\rho(f_j,f_k)=|\hat B(f_j,f_k)|/\sqrt{S_j S_k S_{j+k}}$. For large $m,n$
the real and imaginary parts of $\hat B$ are independent with

```{math}
:label: eq-fft-bispecvar
\operatorname{var}(\operatorname{Re}\hat B)=\operatorname{var}(\operatorname{Im}\hat B)
   \approx \frac{m}{2n}\,S(f_j)\,S(f_k)\,S(f_{j+k}),
\qquad
\operatorname{var}(\hat\rho)\approx \frac{m}{2n}\big[1-\rho^2(f_j,f_k)\big]^2 .
```

The last variance is *identical in form* to the coherence variance {eq}`eq-fft-cohvar` because the
bispectrum estimate at $(f_j,\,f-f_j)$ is precisely the **cross spectrum** of $\{X_t\}$ with the product
process $Z_t(f)=\sum_s X_t X_{t-s}e^{2\pi i s f}$ — whose $j$-th Fourier coefficient is
$A(f_j)A(f-f_j)$. So the entire second-order estimation theory — averaging for consistency, the
resolution–variance trade-off, coherence as normalized magnitude — carries over verbatim from the
spectrum to the bispectrum, one Fourier transform higher.

For economic data the payoff is concrete: asset returns, output growth, and exchange-rate changes are
often close to serially uncorrelated (flat spectrum) yet visibly *dependent* through volatility
clustering and asymmetries. A flat spectrum says nothing about such structure; a non-zero estimated
bicoherence pins down quadratic dependence and non-Gaussianity that the spectrum cannot see — the same
message, in the frequency domain, as the {doc}`nonlinear moving-average theory <39_nonlinear_representation>`
and the moving-spectrum diagnostics of {doc}`complex demodulation <41_comp_demod>`.

## Exercises

**1.** *(Why the periodogram is inconsistent, and what averaging buys.)* Let $\{X_t\}$ be Gaussian white
noise with mean zero and variance $\sigma^2$, so $S(f)=\sigma^2$ (flat), and take $\Delta=1$.

&nbsp;&nbsp;**A.** Using {eq}`eq-fft-Eperio` and {eq}`eq-fft-varperio`, state $E(I_k)$ and
$\operatorname{var}(I_k)$ for $k\ge1$, and explain in one sentence why the periodogram is *not* a
consistent estimator of $S$.

&nbsp;&nbsp;**B.** Form the averaged periodogram $S_k=\tfrac1r\sum_{p=1}^r I_k^{(p)}$ over $r$ independent
length-$m$ segments ($n=rm$). Using that $2 I_k^{(p)}/S(f_k)\sim\chi^2_2$ and the segments are independent,
find the distribution of $2rS_k/S(f_k)$ and hence $E(S_k)$ and $\operatorname{var}(S_k)$.

&nbsp;&nbsp;**C.** With the total sample size $n$ fixed, describe the trade-off you face in choosing the
segment length $m$: what happens to the frequency resolution and to $\operatorname{var}(S_k)$ as $m$
increases? Relate the trade-off to the uncertainty principle of {doc}`05a_uncertainty_principle`.

```{admonition} Solution to Exercise 1
:class: dropdown

**A.** With $S(f)=\sigma^2$ constant, {eq}`eq-fft-Eperio` gives $E(I_k)\approx\sigma^2$ (the Fejér-kernel
convolution of a flat spectrum is flat, so there is no leakage bias here), and {eq}`eq-fft-varperio` gives
$\operatorname{var}(I_k)\approx S(f_k)^2=\sigma^4$. The variance does **not** depend on $n$ and does not go
to zero as $n\to\infty$, so $I_k$ never settles down to $S$: it is asymptotically unbiased but
inconsistent. (Its standard deviation, $\sigma^2$, equals its mean.)

**B.** For each segment $2 I_k^{(p)}/\sigma^2\sim\chi^2_2$, and the $r$ segments are independent, so the
sum of the $r$ chi-squares is $\chi^2_{2r}$:

$$
\frac{2r\,S_k}{\sigma^2}=\sum_{p=1}^r \frac{2 I_k^{(p)}}{\sigma^2}\ \sim\ \chi^2_{2r}.
$$

Since $E\chi^2_{2r}=2r$ and $\operatorname{var}\chi^2_{2r}=4r$,

$$
E(S_k)=\frac{\sigma^2}{2r}\,E\chi^2_{2r}=\sigma^2,
\qquad
\operatorname{var}(S_k)=\Big(\frac{\sigma^2}{2r}\Big)^2\operatorname{var}\chi^2_{2r}=\frac{\sigma^4}{r}
   =\frac{m}{n}\,\sigma^4 ,
$$

confirming {eq}`eq-fft-bartvar`. Averaging $r$ independent ordinates cuts the variance by the factor $r$
and gives an estimator with $2r$ degrees of freedom, which *is* consistent as $r\to\infty$.

**C.** With $n=rm$ fixed, increasing $m$ **improves resolution** (the Fourier frequencies $f_k=k/m$ are
spaced $1/m$ apart, finer) but **reduces the number of segments** $r=n/m$, so
$\operatorname{var}(S_k)=\sigma^4/r=(m/n)\sigma^4$ **rises**. You cannot make both the bandwidth $1/m$ and
the variance $m/n$ small at once — sharpening frequency resolution costs precision. This is the
time–frequency uncertainty principle of {doc}`05a_uncertainty_principle`: a window narrow in frequency is
long in time, leaving fewer independent stretches to average.
```

**2.** *(The bispectrum sees what the spectrum cannot.)* Let $\{\varepsilon_t\}$ be i.i.d. with mean zero,
variance $\sigma^2$, symmetric ($E\varepsilon_t^3=0$), and define the nonlinear moving average

$$
X_t = \varepsilon_t + a\,\varepsilon_{t-1}\varepsilon_{t-2}, \qquad a\neq 0 .
$$

&nbsp;&nbsp;**A.** Show that $\{X_t\}$ is **white** — its autocovariances $\rho_k=E(X_tX_{t+k})$ vanish for
all $k\neq0$ — so its spectrum is flat and reveals nothing about the nonlinearity.

&nbsp;&nbsp;**B.** Compute the third-order cumulant $C(\sigma,\tau)=E(X_tX_{t+\sigma}X_{t+\tau})$ and
exhibit a lag pair $(\sigma,\tau)$ at which it is non-zero. Conclude that the bispectrum
{eq}`eq-fft-bispec` is not identically zero.

&nbsp;&nbsp;**C.** What does a non-zero estimated bicoherence $\hat\rho$ tell you about a process, and how
does this exercise motivate looking beyond the spectrum for economic series such as asset returns?

```{admonition} Solution to Exercise 2
:class: dropdown

**A.** $EX_t=E\varepsilon_t+a\,E\varepsilon_{t-1}E\varepsilon_{t-2}=0$. For $k\neq0$,

$$
\rho_k=E\big[(\varepsilon_t+a\varepsilon_{t-1}\varepsilon_{t-2})(\varepsilon_{t+k}+a\varepsilon_{t+k-1}\varepsilon_{t+k-2})\big].
$$

The four terms vanish: $E\varepsilon_t\varepsilon_{t+k}=0$; the two cross terms are third moments of
independent zero-mean variables (e.g. $aE[\varepsilon_t^2\varepsilon_{t+k-1}]... $ contains a lone factor
with zero mean), hence $0$; and the quartic term
$a^2E[\varepsilon_{t-1}\varepsilon_{t-2}\varepsilon_{t+k-1}\varepsilon_{t+k-2}]$ is non-zero only when the
index sets $\{t-1,t-2\}$ and $\{t+k-1,t+k-2\}$ coincide, i.e. $k=0$. So $\rho_k=0$ for all $k\neq0$ and
$X_t$ is white with $\rho_0=\sigma^2+a^2\sigma^4$. Its spectrum $S(f)=\sigma^2+a^2\sigma^4$ is **flat** —
the quadratic dependence is invisible to second-order analysis.

**B.** The cumulant is non-zero when the linear parts of two factors pair with the quadratic part of the
third. Take $(\sigma,\tau)=(1,2)$: in $E(X_t X_{t+1}X_{t+2})$ the only surviving contribution is
$\varepsilon_t\cdot\varepsilon_{t+1}\cdot a\varepsilon_{t+1}\varepsilon_t$ (the quadratic term of $X_{t+2}=\varepsilon_{t+2}+a\varepsilon_{t+1}\varepsilon_t$), giving

$$
C(1,2)=E(X_tX_{t+1}X_{t+2})=a\,E[\varepsilon_t^2\,\varepsilon_{t+1}^2]=a\,\sigma^4\neq0
$$

(with symmetric images at the permutations of $(1,2)$). Hence the bispectrum
$B(f_1,f_2)=\sum_{\sigma,\tau}C(\sigma,\tau)e^{-2\pi i(f_1\sigma+f_2\tau)}$ contains the term
$a\sigma^4 e^{-2\pi i(f_1+2f_2)}+(\text{permutations})$ and is **not identically zero** whenever $a\neq0$.
The bispectrum detects exactly the quadratic structure that the flat spectrum in part A missed.

**C.** Because a Gaussian or a linear-in-Gaussian process has $B\equiv0$ (so $\rho\equiv0$), a
statistically non-zero estimated bicoherence $\hat\rho$ is evidence of **non-Gaussianity and
nonlinearity**. Many economic series — asset returns, output growth, exchange-rate changes — are close to
serially uncorrelated (nearly flat spectra) yet plainly dependent through volatility clustering and
asymmetric responses. As in part A, the spectrum cannot register such dependence, but the bicoherence
can; this is why the FFT bispectrum, and the broader {doc}`nonlinear theory <39_nonlinear_representation>`,
are worth the extra Fourier transform.
```

## References

```{bibliography}
:labelprefix: FF
:filter: key in {"blackmantukey1959measurement", "cooleytukey1965algorithm", "hinichclay1968application", "jenkins1963cross", "rosenblattvanness1965estimation", "welch1967use"}
```
