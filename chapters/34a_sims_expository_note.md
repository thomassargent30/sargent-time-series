# Sims's Formula Again

The preceding section derived Sims's formula for the discrete time approximation to a continuous time distributed lag. This section presents a self-contained, purely time-domain derivation of the same formula. The derivation avoids the frequency domain calculations that Sims used. This might be useful for  economists who prefer  calculations in the time domain to ones in the frequency domain. An ancillary  benefit of the time domain derivation is that it makes apparent the interpretation of Sims's formula as a version of Theil's specification error theorem.

## The distributed lag model

Sims considered the model

```{math}
:label: eq-note-1
y(t) = \int_{-\infty}^\infty b(s)\, x(t-s)\, ds + u(t)
```

where $b(t)$ is an absolutely integrable function,[^note-1] and where $y(s)$, $x(s)$, and $u(s)$ are continuous time covariance stationary stochastic processes with means of zero and finite variances. The disturbance process $u(s)$ is orthogonal to the $x(s)$ process, that is

```{math}
:label: eq-note-2
E[u(s)\, x(t)] = 0 \qquad \text{for all real } t, s.
```

The specification {eq}`eq-note-2` identifies the convolution $\int_{-\infty}^\infty b(s)\, x(t-s)\, ds$, which is denoted $b * x(t)$, as the projection of $y(t)$ on the entire $x(s)$ process, $s \in (-\infty, \infty)$.[^note-2]

Sims studied the situation where $y(t)$ and $x(t)$ are only observed at the integers. The data thus consist of the sequences of random variables

$$
X(n) = x(n), \qquad Y(n) = y(n), \qquad n = \ldots, -2, -1, 0, 1, 2, \ldots
$$

Here $X(n)$ and $Y(n)$ are discrete time stochastic processes. Consider the least squares distributed lag regression (i.e., the projection) of $Y(n)$ on past, present, and future $X(r)$'s:

```{math}
:label: eq-note-3
Y(n) = \sum_{s=-\infty}^\infty B(s)\, X(n-s) + U(n)
```

where $U(n)$ is the least squares disturbance and where

```{math}
:label: eq-note-4
E[U(n)\, X(r)] = 0 \qquad \text{for all integers } n, r.
```

Sims was interested in studying the relation between $B(s)$ and $b(t)$ and in investigating the conditions under which $B(s)$ well or poorly represents $b(t)$ sampled at the integers.

## Deriving Sims's formula

To obtain Sims's formula, we begin by recalling that the orthogonality condition {eq}`eq-note-4` uniquely determines the $B(s)$'s.[^note-3] Substitute for $U(n)$ from {eq}`eq-note-3` into {eq}`eq-note-4` to obtain

```{math}
:label: eq-note-5
E\!\left[\left(Y(n) - \sum_{s=-\infty}^\infty B(s)\, X(n-s)\right) X(r)\right] = 0 \qquad \text{for all integers } n, r,
```

or

```{math}
:label: eq-note-6
E[Y(n)\, X(r)] - \sum_{s=-\infty}^\infty B(s)\, E[X(n-s)\, X(r)] = 0.
```

These are the least squares "normal equations." We write {eq}`eq-note-6` as

$$
R_{YX}(n-r) - \sum_{s=-\infty}^\infty B(s)\, R_X(n-r-s) = 0
$$

or

```{math}
:label: eq-note-7
R_{YX}(\tau) - \sum_{s=-\infty}^\infty B(s)\, R_X(\tau - s) = 0
```

where we define the covariance sequences

```{math}
:label: eq-note-8
\begin{aligned}
R_{YX}(\tau) &= E[Y(t)\, X(t-\tau)] && \text{for all } t, \text{ integer } \tau, \\
R_X(\tau) &= E[X(t+\tau)\, X(t)] && \text{for all } t, \text{ integer } \tau.
\end{aligned}
```

$R_{YX}$ and $R_X$ are independent of $t$ by virtue of the covariance stationarity of $Y$ and $X$. It is also convenient to define here the covariance functions

```{math}
:label: eq-note-9
\begin{aligned}
R_{yx}(\tau) &= E[y(t)\, x(t-\tau)] && \text{for all real } t, \text{ real } \tau, \\
R_x(\tau) &= E[x(t)\, x(t+\tau)] && \text{for all real } t, \text{ real } \tau.
\end{aligned}
```

Clearly, the covariance sequences $R_{YX}$ and $R_X$ correspond to the covariance functions $R_{yx}$ and $R_x$, respectively, sampled at the integers.

Now from {eq}`eq-note-1` we have that

$$
\begin{aligned}
R_{YX}(\tau) = E[Y(t)\, X(t-\tau)]
&= E\!\left[\left(\int_{-\infty}^\infty b(s)\, x(t-s)\, ds + u(t)\right) x(t-\tau)\right] \\
&= \int_{-\infty}^\infty b(s)\, E[x(t-s)\, x(t-\tau)]\, ds + E[u(t)\, x(t-\tau)]
\end{aligned}
$$

which, applying {eq}`eq-note-2` and {eq}`eq-note-9`, equals

$$
\int_{-\infty}^\infty b(s)\, R_x(\tau - s)\, ds
$$

so that

```{math}
:label: eq-note-10
R_{YX}(\tau) = b * R_x(\tau) \qquad \text{for integer values of } \tau.
```

So equation {eq}`eq-note-7` becomes

```{math}
:label: eq-note-11
b * R_x(\tau) = \sum_{s=-\infty}^\infty B(s)\, R_X(\tau - s).
```

Now let $R_X^{-1}(s)$ denote the sequence which is the inverse under convolution of the sequence $R_X(s)$. This inverse is defined by[^note-4]

```{math}
:label: eq-note-12
\sum_{j=-\infty}^\infty R_X^{-1}(j)\, R_X(n-j) =
\begin{cases}
1 & \text{for } n = 0, \\
0 & \text{for } n \neq 0, \ n \text{ an integer.}
\end{cases}
```

Convoluting the left side of {eq}`eq-note-11` with $R_X^{-1}(j)$ gives

$$
\begin{aligned}
\sum_{j=-\infty}^\infty R_X^{-1}(j)\left[b * R_x(\tau - j)\right]
&= \sum_{j=-\infty}^\infty R_X^{-1}(j) \int_{-\infty}^\infty b(s)\, R_x(\tau - j - s)\, ds \\
&= \int_{-\infty}^\infty b(s)\left(\sum_{j=-\infty}^\infty R_X^{-1}(j)\, R_x(\tau - j - s)\right) ds.
\end{aligned}
$$

Notice that

$$
\sum_{j=-\infty}^\infty R_X^{-1}(j)\, R_x(\tau - j) = R_X^{-1} * R_x(\tau),
$$

which is the convolution of the sequence $R_X^{-1}(j)$ with the function $R_x(t)$.[^note-5] Sims calls this convolution

```{math}
:label: eq-note-13
r_x(\tau) = R_X^{-1} * R_x(\tau) = \sum_{j=-\infty}^\infty R_X^{-1}(j)\, R_x(\tau - j).
```

Notice that it is defined for all real $\tau - s$. Thus, the convoluted left side can be written

$$
\int_{-\infty}^\infty b(s)\left(R_X^{-1} * R_x(\tau - s)\right) ds = b * R_X^{-1} * R_x(\tau).
$$

While this convolution exists for all real $\tau$, we are interested in its values only at integer $\tau$ (refer again to {eq}`eq-note-11`).

Now convolute the right side of equation {eq}`eq-note-11` with $R_X^{-1}(j)$ to get

$$
\begin{aligned}
R_X^{-1} * \sum_{s=-\infty}^\infty B(s)\, R_X(\tau - s)
&= \sum_{j=-\infty}^\infty R_X^{-1}(j) \sum_{s=-\infty}^\infty B(s)\, R_X(\tau - s - j) \\
&= \sum_{s=-\infty}^\infty B(s) \sum_{j=-\infty}^\infty R_X^{-1}(j)\, R_X(\tau - j - s) \\
&= B(\tau).
\end{aligned}
$$

Combining the results of convoluting the left and right sides of {eq}`eq-note-11` with $R_X^{-1}$ gives

```{math}
:label: eq-note-14
B(\tau) = b * R_X^{-1} * R_x(\tau)
```

or more explicitly

```{math}
:label: eq-note-15
B(\tau) = \int_{-\infty}^\infty b(s)\left(\sum_{j=-\infty}^\infty R_X^{-1}(j)\, R_x(\tau - j - s)\right) ds.
```

Equation {eq}`eq-note-14` is Sims's formula. It states that $B(\tau)$ is formed by weighting $b(s)$ by $R_X^{-1} * R_x(\tau - s)$ and then integrating over all real $s$'s. The weighting function $R_X^{-1} * R_x(\tau - s)$ clearly depends on the stochastic structure of the $x$-process. Sims's paper contains a variety of interesting and useful results about the shape of the weighting function for various classes of $x$-processes.

## Relation to Theil's specification theorem

Let the least squares projection of a random variable $z$ on a $1 \times k$ vector of random variables $Z$ be $Z\alpha$ where $\alpha$ is the $k \times 1$ vector of regression coefficients. Partition $Z$ as $Z = (Z_1 \ Z_2)$ where $Z_1$ is $(1 \times k_1)$ and $Z_2$ is $(1 \times k_2)$ with $k_1 + k_2 = k$. Theil's specification theorem states that the projection of $z$ on $Z_1$ is $Z_1 \xi$ where $\xi$ is $(k_1 \times 1)$ and

$$
\underset{k_1 \times 1}{\xi} = \underset{k_1 \times k}{\Gamma}\ \underset{k \times 1}{\alpha}
$$

where the projection of the $i$-th element of $Z$ on $Z_1$ is $\sum_{j=1}^{k_1} \Gamma_{ji} Z_j$; $\Gamma_{ji}$ is the partial regression coefficient of the $i$-th dependent variable on the $j$-th $Z$. The preceding equation can be written

```{math}
:label: eq-note-16
\xi_i = \sum_{j=1}^k \Gamma_{ij}\, \alpha_j,
```

which says that the coefficient on $Z_i$ in the projection of $z$ on $Z_1$ equals the vector product of $\alpha$ with the vector of $i$-th partial regression coefficients in the regressions of all of the variables in $Z$ on the subset $Z_1$.

Now consider the projection of $x(t)$ on the sampled $x$-process, $X(n)$,

$$
\sum_{j=-\infty}^\infty \gamma_j^t\, X(j),
$$

where there is one such projection, and hence one sequence $\gamma_j^t$, for each real $t$. The regression coefficients $\gamma_j^t$ are uniquely determined by the orthogonality requirement

$$
E\!\left[\left(x(t) - \sum_{j=-\infty}^\infty \gamma_j^t\, X(j)\right) X(n)\right] = 0
$$

or

$$
R_x(n-t) = \sum_{j=-\infty}^\infty \gamma_j^t\, R_X(n-j).
$$

Convoluting both sides of the above equation with $R_X^{-1}$ gives

```{math}
:label: eq-note-17
\gamma_n^t = \sum_{j=-\infty}^\infty R_X^{-1}(j)\, R_x(n-t-j) = R_X^{-1} * R_x(n-t).
```

In equation {eq}`eq-note-17`, $\gamma_n^t$ gives the regression coefficient on $X(n)$ in the projection of $x(t)$ on the $X(r)$ sequence.

Theil's specification formula {eq}`eq-note-16` leads us to expect that

```{math}
:label: eq-note-18
B(n) = \int_{-\infty}^\infty b(t)\, \gamma_n^t\, dt.
```

Substituting {eq}`eq-note-17` into the above equation convinces us that the above equation is equivalent with Sims's formula,[^note-6]

$$
B(n) = \int_{-\infty}^\infty b(t)\left(R_X^{-1} * R_x(n-t)\right) dt = b * R_X^{-1} * R_x(n).
$$

## Geweke's formula

The preceding approach can be used to derive Geweke's generalization of Sims's formula, where in {eq}`eq-note-1` $b(s)$ is now interpreted as a $1 \times k$ vector at each $s$, $x(t)$ is a $(k \times 1)$ vector stochastic process, $X(n)$ is the $(k \times 1)$ vector discrete time stochastic process corresponding to $x(t)$ sampled at the integers, and the orthogonality condition {eq}`eq-note-2` is modified to be

$$
\underset{1 \times 1}{E[u(s)\, x'(t)]} = 0 \qquad \text{for all real } t, s,
$$

where $u(s)$ is $1 \times 1$ and $x'(t)$ is $1 \times k$. Now $B(s)$ is $1 \times k$ at each integer $s$. The orthogonality condition implies

$$
E\!\left[\left(Y(n) - \sum_{s=-\infty}^\infty B(s)\, X(n-s)\right) X'(r)\right] = 0, \qquad n, s, r \text{ integers},
$$

which implies the normal equations

```{math}
:label: eq-note-19
R_{YX}(n-r) = \sum_{s=-\infty}^\infty B(s)\, R_X(n-r-s)
```

where

$$
\underset{1 \times k}{R_{YX}(\tau)} = E[Y(t)\, X'(t-\tau)], \qquad
\underset{k \times k}{R_X(\tau)} = E[X(t)\, X'(t-\tau)].
$$

Let $R_X^{-1}(j)$ be the $(k \times k)$ inverse under convolution of the matrix $R_X$, where $R_X^{-1}(j)$ is the sequence defined by

```{math}
:label: eq-note-20
\sum_{j=-\infty}^\infty R_X^{-1}(j)\, R_X(s-j) =
\begin{cases}
I_{k \times k} & s = 0, \\
0 & s \neq 0.
\end{cases}
```

Convoluting both sides of {eq}`eq-note-19` with $R_X^{-1}$ then gives equation {eq}`eq-note-14`, only where all quantities are now interpreted as the matrices defined above,

```{math}
:label: eq-note-21
B(\tau) = b * R_X^{-1} * R_x(\tau).
```

This is Geweke's formula.

## References

```{bibliography}
:labelprefix: SN
:filter: key in {"ash1972real", "geweke1975employment", "nerlove1967distributed", "papoulis1962fourier", "sims1971approximate", "sims1971discrete", "theil1971principles"}
```

[^note-1]: The results do not actually require that $b(t)$ be absolutely integrable. They will remain true if $b(t)$ is viewed as a generalized function, for example, a train of delta functions or derivatives of delta functions. See {cite:t}`sims1971discrete`. For an introductory discussion of the properties of delta functions and other generalized functions, see {cite:t}`papoulis1962fourier`.
[^note-2]: That condition {eq}`eq-note-2` uniquely determines the projection is proved, for example, by Ash (1972, p. 121).
[^note-3]: Again, see Ash (1972, p. 121).
[^note-4]: Readers familiar with lag operators may find the following helpful. The covariance generating function or $z$-transform of $R_X$ is $\rho_X(z) = \sum_{\tau=-\infty}^\infty R_X(\tau)\, z^\tau$, so that the coefficient on $z^\tau$ is the covariance at lag $\tau$. The $z$-transform of the inverse under convolution of $R_X$, namely $\rho_X^{-1}(z)$, must satisfy $\rho_X(z)\, \rho_X^{-1}(z) = 1$. Suppose, for example, that $x_t$ follows the moving average process $x_t = B(L)\varepsilon_t$, with $\varepsilon_t$ white noise of variance $\sigma_\varepsilon^2$ and $B(L) = (1 - b_1 L - b_2 L^2 - \cdots - b_p L^p)$, where $L$ is the lag operator, $L^n x_t = x_{t-n}$. It is easy to show that $\rho_X(z) = \sigma_\varepsilon^2\, B(z)\, B(z^{-1})$ (e.g., see Nerlove 1967). Then the inverse under convolution of $R_X$ has $z$-transform $\rho_X^{-1}(z) = \frac{1}{\sigma_\varepsilon^2}\frac{1}{B(z)\, B(z^{-1})}$. For example, suppose $B(L) = (1 - b_1 L)^{-1}$, so that $x$ is first-order Markov. Then $\rho_X(z) = \sigma_\varepsilon^2\frac{1}{1 - b_1 z}\frac{1}{1 - b_1 z^{-1}}$ and $\rho_X^{-1}(z) = \frac{1}{\sigma_\varepsilon^2}(1 - b_1 z)(1 - b_1 z^{-1}) = \frac{1}{\sigma_\varepsilon^2}\big({-b_1} z^{-1} + (1 + b_1^2) - b_1 z\big)$. The value of $R_X^{-1}(n)$ is the coefficient on $z^n$ in the above expression.
[^note-5]: Although $R_X^{-1} * R_x(\tau)$ is well defined in the preceding equation of the text, naming it the convolution of a sequence with a function is a slight abuse. More precisely, $R_X^{-1} * R_x(\tau)$ is the convolution of the generalized function, say $R_{Xg}^{-1}$, corresponding to the sequence $R_X^{-1}$ with $R_x$. That is, define $R_{Xg}^{-1}(t) = \sum_{n=-\infty}^\infty R_X^{-1}(n)\, \delta(t-n)$, so that $R_{Xg}^{-1}(t)$ is the generalized function with "mass" $R_X^{-1}(n)$ at integer $n$ and value zero elsewhere. Then $R_{Xg}^{-1} * R_x(t) = \int_{-\infty}^\infty R_x(t-\tau) \sum_{n} R_X^{-1}(n)\, \delta(\tau - n)\, d\tau = \sum_{n=-\infty}^\infty R_X^{-1}(n)\, R_x(t-n)$, which agrees with the definition in the text. In line with the pedagogical purpose of this note, generalized functions have been kept out of the text.
[^note-6]: Notice that for $t$ an integer we must have $\gamma_n^t = 1$ for $n = t$, and $\gamma_n^t = 0$ for all $n \neq t$. This follows because the projection of $x(t)$ on the sequence $X(s)$ is simply $X(t) = X(n)$ for $t = n$. Since $\gamma_0^t = R_X^{-1} * R_x(-t)$, this shows that the weighting function $R_X^{-1} * R_x(t)$ must be unity at $t = 0$ and zero at all other integers.
