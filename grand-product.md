
## Terminology \& Notation

We consider $\mathbb{F}$ a finite field of prime order $p$, i.e., $\mathbb{F} = \mathbb{Z}_p$.
For a given integer $n$, we denote by $[n]$ the set of integers $\{1,...,n\}$. 
We explicitly define the multiplicative subgroup $H\subset\mathbb{F}$ as the subgroup containing the $n$-th roots of unity in $\mathbb{F}$, where $\omega$ is a primitive $n$-th root of unity and a generator of $H$. That is,
$$
H = \{\omega, \dots, \omega^{n-1}, \omega^n = 1\}.
$$

We denote by $\mathbb{F}_{<n}[X]$ the set of univariate polynomials over $\mathbb{F}$ of degree strictly smaller than $n$. 

For a polynomial $f(x) \in \mathbb{F}_{<n}[X]$ and $i \in [n]$, we sometimes denote $f_i := f(\omega^i)$. For a vector $\boldsymbol{f} \in \mathbb{F}^n$, we also denote by $f(x)$ the polynomial in $\mathbb{F}_{<n}[X]$ with $f(\omega^i) = f_i$.

The Lagrange polynomial $L_i(X) \in \mathbb{F}_{<n}[X]$ for $i \in [n]$ has the form
$$
L_i(X) = \frac{\omega^i\:(X^n - 1)}{n\:(X - \omega^i)}.
$$

Thus, the zero polynomial $Z_H(X) \in \mathbb{F}_{<n+1}[X]$ is defined as
$$
Z_H(X) = (X - \omega)  \cdots  (X - \omega^{n-1})(X - \omega^n) = X^n - 1.
$$
 
## (Non-Interactive) Protocol

**Common Preprocessed Input**

\begin{align}
&n, \quad \left([1]_1, [x]_1, \dots, [x^{n+4}]_1, [1]_2, [x]_2\right).
\end{align}

### Prover Algorithm

**Prover Input**: Multisets $\boldsymbol{f} = \{f_1, f_2, \dots, f_n\}$  and $\boldsymbol{t} = \{t_1, t_2, \dots, t_n\}$. The prover aims to demonstrate that these multisets are the same (or, equivalently, that there exists a permutation between their values), i.e., $\boldsymbol{f} = \boldsymbol{t}$, by reducing it to a grand-product check. This is often called a *multiset equality check*.

**Round 1:**

1. Compute the polynomials $f,t \in \mathbb{F}_{<n}[X]$ resulting from the interpolation of $\boldsymbol{f}$ and $\boldsymbol{t}$, respectively.
     
1. Compute $[f(x)]_1$, $[t(x)]_1$, the commitments to $f$ and $t$.

The first output of the prover is $\left([f(x)]_1, [t(x)]_1\right)$.

**Round 2:**

1. Compute the (reduction to grand-product check) challenge $\gamma \in \mathbb{F}_p$:
$$
\gamma = \mathsf{Hash}(\mathsf{transcript}).
$$

1. Compute the grand-product polynomial $z \in \mathbb{F}_{<n}[X]$:
\begin{align*}
z(X) = &~L_1(X) + \sum_{i=1}^{n-1} \left( L_{i+1}(X)\prod_{j=1}^{i} \frac{(f_j + \gamma)}{(t_j + \gamma)} \right).
\end{align*}

1. Compute $[z(x)]_1$.

The second output of the prover is $[z(x)]_1$.

**Round 3:**

1. Compute the (reduction to quotient check) challenge $\alpha \in \mathbb{F}_p$
$$
\alpha = \mathsf{Hash}(\mathsf{transcript}).
$$

2. Compute the quotient polynomial $q \in \mathbb{F}_{<n-1}[X]$:
\begin{align}
q(X) = 
\frac{1}{Z_H(X)}
\left(
\begin{array}{l}
(z(X)-1)L_1(X)~+ 
\\[0.1cm] +~\alpha\big[z(X\omega)(t(X)+\gamma) - z(X)(f(X)+\gamma)\big]
\end{array}
\right)
\end{align}

5. Compute $[q(x)]_1$.

The third output of the prover is $[q(x)]_1$.

**Round 4:**

1. Compute the evaluation challenge $\displaystyle \mathfrak{z} \in \mathbb{F}_p$:
$$
\mathfrak{z} = \mathsf{Hash}(\mathsf{transcript}),
$$

1. Compute the opening evaluations:
\begin{align*}
f(\mathfrak{z}), \quad z(\mathfrak{z}\omega).
\end{align*}

The fourth output of the prover is $\left(f(\mathfrak{z}), z(\mathfrak{z}\omega)\right)$.

**Round 5:**

1. Compute the opening challenges $v \in \mathbb{F}_p$:
$$
v = \mathsf{Hash}(\mathsf{transcript})
$$

1. Compute the linearisation polynomial $r \in \mathbb{F}_{<n}[X]$:
\begin{align}
r(X) :=~ &(z(X)-1)L_1(\mathfrak{z})~+ \\[0.2cm] &+\alpha\big[z(\mathfrak{z}\omega) (t(X)+\gamma) - z(X)(f(\mathfrak{z})+\gamma)\big] - Z_H(\mathfrak{z})q(X).
\end{align}

3. Compute the opening proof polynomials $W_{\mathfrak{z}},W_{\mathfrak{z}\omega} \in \mathbb{F}_{<n}[X]$:
\begin{align*}
W_{\mathfrak{z}}(X) &= \frac{1}{X - \mathfrak{z}}\left[r(
X) + v(f(X) - f(\mathfrak{z}))\right] \\[0.2cm]
W_{\mathfrak{z}\omega}(X) &= \frac{z(X) - z(\mathfrak{z}\omega)}{X - \mathfrak{z}\omega}
\end{align*}

4. Compute $[W_{\mathfrak{z}}(x)]_1,[W_{\mathfrak{z}\omega}(x)]_1$.

The fifth output of the prover is $([W_{\mathfrak{z}}(x)]_1,[W_{\mathfrak{z}\omega}(x)]_1)$.

The complete proof is:
$$\pi = 
\left(
\begin{align*}
&[f(x)]_1, [t(x)]_1, [z(x)]_1, [q(x)]_1, [W_{\mathfrak{z}}(x)]_1, [W_{\mathfrak{z}\omega}(x)]_1 \\[0.2cm]
&\qquad\qquad\qquad\qquad f(\mathfrak{z}), z(\mathfrak{z}\omega)
\end{align*}
\right)
$$

Compute the multipoint evaluation batcher $u\in \mathbb{F}_p$:
$$
u = \mathsf{Hash}(\mathsf{transcript}),
$$


### Verifier Algorithm

1. Validate $[f(x)]_1$, $[t(x)]_1$, $[z(x)]_1$, $[q(x)]_1$, $[W_\mathfrak{z}(x)]_1$, $[W_{\mathfrak{z}\omega}(x)]_1$ $\in \mathbb{G}_1$.

1. Validate that $f(\mathfrak{z}), z(\mathfrak{z}\omega) \in \mathbb{F}$.

1. Compute the challenges $\gamma, \alpha, \mathfrak{z}, v, u \in \mathbb{F}$ as in prover description, from the inputs and the elements of $\pi$.

1. Compute the zerofier evaluation $Z_H(\mathfrak{z}) = \mathfrak{z}^n-1$.

1. Compute the Lagrange polynomial evaluations $L_1(\mathfrak{z}) = \frac{\omega\:(\mathfrak{z}^n - 1)}{n\:(\mathfrak{z} - \omega)}$.

1. To save a verifier scalar multiplication, we split $r$ into its constant and non-constant terms. Compute $r$'s constant term:
\begin{align}
r_0 := \alpha \cdot \gamma \cdot z(\mathfrak{z}\omega) - L_1(\mathfrak{z}),
\end{align}
and let $r'(X) := r(X) - r_0$.

1. Compute the first part of the KZG check $\left[D\right]_1 := [r'(x)]_1 + u [z(x)]_1$:
$$[D]_1 := (L_1(\mathfrak{z}) - \alpha(f(\mathfrak{z}) + \gamma) + u)[z(x)]_1 + \alpha \cdot z(\mathfrak{z}\omega)[t(x)]_1 - Z_H(\mathfrak{z})[q(x)]_1.$$

8. Compute the full batched polynomial commitment $[F]_1$:
$$
[F]_1 := [D]_1 + v \cdot [f(x)]_1.
$$

9. Compute the group-encoded batch evaluation $[E]_1$:
$$
[E]_1 := \bigg[- r_0 + vf(\mathfrak{z}) + u \cdot z(\mathfrak{z}\omega)\bigg][1]_1
$$

10. Finally, batch validate all evaluations:
$$
e\left([W_\mathfrak{z}(x)]_1+u\cdot[W_{\mathfrak{z}\omega}(x)]_1, [x]_2\right) \stackrel{?}{=} e\left(\mathfrak{z}\cdot [W_\mathfrak{z}(x)]_1+u\mathfrak{z}\omega\cdot[W_{\mathfrak{z}\omega}(x)]_1 + [F]_1 - [E]_1, [1]_2\right).
$$
which is the same as checking that:
$$
e\left(-[W_\mathfrak{z}(x)]_1-u\cdot[W_{\mathfrak{z}\omega}(x)]_1, [x]_2\right) \cdot e\left(\mathfrak{z}\cdot [W_\mathfrak{z}(x)]_1+u\mathfrak{z}\omega\cdot[W_{\mathfrak{z}\omega}(x)]_1 + [F]_1 - [E]_1, [1]_2\right) \stackrel{?}{=} 1.
$$

###### tags: `Plookup`