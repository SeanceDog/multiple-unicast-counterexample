# The 182-vertex causal construction

This document defines an undirected graph with 182 vertices and 735 physical edges, together with 157 independent unicast messages over $\mathbb F_9$. A deterministic causal linear code delivers all messages using exactly one field symbol on each physical edge. Every source–sink distance is five, so shortest-path routing uses 785 field-symbol transmissions.

All finite choices are specified in [data/182-vertex-code.json](data/182-vertex-code.json). The [verifier](verify.py) reconstructs the graph and the code from those choices and checks correctness symbolically for every possible input. See the [README](README.md) for commands. The generated [transcript](data/reference-transcript.json) lists every encoder and decoder.

The communication alphabet and the unit of transmission in this construction are field symbols over $\mathbb F_9$. Vertices may hold several independent source messages or be destinations for several messages. Each demand has distinct source and destination vertices. Edges are undirected, have unit transmission cost, and have no loops or parallel copies.

## 1. Field elements, vertices, and candidate edges

Work in

$$
\mathbb F=\mathbb F_9=\mathbb F_3[\theta]/(\theta^2+1).
$$

Thus $\theta^2=-1$, and

$$
(a+b\theta)(c+d\theta)=(ac-bd)+(ad+bc)\theta,
$$

with coefficient arithmetic modulo three. In every data file, integer $r\in\{0,\ldots,8\}$ represents

$$
\eta_r=(r\bmod3)+\lfloor r/3\rfloor\theta.
$$

The integers are field-element labels; arithmetic is not arithmetic modulo nine. For example, the label $3$ represents $\theta$, and its square has label $2$.

Use these 91 triples in precisely this order:

$$
v_0=(0,0,1),\qquad
v_{1+r}=(0,1,\eta_r),\qquad
v_{10+9r+s}=(1,\eta_r,\eta_s)
\quad(0\le r,s\le8).
$$

Their set $T$ contains one representative for each one-dimensional subspace of $\mathbb F^3$: the first nonzero coordinate of each representative is one. The vertices consist of two disjoint copies of $T$. Vertex $i$ is the point vertex with representative $v_i$; vertex $91+j$ is the line vertex with representative $v_j$.

Define the candidate graph $H$ by

$$
\{i,91+j\}\in E(H)\quad\Longleftrightarrow\quad v_i\cdot v_j=0.
$$

For each point index $i$, list its line indices increasingly as

$$
j_{i,0}<j_{i,1}<\cdots<j_{i,9}.
$$

The zero-based index of the edge $\{i,91+j_{i,r}\}$ is $10i+r$. Hence candidate edge indices range from 0 to 909. We write $e_q$ for the edge with index $q$.

To verify the counts, the vectors orthogonal to any nonzero triple form a two-dimensional vector space. It has $9^2-1=80$ nonzero vectors, with eight scalar multiples per projective representative. Every vertex therefore has ten neighbors. Consequently,

$$
|V(H)|=182,\qquad |E(H)|=910.
$$

The graph is bipartite. Two distinct point representatives are linearly independent, so their common orthogonal space has dimension one and gives exactly one common line neighbor. In particular, $H$ has no four-cycles.

## 2. The finite order, deleted edges, and demands

The parameter file has four keys:

| Key | Meaning |
| --- | --- |
| `order` | A permutation $\pi=(\pi_0,\ldots,\pi_{181})$ of the vertex IDs. |
| `backward` | The set $R$ of 175 deleted candidate-edge indices. |
| `C` | Nine deleted coordinates whose values are publicly fixed at zero. |
| `J` | Nine other deleted coordinates, carrying no initial information. |

All indices start at zero. Write $\operatorname{pos}(\pi_a)=a$. An endpoint is *earlier* if its position is smaller. The array $\pi$ lists vertices in processing order; it is not the position function.

The two small sets are

$$
\begin{aligned}
C&=\{102,163,164,341,386,406,476,513,862\},\\
J&=\{103,349,470,517,603,606,674,716,861\}.
\end{aligned}
$$

They are disjoint subsets of $R$. Set

$$
B=R\setminus(C\cup J),\qquad |B|=157.
$$

The physical communication graph is

$$
G=H\setminus\{e_q:q\in R\}.
$$

It has 735 edges. None of the edges indexed by $B$, $C$, or $J$ is a communication edge.

For use in the proof, give every candidate edge an auxiliary direction: edges outside $R$ point from the earlier endpoint to the later endpoint; edges in $R$ point from the later endpoint to the earlier endpoint. Thus the latter edges point backward in the processing order. The finite parameters satisfy

$$
d^-(v)=d^+(v)=5,
\qquad d_C^+(v)=d_J^-(v)
\quad\text{for every vertex }v.
$$

Here all degrees refer to the auxiliary orientation, and subscripts restrict the counted edges. The verifier checks these identities directly. The second identity can also be read from the following complete table; its two sides are zero at every unlisted vertex.

| Vertex | Outgoing $C$ edge | Incoming $J$ edge |
| ---: | ---: | ---: |
| 10 | 102 | 103 |
| 34 | 341 | 349 |
| 47 | 476 | 470 |
| 51 | 513 | 517 |
| 86 | 862 | 861 |
| 125 | 163 | 603 |
| 134 | 164 | 674 |
| 147 | 406 | 716 |
| 149 | 386 | 606 |

For each $q\in B$, place an arbitrary independent input $x_q\in\mathbb F$ at the earlier endpoint $s_q$ of $e_q$, and demand it at the later endpoint $t_q$. Thus $s_q$ is the auxiliary head and $t_q$ the auxiliary tail of this backward edge. Integer message IDs $0,\ldots,156$ correspond to the increasing enumeration of $B$.

## 3. The local interpolation rules

At a vertex represented by $z\in T$, its neighboring triples lie in

$$
z^\perp=\{y\in\mathbb F^3:z\cdot y=0\}.
$$

Choose the first two triples of $T\cap z^\perp$, in the ordering of Section 1, as a basis $u,w$. They are distinct normalized representatives and hence linearly independent. Express each neighboring normalized triple uniquely as

$$
y=\lambda u+\mu w.
$$

Use these exact coordinates. Do not normalize $(\lambda,\mu)$ again.

At a point vertex, allowed incident values are the evaluations of

$$
Q(\lambda,\mu)=\sum_{r=0}^4a_r\lambda^{4-r}\mu^r.
$$

At a line vertex, allowed incident values are the evaluations of

$$
Q(\lambda^3,\mu^3)=\sum_{r=0}^4a_r\lambda^{3(4-r)}\mu^{3r}.
$$

The coefficients $a_0,\ldots,a_4$ may differ at different vertices. Any five distinct incident evaluations uniquely determine them. Indeed, a nonzero homogeneous binary quartic has at most four distinct projective zeros. Cubing is an automorphism of $\mathbb F$, so it preserves distinctness of the projective coordinate pairs.

Let $W\subseteq\mathbb F^{910}$ consist of assignments to all candidate edges satisfying both endpoint rules simultaneously. Section 7 gives an elementary proof that

$$
\dim W\ge166.
$$

## 4. Why every input is encodable

Set

$$
W_C=\{w\in W:w_q=0\text{ for every }q\in C\}.
$$

Since there are only nine imposed coordinate equations,

$$
\dim W_C\ge166-9=157.
$$

At vertex $v$, define the set $S_v$ of known ports to consist of its incoming physical edges, incoming $B$ edges, and all incident $C$ edges. These are distinct ports. Their number is

$$
\begin{aligned}
|S_v|
&=d^-_{E(G)}(v)+d^-_B(v)+d^-_C(v)+d^+_C(v)\\
&=5-d^-_J(v)+d^+_C(v)=5.
\end{aligned}
$$

Restriction to $B$ defines a linear map

$$
\rho:W_C\longrightarrow\mathbb F^B.
$$

It is injective. Suppose $w\in W_C$ has all $B$ coordinates zero, and process vertices in order. At each vertex, its five ports in $S_v$ are zero: incoming physical-edge values have already been proved zero at earlier vertices, while the $B$ and $C$ values are zero by assumption. Five zero evaluations force the entire local polynomial to vanish. Induction proves every candidate-edge value zero, including the unused $J$ coordinates.

Both the codomain dimension and the lower bound for the domain dimension are 157. Therefore $\rho$ is an isomorphism. Every arbitrary input tuple $(x_q)_{q\in B}$ has a unique compatible extension in $W_C$. This does not impose a relation among the source messages.

## 5. The causal protocol and its cost

Process vertices in the order $\pi$. When processing $v$:

1. Read the five values at the ports in $S_v$. Incoming physical-edge packets have already arrived. Incoming $B$ values are source messages initially held at $v$. Incident $C$ values are the public constant zero.
2. Interpolate the unique local polynomial using the rule of Section 3.
3. Send its evaluation along every outgoing physical edge, once, in increasing edge-index order.
4. For every $q\in B$ whose destination is $v$, output the polynomial's evaluation at port $q$.

The unique compatible extension from Section 4 proves correctness by induction: the five supplied values are its true values, so interpolation recovers its other incident values. Each destination therefore outputs $x_q$. Every operation is linear over $\mathbb F$, and every encoder uses only the sender's own inputs and packets received earlier. The compatible extension is used in the proof; computing it globally is not an operation of the protocol.

Each physical edge is used exactly once, from earlier to later. Hence

$$
\operatorname{cost}_{\mathrm{code}}=735.
$$

For a demand $q\in B$, its endpoints are incident in $H$ but their edge is deleted from $G$. A path between them in the bipartite graph $G$ has odd length. Its length cannot be one. If its length were three, adjoining $e_q$ would form a four-cycle in $H$, which is impossible. Thus every demand distance is at least five. Breadth-first search in the explicitly defined graph verifies that all 157 distances equal five. Consequently,

$$
\operatorname{cost}_{\mathrm{routing}}=157\cdot5=785>735.
$$

The code saves exactly 50 field-symbol transmissions relative to shortest-path routing.

### Throughput compared with fractional flow

Give each undirected edge total capacity one field symbol per time slot, shared between its two directions. No additional constraint is imposed on simultaneous transmissions on different edges. Let $R_{\mathrm{code}}$ and $R_{\mathrm{flow}}$ be the supremal common rates per session for causal network coding and fractional multicommodity flow.

Pipeline independent generations of the code. For generation $b$, process vertex $v$ at time $b+\operatorname{pos}(v)$. A packet from an earlier vertex arrives before the receiving vertex processes that generation, using unit transmission delay. Each fixed edge sends at most one symbol per time slot. The startup delay is bounded independently of the number of generations, so

$$
R_{\mathrm{code}}\ge1.
$$

For fractional routing at common rate $r$, let $f_{q,P}$ be the flow of session $q$ on path $P$. Every such path has at least five edges. Summing the edge capacities therefore gives

$$
735\ge\sum_{q,P}|P|f_{q,P}
\ge5\sum_{q,P}f_{q,P}=785r.
$$

This permits arbitrary splitting among paths and both directions of an edge, subject to its shared capacity. It follows that

$$
R_{\mathrm{flow}}\le\frac{147}{157}<1\le R_{\mathrm{code}}.
$$

The number $147/157$ is an upper bound on flow throughput; no claim that this upper bound is attained is needed.

## 6. What the symbolic verifier proves

The verifier provides a second, direct correctness proof for these finite parameters. It does not assume the dimension bound of Section 7 and does not rely on sampled inputs or on the saved transcript.

It first reconstructs field arithmetic, all 91 representatives, all 910 candidate edges, the auxiliary orientation, the 735 physical edges, the 157 terminal pairs, and the local evaluation matrices. It checks the parameter counts, balanced degrees, the $C$–$J$ identity, and the five-port condition.

It represents input $x_i$ by the $i$th standard basis vector in $\mathbb F^{157}$. Every packet is then a coefficient vector for a linear form in the 157 independent variables. The verifier runs the actual chronological interpolation procedure on these vectors. For every nonzero term of every encoder and decoder, it checks that the referenced input belongs to that vertex or that the referenced packet was previously received at that vertex. It checks each decoder vector equals exactly the required standard basis vector. Therefore every decoder is correct for all $9^{157}$ possible input tuples.

It also checks that the locally reconstructed values agree at the two endpoints of every one of the 910 candidate edges, including $J$, and computes the shortest-path distances in the physical graph. The transcript records each local encoder, its global coefficient vector, each decoder, and all demand distances. It is an inspectable output of verification, not an assumed source of correct coefficients.

## 7. An elementary analytic proof that $\dim W\ge166$

This section proves the dimension bound without finite-matrix rank calculations or algebraic geometry.

### Compatible assignments from polynomials

Let $X=(X_0,X_1,X_2)$ and $Y=(Y_0,Y_1,Y_2)$. A polynomial $h(X,Y)$ is bihomogeneous of bidegree $(4,4)$ if each of its monomials has total degree four in $X$ and total degree four in $Y$. For every such polynomial, define the edge assignment

$$
w_h(p,\ell)=h(p^{(3)},\ell),
\qquad p^{(3)}=(p_0^3,p_1^3,p_2^3).
$$

This assignment belongs to $W$. At a point vertex $p$, it is a quartic in the coordinates of $\ell$. At a line vertex, write $p=\lambda u+\mu w$. Characteristic three gives

$$
p^{(3)}=\lambda^3u^{(3)}+\mu^3w^{(3)},
$$

so the value is a quartic in $(\lambda^3,\mu^3)$, as required.

There are $\binom{6}{2}^2=225$ monomials $X^aY^b$ with $|a|=|b|=4$. Retain those divisible by neither

$$
X_0Y_0^3\qquad\text{nor}\qquad X_1^3Y_1.
$$

The first divisibility condition excludes $10\cdot3=30$ monomials, and the second excludes $3\cdot10=30$. Their intersection consists only of $X_0X_1^3Y_0^3Y_1$. Thus there are

$$
225-30-30+1=166
$$

selected monomials. We prove that their edge assignments are linearly independent.

### From edge evaluations to a polynomial identity

Let $h$ be a linear combination of the selected monomials, and suppose $w_h$ is zero on every candidate edge. Work temporarily over an algebraic closure $K$ of $\mathbb F$.

For each normalized $x\in T$, take $p=x^{(3)}$. Its normalization is unchanged, and $p^{(3)}=x$ because every element of $\mathbb F_9$ satisfies $a^9=a$. Hence $h(x,Y)$ vanishes at the ten rational projective points of the line

$$
C_x:\quad x^{(3)}\cdot Y=0.
$$

Its restriction to this line is a homogeneous binary quartic. Having ten distinct projective zeros forces that restriction to vanish identically, over $K$ as well as over $\mathbb F$.

Introduce three formal variables $Z=(Z_0,Z_1,Z_2)$ and the cross product

$$
\begin{aligned}
A(Z)&=Z\times Z^{(9)}\\
&=(Z_1Z_2^9-Z_2Z_1^9,\;
Z_2Z_0^9-Z_0Z_2^9,\;
Z_0Z_1^9-Z_1Z_0^9).
\end{aligned}
$$

Define the formal polynomial

$$
Q_h(Z)=h(A(Z),Z^{(3)}).
$$

Its degree is at most $4\cdot10+4\cdot3=52$. Here $Z_i^9$ and $Z_i$ are different formal polynomials; no reduction modulo $Z_i^9-Z_i$ is made.

Fix $x\in T$ and consider the entire $K$-line $x\cdot Z=0$. Since $x^{(9)}=x$,

$$
x\cdot Z^{(9)}=(x\cdot Z)^9=0.
$$

Thus $Z$ and $Z^{(9)}$ lie in $x^\perp$. Their cross product is a scalar multiple $\lambda x$, including the case $\lambda=0$. Also $x^{(3)}\cdot Z^{(3)}=(x\cdot Z)^3=0$. The previously established identity on $C_x$ gives

$$
Q_h(Z)=h(\lambda x,Z^{(3)})
=\lambda^4h(x,Z^{(3)})=0.
$$

Each of the 91 distinct linear polynomials $x\cdot Z$, for $x\in T$, therefore divides $Q_h$. Their product has degree 91, whereas $Q_h$ has degree at most 52. It follows that $Q_h$ is the zero polynomial.

### Distinct leading terms imply independence

Use lexicographic monomial order $Z_0>Z_1>Z_2$. The leading terms of the components of $A$ are

$$
\operatorname{LT}(A_0)=-Z_1^9Z_2,\qquad
\operatorname{LT}(A_1)=Z_0^9Z_2,\qquad
\operatorname{LT}(A_2)=-Z_0^9Z_1.
$$

For a selected monomial $X^aY^b$, the leading exponent vector of its pullback $A(Z)^a(Z^{(3)})^b$ is

$$
\bigl(9(a_1+a_2)+3b_0,\;
9a_0+a_2+3b_1,\;
a_0+a_1+3b_2\bigr).
$$

Its leading coefficient is $(-1)^{a_0+a_2}$ and is nonzero.

Suppose two selected monomials, with exponent pairs $(a,b)$ and $(a',b')$, have the same leading exponent vector. Equality of the first coordinates and $|a|=|a'|=4$ imply

$$
b'_0-b_0=3(a'_0-a_0).
$$

If $a'_0>a_0$, then $a'_0\ge1$ and $b'_0\ge3$, contradicting the exclusion of monomials divisible by $X_0Y_0^3$. The reverse inequality gives the same contradiction for the unprimed monomial. Thus $a'_0=a_0$ and $b'_0=b_0$.

Equality of the second leading coordinates now gives

$$
a'_1-a_1=3(b'_1-b_1).
$$

If $b'_1>b_1$, then $b'_1\ge1$ and $a'_1\ge3$, contradicting the exclusion of $X_1^3Y_1$. Again the reverse inequality gives the same contradiction. Therefore $a'_1=a_1$ and $b'_1=b_1$. The fixed total degrees force $a'_2=a_2$ and $b'_2=b_2$ as well.

The 166 selected pullbacks consequently have distinct leading monomials and are linearly independent. Since $Q_h=0$, all coefficients of $h$ vanish. Their 166 compatible edge assignments are therefore linearly independent, proving

$$
\boxed{\dim W\ge166.}
$$

Together with the injectivity argument of Section 4, this also gives $\dim W_C=157$ and $\dim W=166$: the nine coordinate conditions imply $\dim W\le\dim W_C+9=166$.
