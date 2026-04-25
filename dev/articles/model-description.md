# tinyVAST model description

## Bivariate generalized linear mixed model structure

tinyVAST (J. T. Thorson et al. 2025) is a bivariate extension of a
generalized linear mixed model (see Tables 1 and 2 for notation), which
includes two linear predictors that each include four additive
components:

1.  *Spatial interactions among variables*: The user can specify
    interactions among variables at a given site (or for spatially
    correlated latent variables) using arrow notation derived from path
    analysis, based on the interface from R package sem;

2.  *Temporal interaction among variables*: The user can specify
    simultaneous and lagged interactions among variables over time using
    an expanded arrow-and-lag notation that is derived from R package
    dsem, where these interactions the annual intercept for a given
    variable therefore apply uniformly for all locations.

3.  *Spatio-temporal interactions among variables*: The user can specify
    simultaneous and lagged interactions among variables over time,
    where these interactions occur on a site-by-site basis.

4.  *Generalized additive model*: The user specifies a formula for a
    generalized additive model (GAM) that is interpreted by R package
    mgcv. If other model components are missing, tinyVAST estimates
    parameters that are similar to mgcv, with small differences
    resulting from different methods for parameter estimation;

These four components are assembled into two linear predictors:

$$\begin{aligned}
p_{1,i} & {= \underset{\text{formula}}{\underbrace{\mathbf{X}_{1,i}{\mathbf{α}}_{1} + \mathbf{Z}_{1,i}{\mathbf{γ}}_{1}}} + \underset{\text{space\_term}}{\underbrace{\mathbf{A}_{i}\mathbf{\Omega}_{1,c{\lbrack i\rbrack}}}} + \underset{\text{time\_term}}{\underbrace{\mathbf{D}_{1,c{\lbrack i\rbrack},t{\lbrack i\rbrack}}}} + \underset{\text{spacetime\_term}}{\underbrace{\mathbf{A}_{i}\mathbf{E}_{1,c{\lbrack i\rbrack},t{\lbrack i\rbrack}}}} + \underset{\text{spatially\_varying}}{\underbrace{\mathbf{W}_{1,i}\left( \mathbf{A}_{i}\mathbf{\Xi}_{1} \right)^{T}}}} \\
p_{2,i} & {= \underset{\text{formula}}{\underbrace{\mathbf{X}_{2,i}{\mathbf{α}}_{2} + \mathbf{Z}_{2,i}{\mathbf{γ}}_{2}}} + \underset{\text{space\_term}}{\underbrace{\mathbf{A}_{i}\mathbf{\Omega}_{2,c{\lbrack i\rbrack}}}} + \underset{\text{time\_term}}{\underbrace{\mathbf{D}_{2,c{\lbrack i\rbrack},t{\lbrack i\rbrack}}}} + \underset{\text{spacetime\_term}}{\underbrace{\mathbf{A}_{i}\mathbf{E}_{2,c{\lbrack i\rbrack},t{\lbrack i\rbrack}}}} + \underset{\text{spatially\_varying}}{\underbrace{\mathbf{W}_{2,i}\left( \mathbf{A}_{i}\mathbf{\Xi}_{2} \right)^{T}}}}
\end{aligned}$$

where

- $p_{1,i}$ is the first linear predictor;
- $\mathbf{X}_{1,i}{\mathbf{α}}_{1} + \mathbf{Z}_{1,i}{\mathbf{γ}}_{1}$
  is the contribution of the GAM component specified using a formula
  parsed by package `mgcv`;
- ${\mathbf{α}}_{1}$ is the GAM fixed effects with associated design
  matrix $\mathbf{X}_{1}$, where $\mathbf{X}_{1,i}$ is the row of the
  design matrix for sample $i$;
- ${\mathbf{γ}}_{1}$ is the GAM random effects with associated design
  matrix $\mathbf{Z}$, where $\mathbf{Z}_{1,i}$ is the row of the design
  matrix;
- $\mathbf{A}_{i}\mathbf{\Omega}_{1,c{\lbrack i\rbrack}}$ is
  contribution of the space-veriable interaction $\mathbf{\Omega}$
  projected to sample $i$, where $\mathbf{A}$ is the interpolation
  matrix with dimension $I \times S$ such that row $\mathbf{A}_{i}$
  interpolates spatial vector $\mathbf{\Omega}_{1,c{\lbrack i\rbrack}}$
  for sample $i$ and variable $c\lbrack i\rbrack$;
- $\mathbf{D}_{1,c{\lbrack i\rbrack},t{\lbrack i\rbrack}}$ is the
  contribution of the time-variable interaction $\mathbf{D}$ for sample
  $i$;
- $\mathbf{A}_{i}\mathbf{E}_{1,c{\lbrack i\rbrack},t{\lbrack i\rbrack}}$
  is the contribution of the space-variable-time interaction
  $\mathbf{E}$ projected to sample $i$;
- $\mathbf{W}_{1,i}\mathbf{A}_{i}\mathbf{\Xi}_{1,c{\lbrack i\rbrack}}$
  is the contribution of spatially varying coefficients, i.e., a
  zero-centered and spatially varying slope
  $\Xi_{1,c{\lbrack i\rbrack}}$ that is projected to the location of
  samples, where these local slopes are then multiplied by user-supplied
  covariates in design matrix $W_{1,i}$

and terms are defined similarly for the second linear predictor
$p_{2,i}$ except using subscript-2. The linear predictors are then
passed through a bivariate inverse-link function to specify the
distribution for errors:

$$y_{i} \sim f_{e{\lbrack i\rbrack}}\left( g_{e{\lbrack i\rbrack}}^{- 1}\left( p_{1,i},p_{2,i} \right),\theta_{e{\lbrack i\rbrack}} \right)$$

where

- $f_{e{\lbrack i\rbrack}}$ probability density or mass function for
  sample $i$;
- $g_{e{\lbrack i\rbrack}}^{- 1}\left( p_{1,i},p_{2,i} \right)$ is the
  bivariate inverse-link function that transforms linear predictors to
  the central tendancy parameters of the distribution;
- $\theta_{e{\lbrack i\rbrack}}$ is the dispersion parameters for
  distribution $e\lbrack i\rbrack$ for sample $i$;

In the simple case, the distribution only requires a single linear
predictor such that $\mathbf{p}_{2} = \mathbf{0}$ by construction and
drops out of the model. In this case, the model collapses to a
generalized linear mixed model. For example we might have a log-link and
Poisson distribution, such that it collapses to a log-linked Poisson
GLMM, $y_{i} \sim {Poisson}\left( e^{p_{1,i}} \right)$.

However, tinyVAST can also handle a delta-model using either logit-log
or Poisson-linked link functions:

$$g_{e{\lbrack i\rbrack}}^{- 1}\left( p_{1,i},p_{2,i} \right) = \left( \mu_{1,i},\mu_{2,i} \right)$$

where $\mu_{1,i}$ and $\mu_{2,i}$ are the two linear predictors, such
that
${\mathbb{E}}\left\lbrack y_{i} \right\rbrack = \mu_{1,i}\mu_{2,i}$;

For the conventional logit-log bivariate link function we use a
logit-link for encounter probabilities and a log-link for positive catch
rates:

$$\begin{aligned}
\mu_{1,i} & {= \frac{e^{p_{1,i}}}{1 + e^{p_{1,i}}}} \\
\mu_{2,i} & {= e^{p_{2,i}}}
\end{aligned}$$

while for the Poisson-linked link function (J. T. Thorson 2018) we
specify a complemetary log-log link for encounter probabilities and
define the second link function such that
$\mu_{1,i}\mu_{2,i} = e^{p_{1,i}}e^{p_{2,i}}$:

$$\begin{aligned}
\mu_{1,i} & {= 1 - e^{- e^{p_{1,i}}}} \\
\mu_{2,i} & {= \frac{e^{p_{1,i}}}{\mu_{1,i}}e^{p_{2,i}}}
\end{aligned}$$ where $e^{p_{1,i}}$ is the density for an underlying
point process, and $e^{p_{2,i}}$ is biomass per point (i.e., animal or
group).

In either the conventional or Poisson-linked delta model,
$\mu_{1,i} = {Pr}(Y > 0)$ is the encounter probability (and is fitted
with a Bernoulli distribution), while $\mu_{2,i}$ is the central
tendancy for positive values, i.e.,
$\mu_{2,i} = {\mathbb{E}}\left( Y|Y > 0 \right)$.

Usefully, the interface allows analysts to specify a different
distribution $e_{i}$ for different partitions of the data. This allows
the model to jointly analyze different data types (Grüss and Thorson
2019), or data with different magnitudes of sampling error.

## Spatial domains

Linear predictors include spatially autocorrelated latent variables.
These variables are treated as Gaussian Markov random fields (GMRFs),
and evaluating the probability density of GMRFs involves calculating the
precision matrix $\mathbf{Q}_{domain} = \mathbf{\Sigma}^{- 1}$ as the
inverse of the spatial covariance matrix (J. Thorson and Kristensen
2024). tinyVAST involves three options for specifying this spatial
precision:

#### Stochastic partial differentiaul equation (SPDE)

The analyst can approximate spatial variation over a continuous
two-dimensional surface by constructing finite element mesh (FEM),
treating the value at vertices as a GMRF, and then using bilinear
interpolation (i.e., a piecewise linear approximation) to interpolate
from vertices to the spatial domain. This method was developed by (F.
Lindgren, Rue, and Lindström 2011), as popularized in the software
R-INLA (Lindgren 2012), first implemented in TMB by (J. T. Thorson et
al. 2014), and using elements constructed by R-package `fmesher` (F.
Lindgren 2023). In this case, the precision is constructed as:

$$\mathbf{Q}_{domain} = \tau^{2}\left( \kappa^{4}\mathbf{M}_{0} + 2\kappa^{2}\mathbf{M}_{1} + \mathbf{M}_{2} \right)$$
where every row $\mathbf{A}_{i}$ of the interpolation matrix
$\mathbf{A}$ is nonzero for only the three vertices of the triangle that
contains sample $i$

#### Simultaneous autoregressive (SAR):

Alternatively, the analyst can apply a simultaneous autoregressive
process (Ver Hoef, Hanks, and Hooten 2018), specifying an areal model
that represents the value within spatial strata:

$$\mathbf{Q}_{domain} = \tau^{2}\left( \mathbf{I} - \kappa\mathbf{A}^{*} \right)^{2}$$
where $\mathbf{A}^{*}$ is the adjacency matrix of the graph specified by
the analyst, $\kappa$ is the estimated partial correlation between
adjacent areas, and each row $\mathbf{A}_{i}$ of interpolation matrix
$\mathbf{A}$ is nonzero only for the single spatial stratum containing
sample $i$ (and noting that adjacency matrix $\mathbf{A}^{*}$ is
different from interpolation matrix $\mathbf{A}$, but we use the same
notation for both due to a collision in notational standards).

#### Stream networks using a tail-down exponential correlation function

Finally, the analyst can specify that sites are partially correlated if
they are adjacent along a stream network (while ignoring flow
direction). This results in an Onstein-Uhlenbeck process along the
stream network (Charsley et al. 2023), or an *exponential tail-down*
model:

$$\mathbf{Q}_{domain} = (\mathbf{I} + \mathbf{D} - \mathbf{P})^{T}{\mathbf{I} + \mathbf{D}}^{- 1}(\mathbf{I} + \mathbf{D} - \mathbf{P})$$
where $\mathbf{D}$ is a sparse diagonal matrix with diagonal elements

$$d_{i,j} = {exp}\left( - 2\theta\left| \mathbf{s}_{i},\mathbf{s}_{j} \right| \right)/\left( 1 - {exp}\left( - 2\theta\left| \mathbf{s}_{i},\mathbf{s}_{j} \right| \right) \right)$$
where $\left| \mathbf{s}_{i},\mathbf{s}_{j} \right|$ is the distance
from downstream (parent) node $j$ to upstream (child) node $i$, and
$\theta$ is the O-U parameter governing the decorrelation rate.
Similarly, $\mathbf{P}$ is a sparse matrix containing values
$\rho_{i,j}$, where:
$$\rho_{i,j} = {exp}\left( - \theta\left| \mathbf{s}_{i},\mathbf{s}_{j} \right| \right)/\left( 1 - {exp}\left( - 2\theta\left| \mathbf{s}_{i},\mathbf{s}_{j} \right| \right) \right)$$
such that $\rho_{i,j}$ the regression slope predicting upstream node $i$
from downstream node $j$. The spatial interpolation matrix $\mathbf{A}$
again has row $\mathbf{A}_{i}$ for each sampling or predictive location,
and $\mathbf{A}_{i}$ is nonzero for the two nodes immediately up and
downstream from a given location, with values such that a given location
is predicted as a weighted average based upon the distance from those
two nodes.

For the SPDE and SAR spatial domains, the term $\tau$ is defined such
that precision $\mathbf{Q}_{domain}$ has unit variance. This is done
because the spatial domain always arises in combination with other
parameters, which are used to define the variance of the associated
spatial variable.

## Structural equation models

tinyVAST also involves specifying a structural equation model (SEM).
This SEM be viewed either:

1.  *Weak interpretation*: as an expressive interface to parameterize
    the correlation among variables, using as many or few parameters as
    might be appropriate; or
2.  *Strong interpretation*: as a structural causal model, allowing
    predictions about the consequence of counterfactual changes to the
    system.

To specify a SEM, the user uses *arrow notation* derived from package
`sem` (Fox 2006), as described in TMB by (J. Thorson and Kristensen
2024). For example, to specify a linear model this involves:

``` r
w1 -> w2, b_12
w1 <-> w1, sd_1
w2 <-> w2, sd_2
```

This then estimates a single slope parameter (represented with a
one-headed arrow), as well as the variance of $W_{1}$ and $W_{2}$
(specified with two-headed arrows). In a more complicated case, $W_{1}$
might cause $W_{2}$, which in turn causes $W_{3}$. This is then
represented as:

``` r
# Path, parameter_name, start value
w1 -> w2, b_12, 0
w2 -> w3, b_23, 0
w1 <-> w1, s_1, 1
w2 <-> w2, s_2, 1
w3 <-> w3, s_3, 1
```

SEM interactions can be as complicated or simple as desired, and can
include:

1.  Latent variables and loops (i.e., they are not restricted to
    directed acyclic graphs);
2.  Values that are fixed a priori, where the `parameter_name` is
    provided as `NA` and the starting value that follows is the fixed
    value;
3.  Values that are mirrored among path coefficients, where the same
    `parameter_name` is provided for multiple rows of the text file.

In this preceding example, path coefficients for one-headed arrows then
define path matrix $\mathbf{P}$:

$$\mathbf{P} = \begin{pmatrix}
0 & 0 & 0 \\
b_{12} & 0 & 0 \\
0 & b_{23} & 0
\end{pmatrix}$$

and coefficents for two-headed arrows define the Cholesky $\mathbf{G}$
of the exnogenous covariance matrix $\mathbf{G}^{T}\mathbf{G}$:

$$\mathbf{G} = \begin{pmatrix}
s_{1} & 0 & 0 \\
0 & s_{2} & 0 \\
0 & 0 & s_{3}
\end{pmatrix}$$

These matrices are define a simultaneous equation model:

\$\$ \mathbf{ w = P w + \epsilon} \\ \mathbf\epsilon \sim \mathrm{MVN}(
\mathbf 0, \mathbf G^T \mathbf G ) \$\$ where the variance
${Var}(\mathbf{w}) = \left( {\mathbf{I} - \mathbf{P}} \right)^{- 1}\mathbf{G}^{2}\left( {\mathbf{I} - \mathbf{P}}^{T} \right)^{- 1}$.
This then results in a sparse precision matrix:

$$\mathbf{Q} = \left( {\mathbf{I} - \mathbf{P}}^{T} \right)\mathbf{G}^{- 1}\mathbf{G}^{- T}\left( {\mathbf{I} - \mathbf{P}} \right)$$

where this precision matrix $\mathbf{Q}$ is then used as a modular
component in the larger tinyVAST model.

## Dynamic structural equation models

Similarly, tinyVAST involves specifying dynamic structural equation
models (DSEM) (J. T. Thorson et al. 2024). To specify a DSEM, the user
uses *arrow-and-lag notation*. For example, to specify a univariate
first-order autoregressive process:

``` r
w1 -> w1, 1, rho, 0.8
```

If there were four time-intervals ($T = 4$) this would then result in
the path matrix:

$$\mathbf{P} = \begin{pmatrix}
0 & 0 & 0 & 0 \\
\rho & 0 & 0 & 0 \\
0 & \rho & 0 & 0 \\
0 & 0 & \rho & 0
\end{pmatrix}$$

and when the DSEM involves multiple times and variables, the sparse
precision is formed by summing across the Kronecker product of time-lag
and interaction matrices. This DSEM defines a GMRF over a nonseparable
interaction of time and variables, represented by a matrix
$\mathbf{Q}_{time\_ term}$ with dimension $CT \times CT$. The user can
specify a separate *arrow-and-lag* notation to define the precision
matrix for the time-variable interaction $\mathbf{Q}_{time\_ term}$ and
for the space-time-variable interaction $\mathbf{Q}_{spacetime\_ term}$.

The precision matrix $\mathbf{Q}_{time\_ term}$ for the time term
$\mathbf{D}$ is used to define the time-varying intercept for each
variable:

$${vec}(\mathbf{D}) \sim {MVN}\left( \mathbf{0},\mathbf{Q}_{time\_ term} \right)$$
Meanwhile, the space-time term $\mathbf{Q}_{spacetime\_ term}$ is
combined with the spatial precision $\mathbf{Q}_{space\_ term}$ as we
explain next.

## Spatial interactions for SEM and DSEM

tinyVAST uses the SEM and DSEM notation to construct the joint precision
for the space-variable interaction $\mathbf{\Omega}$ with dimension
$S \times C$, and the space-time-variable interaction $\mathbf{E}$ with
dimension $S \times C \times T$. To do so, it constructs the separable
precision for each process:

$${vec}(\mathbf{E}) \sim {MVN}\left( \mathbf{0},\mathbf{Q}_{spacetime\_ term} \otimes \mathbf{Q}_{domain} \right)$$

where the precision matrix
$\mathbf{Q}_{spacetime\_ term} \otimes \mathbf{Q}_{domain}$ has
dimension $STC \times STC$ to match the length of ${vec}(\mathbf{E})$,
and
$${vec}(\mathbf{\Omega}) \sim {MVN}\left( \mathbf{0},\mathbf{Q}_{space\_ term} \otimes \mathbf{Q}_{domain} \right)$$

where the precision matrix
$\mathbf{Q}_{space\_ term} \otimes \mathbf{Q}_{domain}$ has dimension
$SC \times SC$, to match the length of ${vec}(\mathbf{\Omega})$. This
specification generalizes spatial factor analysis (Thorson et al. 2015)
and spatial dynamic factor analysis (J. T. Thorson et al. 2016).

## Generalized additive model

The analyst can specify a generalized additive model using syntax from
package `mgcv` (**wood_generalized_2006?**), and parsed into TMB using
an interface developed by `sdmTMB` (Anderson et al. n.d.). For example
this might involve:

``` r
count ~ year + offset(log_area) + s(depth) + s(species, bs="re")
```

If `year` and `species` are factors and `depth` and `log_area` are
continuous, then this would specify a fixed effect for each level
`year`, a spline smoother for `depth`, using `log_area` as an offset,
and estimating a random intercept for each level of `species`. This
formula is parsed internally to assemble fixed effects in design matrix
$\mathbf{X}$ and the basis functions for spline smoothers and random
effects in design matrix $\mathbf{Z}$. The coefficients $\mathbf{γ}$
associated with smoothers and random effects are then specified to
follow a GMRF:

$${\mathbf{γ}} \sim {GMRF}\left( \mathbf{0},\mathbf{Q}_{gam} \right)$$
where $\mathbf{Q}_{gam}$ is a blockwise diagonal matrix, assembled from
estimated variance parameters and matrices constructed by mgcv.

## Spatially varying coefficients

Finally, the analyst can account for spatial variation in the
relationship between a specified covariate and the estimated response.
This is done by specifying a spatially varying coefficient (SVC) (J. T.
Thorson et al. 2023), and this can be used e.g. to account for spatially
varying differences in catchability/detectability (Grüss et al. 2023) or
to account for nonlocal impacts of regional oceanographic indices (J. T.
Thorson 2019). This involves estimating each SVC $l \in \{ 1,2,...,L\}$
as a zero-centered Gaussian Markov random field while estimating its
corresponding variance $\sigma_{\xi,l}^{2}$:

$$\Xi_{1,l} \sim {GMRF}\left( \mathbf{0},\sigma_{\xi,l}^{- 2}\mathbf{Q}_{domain} \right)$$

Any covariate specified as an SVC it is also typically specified in the
`formula` to estimate a non-zero mean for the SVC. If the estimated
variance of the SVC approaches zero, it then suggests that the covariate
slope does not vary spatially.

## Notation summary

| Symbol              | Description                                             |
|:--------------------|:--------------------------------------------------------|
| $i$                 | Index for each sample, $i$ in $(1,2,...,I)$             |
| $s\lbrack i\rbrack$ | spatial coordinate for sample $i$, $s$ in $(1,2,...,S)$ |
| $t\lbrack i\rbrack$ | time-interval for sample $i$, $t$ in $(1,2,...,T)$      |
| $c\lbrack i\rbrack$ | category for sample $i$, $c$ in $(1,2,...,C)$           |
| $e\lbrack i\rbrack$ | error distribution and link function for sample $i$     |

Table 1: Subscript notation

| Symbol  | Code   | Description             |
|:--------|:-------|:------------------------|
| $y$     | `y_i`  | Observed response data  |
| $p_{1}$ | `p_i`  | first linear predictor  |
| $p_{2}$ | `p2_i` | second linear predictor |

Table 2: Symbol notation, code representation (in model output or in
model template code), and descriptions.

## Works cited

Anderson, Sean C., Eric J. Ward, Philina A. English, Lewis A. K.
Barnett, and James T. Thorson. n.d. “sdmTMB: An R Package for Fast,
Flexible, and User-Friendly Generalized Linear Mixed Effects Models with
Spatial and Spatiotemporal Random Fields.” *Journal of Open Source
Software*, 2022.03.24.485545. Accessed March 10, 2025.
<https://doi.org/10.1101/2022.03.24.485545>.

Charsley, Anthony R., Arnaud Grüss, James T. Thorson, Merrill B. Rudd,
Shannan K. Crow, Bruno David, Erica K. Williams, and Simon D. Hoyle.
2023. “Catchment-Scale Stream Network Spatio-Temporal Models, Applied to
the Freshwater Stages of a Diadromous Fish Species, Longfin Eel
(Anguilla Dieffenbachii).” *Fisheries Research* 259 (March): 106583.
<https://doi.org/10.1016/j.fishres.2022.106583>.

Fox, John. 2006. “Structural Equation Modeling with the Sem Package in
R.” *Structural Equation Modeling-a Multidisciplinary Journal* 13:
465–86.

Grüss, Arnaud, and James T. Thorson. 2019. “Developing Spatio-Temporal
Models Using Multiple Data Types for Evaluating Population Trends and
Habitat Usage.” *ICES Journal of Marine Science* 76 (6): 1748–61.
<https://doi.org/10.1093/icesjms/fsz075>.

Grüss, Arnaud, James T. Thorson, Owen F. Anderson, Richard L.
O’Driscoll, Madison Heller-Shipley, and Scott Goodman. 2023. “Spatially
Varying Catchability for Integrating Research Survey Data with Other
Data Sources: Case Studies Involving Observer Samples,
Industry-Cooperative Surveys, and Predators as Samplers.” *Canadian
Journal of Fisheries and Aquatic Sciences* 80 (10): 1595–615.
<https://doi.org/10.1139/cjfas-2023-0051>.

Lindgren. 2012. “Continuous Domain Spatial Models in R-INLA.” *The ISBA
Bulletin* 19 (4): 14–20.

Lindgren, Finn. 2023. “Fmesher: Triangle Meshes and Related Geometry
Tools.” <https://CRAN.R-project.org/package=fmesher>.

Lindgren, Finn, Håvard Rue, and Johan Lindström. 2011. “An Explicit Link
Between Gaussian Fields and Gaussian Markov Random Fields: The
Stochastic Partial Differential Equation Approach.” *Journal of the
Royal Statistical Society: Series B (Statistical Methodology)* 73 (4):
423–98. <https://doi.org/10.1111/j.1467-9868.2011.00777.x>.

Thorson, James T. 2018. “Three Problems with the Conventional
Delta-Model for Biomass Sampling Data, and a Computationally Efficient
Alternative.” *Canadian Journal of Fisheries and Aquatic Sciences* 75
(9): 1369–82. <https://doi.org/10.1139/cjfas-2017-0266>.

———. 2019. “Measuring the Impact of Oceanographic Indices on Species
Distribution Shifts: The Spatially Varying Effect of Cold-Pool Extent in
the Eastern Bering Sea.” *Limnology and Oceanography* 64 (6): 2632–45.
<https://doi.org/10.1002/lno.11238>.

Thorson, James T., Sean C. Anderson, Pamela Goddard, and Christopher N.
Rooper. 2025. “tinyVAST: R Package With an Expressive Interface to
Specify Lagged and Simultaneous Effects in Multivariate Spatio-Temporal
Models.” *Global Ecology and Biogeography* 34 (4): e70035.
<https://doi.org/10.1111/geb.70035>.

Thorson, James T., Alexander G. Andrews III, Timothy E. Essington, and
Scott I. Large. 2024. “Dynamic Structural Equation Models Synthesize
Ecosystem Dynamics Constrained by Ecological Mechanisms.” *Methods in
Ecology and Evolution* 15 (4): 744–55.
<https://doi.org/10.1111/2041-210X.14289>.

Thorson, James T., Cheryl L. Barnes, Sarah T. Friedman, Janelle L.
Morano, and Margaret C. Siple. 2023. “Spatially Varying Coefficients Can
Improve Parsimony and Descriptive Power for Species Distribution
Models.” *Ecography* 2023 (5): e06510.
<https://doi.org/10.1111/ecog.06510>.

Thorson, James T., James N. Ianelli, Elise A. Larsen, Leslie Ries, Mark
D. Scheuerell, Cody Szuwalski, and Elise F. Zipkin. 2016. “Joint Dynamic
Species Distribution Models: A Tool for Community Ordination and
Spatio-Temporal Monitoring.” *Global Ecology and Biogeography* 25 (9):
1144–58. <https://doi.org/10.1111/geb.12464>.

Thorson, James T., Hans J. Skaug, Kasper Kristensen, Andrew O. Shelton,
Eric J. Ward, John H. Harms, and James A. Benante. 2014. “The Importance
of Spatial Models for Estimating the Strength of Density Dependence.”
*Ecology* 96 (5): 1202–12. <https://doi.org/10.1890/14-0739.1>.

Thorson, James, and Kasper Kristensen. 2024. *Spatio-Temporal Models for
Ecologists*. 1st edition. Boca Raton, FL: Chapman; Hall/CRC.

Thorson, Mark D. Scheuerell, Andrew O. Shelton, Kevin E. See, Hans J.
Skaug, and Kasper Kristensen. 2015. “Spatial Factor Analysis: A New Tool
for Estimating Joint Species Distributions and Correlations in Species
Range.” *Methods in Ecology and Evolution* 6 (6): 627–37.
<https://doi.org/10.1111/2041-210X.12359>.

Ver Hoef, Jay M., Ephraim M. Hanks, and Mevin B. Hooten. 2018. “On the
Relationship Between Conditional (CAR) and Simultaneous (SAR)
Autoregressive Models.” *Spatial Statistics* 25 (June): 68–85.
<https://doi.org/10.1016/j.spasta.2018.04.006>.
