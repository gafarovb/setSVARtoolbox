# setSVARtoolbox

MATLAB toolbox for frequentist inference in set-identified Structural Vector Autoregressions (SVARs), with an analytical implementation of the methods in:

> Bulat Gafarov, Matthias Meier, and José Luis Montiel Olea (2018), **“Delta-Method Inference for a Class of Set-Identified SVARs,”** *Journal of Econometrics* 203(2), 316–327.  
> DOI: https://doi.org/10.1016/j.jeconom.2017.12.004

The published paper is the main methodological reference for the analytical part of this repository. The toolbox also contains a legacy moment-inequality/Bonferroni implementation used for comparison, Monte Carlo utilities, output classes, and older scripts. These auxiliary components should not be confused with the core delta-method procedure in the published paper.

## 1. What the toolbox does

The core problem is inference on impulse responses in an SVAR when equality and/or inequality restrictions identify a **set** of admissible responses rather than a single structural model.

For a reduced-form VAR

\[
Y_t=A_1Y_{t-1}+\cdots+A_pY_{t-p}+B\varepsilon_t,
\qquad E(\varepsilon_t\varepsilon_t')=I,
\]

let

\[
\Sigma=BB',
\qquad
\lambda_{k,i,j}=e_i'C_k(A)Be_j
\]

be the horizon-\(k\) response of variable \(i\) to structural shock \(j\). The paper considers equality and inequality restrictions on a **single column** \(B_j=Be_j\):

\[
Z(\mu)'B_j=0,
\qquad
S(\mu)'B_j\ge 0,
\]

where the reduced-form parameter is \(\mu=(\operatorname{vec}(A)',\operatorname{vec}(\Sigma)')'\).

The toolbox computes:

1. OLS estimates of the reduced-form VAR and its covariance matrix.
2. Reduced-form VMA coefficients \(C_k(A)\) and their derivatives.
3. Sign and zero restrictions, including restrictions at nonzero horizons and cumulative restrictions.
4. The lower and upper endpoints of the identified set for each scalar IRF coefficient.
5. Worst-case delta-method standard errors across candidate active sets.
6. One-sided and two-sided analytical confidence intervals.
7. A legacy moment-inequality/Bonferroni confidence-set procedure for comparison.
8. Asymptotic-normal Monte Carlo experiments for coverage calculations.

## 2. Scope of the published method

The analytical method is for models in which equality and inequality restrictions are imposed on **one structural shock**. This restriction is important: the closed-form active-set calculation is not a general multiple-shock sign-restriction solver.

The final published paper organizes the main methodological results as:

- **Section 4.1 / Theorem 1:** finite active-set algorithm for the maximum and minimum response.
- **Section 4.2 / Theorem 2:** directional differentiability of the identified-set endpoints.
- **Section 4.3 / Theorem 3:** delta-method interval, pointwise frequentist validity, and asymptotic robust-Bayesian credibility under the paper’s additional conditions.
- **Section 4.3.3:** Monte Carlo evidence based on draws of the reduced-form estimator from its asymptotic normal distribution.
- **Section 5:** unconventional monetary policy application.

The repository was developed around earlier drafts, so some comments refer to “GMM 2014” or “2017.” The published reference is the 2018 *Journal of Econometrics* article above.

## 3. Paper-to-code map

| Paper object/result | Code | Role |
|---|---|---|
| Reduced-form OLS estimates \(\widehat A,\widehat\Sigma\) | `estimatedVecAR.m`, `LSestimatesVAR.m` | Estimates the VAR, residual covariance, and asymptotic covariance of reduced-form parameters. |
| VMA coefficients \(C_k(A)\) | `VecAR.getVMAfromAL` | Recursively constructs the reduced-form moving-average representation. |
| Derivatives \(\partial C_k(A)/\partial A\) | `VecAR.getVMAderivatives` | Implements the derivative calculations used in the delta method. |
| Equality and sign restrictions \(Z(\mu),S(\mu)\) | `IDassumptions.m` | Translates user restriction rows into linear constraints on the structural impact vector. |
| Derivatives of restrictions | `IDassumptions.getLinearConstraintsAndDerivatives` | Supplies the derivative terms needed in Theorem 2 / Lemma 2 of the paper. |
| Theorem 1: enumerate possible active sets | `optimizationProblems.initializeSubproblems` | Enumerates combinations of binding sign restrictions, subject to the degrees-of-freedom bound. |
| Closed-form KKT solution for a fixed active set | `subproblemsGivenActiveSet.getProjectedSigma`, `computeKKTpointsAndValues` | Implements the projected covariance matrix and the \(\pm\) closed-form candidate extrema. |
| Feasibility check for inactive inequalities | `subproblemsGivenActiveSet.computeFeasibilityPenalty` | Rejects candidate KKT points violating non-active sign restrictions. |
| Global endpoint over all active sets | `optimizationProblems.getMaxBounds`, `getMinBounds` | Takes the maximum/minimum across feasible active-set candidates. |
| Theorem 2 / derivative of candidate value function | `subproblemsGivenActiveSet.asymptoticStandardDeviationForActiveSet` | Forms the derivative with respect to VAR coefficients and \(\Sigma\), then applies the reduced-form covariance matrix. |
| Equation (4.3): worst-case standard error | `optimizationProblems.getWorstCaseStdMat` | Maximizes the standard error over candidate active sets rather than estimating the unique maximizing active set. |
| Theorem 3: delta-method interval | `SVARanalyticFramework.onesidedUpperIRFCS`, `onesidedLowerIRFCS`, `twoSidedIRFCS` | Adds/subtracts normal critical values times the worst-case standard error from estimated endpoints. |
| Identified-set estimate | `SVARanalyticFramework.identifiedSet` | Returns estimated upper and lower endpoint IRFs. |
| Section 4.3.3 Monte Carlo | `simulatedVecAR.m`, `mcExperiment.m` | Draws reduced-form parameters from the asymptotic Gaussian approximation and evaluates coverage. |
| MSG/GMS-style comparison procedure | `SVARMomentInequalitiesFramework.m`, `stochasticInequalities.m` | Legacy Bonferroni/moment-inequality comparison; not the analytical method proposed in the paper. |
| Paper’s UMP application | `Example2.m`, `integreationTest1.m` | Legacy scripts intended to reproduce/compare the unconventional monetary policy application. See the caveats below. |

## 4. How the analytical algorithm works

### 4.1 Reduce the structural problem to a vector problem

For a particular shock \(j\), the relevant object is \(x=B_j\). Because \(BB'=\Sigma\), an admissible column satisfies

\[
x'\Sigma^{-1}x=1.
\]

For a fixed IRF coefficient, the objective is linear in \(x\):

\[
c'x,
\qquad c=C_k(A)'e_i,
\]

subject to the unit-ellipsoid condition and the equality/sign restrictions.

### 4.2 Enumerate candidate active sets — Theorem 1

The code treats zero restrictions as always active and enumerates subsets of sign restrictions that may bind. This is done in:

- `optimizationProblems.initializeSubproblems`
- `optimizationProblems.getAllsubsetsSize_k_from_N`

For each active set, `subproblemsGivenActiveSet` constructs a matrix `linearEqualityConstraints` that stacks:

- all zero restrictions, and
- the sign restrictions hypothesized to bind.

The number of active sign restrictions is limited by

```matlab
degreesOfFreedom = n - nEqualities - 1;
```

which reflects the dimension lost to equality restrictions plus the ellipsoid normalization.

### 4.3 Closed-form KKT candidate — Lemma 1 / Theorem 1

For an active restriction matrix \(r\), the paper defines a projection matrix of the form

\[
M_{\Sigma^{1/2}r}
=I-\Sigma^{1/2}r(r'\Sigma r)^{-1}r'\Sigma^{1/2}.
\]

The code implements the corresponding effective covariance matrix in

```matlab
subproblemsGivenActiveSet.getProjectedSigma
```

and then evaluates the closed-form candidate value

\[
v(\mu;r)
=\left(c'\Sigma^{1/2}M_{\Sigma^{1/2}r}\Sigma^{1/2}c\right)^{1/2}
\]

in `computeKKTpointsAndValues`.

Both signs, \(+v\) and \(-v\), are considered. `computeFeasibilityPenalty` checks the inactive sign restrictions. Infeasible candidates receive the finite penalty configured by `SVARconfiguration.largeNumber`. Finally:

- `getMaxBounds` selects the largest feasible candidate across active sets;
- `getMinBounds` selects the smallest feasible candidate.

This is the computational implementation of the paper’s finite global active-set algorithm. Unlike random rotation sampling or a generic nonlinear solver, the core analytical routine does not depend on initial conditions and does not approximate the identified set by a finite random grid.

### 4.4 Differentiation — Theorem 2

The derivatives required by the delta method enter from two sources:

1. the objective \(C_k(A)\), and
2. the identifying restrictions when they depend on reduced-form parameters.

`VecAR.getVMAderivatives` computes derivatives of the VMA coefficients with respect to the autoregressive coefficients. `IDassumptions.getLinearConstraintsAndDerivatives` constructs derivatives of the equality/sign restrictions.

For each candidate active set, `subproblemsGivenActiveSet.asymptoticStandardDeviationForActiveSet` combines:

- the derivative of the objective,
- the derivative of the active restrictions weighted by the implied Lagrange multipliers,
- the derivative with respect to \(\Sigma\), and
- the estimated covariance matrix of reduced-form parameters.

The code parameterizes the covariance matrix using `vech(Sigma)` rather than the paper’s `vec(Sigma)` notation. `converter.getVechFromVec` and `converter.getVecFromVech` provide the necessary duplication/elimination mappings.

### 4.5 Worst-case standard error — equation (4.3)

A key feature of the paper is that inference does not require consistently selecting the active set that generates the endpoint. The proposed standard error protects against multiple candidate active sets by taking a maximum.

That logic is implemented by

```matlab
optimizationProblems.getWorstCaseStdMat
```

which takes the largest candidate standard deviation over all enumerated active sets.

### 4.6 Delta-method interval — Theorem 3

`SVARanalyticFramework` is the high-level implementation of the paper’s analytical inference procedure.

- `onesidedUpperIRFHat`: estimated upper endpoint.
- `onesidedLowerIRFHat`: estimated lower endpoint.
- `identifiedSet`: both estimated endpoints.
- `asymptoticStdDeviations`: worst-case delta-method standard errors.
- `onesidedUpperIRFCS`: upper one-sided confidence bound.
- `onesidedLowerIRFCS`: lower one-sided confidence bound.
- `twoSidedIRFCS`: two-sided confidence set using Bonferroni splitting of the two endpoint errors.

The paper writes the interval as

\[
\left[
\widehat v_L-z_{1-\alpha/2}\widehat\sigma/\sqrt T,
\widehat v_U+z_{1-\alpha/2}\widehat\sigma/\sqrt T
\right].
\]

In the code, `estimatedVecAR.getCovarianceOfThetaT` already divides the reduced-form asymptotic covariance by \(T\), so `asymptoticStandardDeviationForActiveSet` returns the finite-sample-scale standard deviation. The confidence-set methods therefore multiply directly by the normal critical value rather than dividing by `sqrt(T)` again.

Theorem 3’s robust-Bayesian statement is a **large-sample property of this delta-method interval** under the paper’s assumptions. It is not evidence that the current repository contains a separate full robust-Bayesian posterior implementation.

## 5. Main object graph

```text
SVAR                                  top-level facade
├── VecARmodel
│   ├── estimatedVecAR                estimated reduced-form VAR
│   │   └── LSestimatesVAR            OLS, Sigma, Omega, VMA, derivatives
│   └── simulatedVecAR                asymptotic-normal reduced-form draw
├── IDassumptions                     zero/sign restrictions and derivatives
├── analytic : SVARanalyticFramework  GMM analytical inference
│   └── optimizationProblems          all candidate active sets
│       └── subproblemsGivenActiveSet closed-form KKT candidates + gradients
└── MIframework : SVARMomentInequalitiesFramework
    └── stochasticInequalities        legacy Bonferroni/GMS comparison
```

`IRF` and `IRFcollection` are output/plotting classes used by both inference paths.

## 6. Basic usage

### 6.1 Default repository example

With the repository root on the MATLAB path, the no-argument constructor uses the configuration in `SVARconfiguration.m`, the data in `MSG/data.csv`, and restrictions in `MSG/restMat.dat`:

```matlab
model = SVAR;

% Inspect identifying restrictions
model.tableForm()

% Estimated identified set: upper and lower endpoint IRFs
idSet = model.estimatedIdentifiedSet();

% 68% analytical delta-method confidence set
cs68 = model.IRFtwoSidedCS(0.68, 'Analytic');

% Plot endpoints and confidence bounds
panel = join(idSet, cs68);
plotPanel(panel);
```

The default `MSG` configuration is a legacy comparison design. It is not the four-variable UMP application in the 2018 paper.

### 6.2 Creating a custom model

```matlab
config = SVARconfiguration;
config.nLags = 12;
config.nNoncontemoraneousHorizons = 36;
config.scedasticity = 'homo';
config.isCumulativeIRF = 'yes';
config.SVARlabel = 'MySVAR';

% data is T x n
% names and units are 1 x n cell arrays
TSdescription = [names; units];
dataset = multivariateTimeSeries(data, TSdescription);
rf = estimatedVecAR(config, dataset);

% Restriction columns:
% [variable_index, horizon, sign_or_zero, cumulative]
% sign_or_zero = +1, 0, or -1
restMat = [
    1  0   1  0;
    2  0   1  0;
    3  0  -1  0;
    4  0   0  0
];

ID = IDassumptions(restMat, 'shock label');
model = SVAR(rf, ID);

idSet = model.estimatedIdentifiedSet();
cs = model.IRFtwoSidedCS(0.68, 'Analytic');
```

The restriction pattern above matches Table 1 of the published UMP application if the variable ordering is CPI, industrial production, 2-year Treasury rate, and federal funds rate.

## 7. Restriction syntax

`IDassumptions` accepts rows documented in the source as

```text
Var   Hor   S/Z   Cum   Shk
```

The implemented fields are:

- `Var`: 1-based index of the variable whose response is restricted.
- `Hor`: horizon, with impact equal to `0`.
- `S/Z`: `+1` for nonnegative, `-1` for nonpositive, `0` for equality to zero.
- `Cum`: intended indicator for a cumulative restriction (`1`) versus a point response (`0`).
- `Shk`: accepted in five-column input but not used by the present `IDassumptions` calculations; the analytical theory and current implementation treat the identifying scheme as restrictions on a single shock.

Four-column input is automatically extended with a fifth column of ones by `IDassumptions.longForm`.

### Important implementation caveat for cumulative restrictions

The current version of

```matlab
IDassumptions.aRestrictionIsCumulative
```

returns `assumptionsMatrixInput(aRestriction,2) + 1`, i.e. it reads the **horizon column**, not the cumulative-indicator column. This appears to be a code defect. Impact restrictions are unaffected because cumulative and non-cumulative responses coincide at horizon zero, but restrictions at positive horizons can be constructed incorrectly. See “Current repository audit notes” below.

## 8. Reduced-form estimation

### `multivariateTimeSeries`

Container for a `T x n` matrix and series metadata.

Public methods:

- `multivariateTimeSeries(tsInColumns, TSdescription)` — constructor.
- `countTS` — number of variables.
- `countTimePeriods` — sample size.
- `getYX(nLags)` — creates dependent-variable and regressor matrices with intercept and lags.
- `getNames`, `getUnitsOfMeasurement` — metadata accessors.

### `estimatedVecAR`

Concrete estimated reduced-form VAR, derived from abstract `VecAR`.

Public methods:

- constructor `estimatedVecAR(config,dataset)`;
- `readDataFromFile`;
- `getNames`, `getUnitsOfMeasurement`, `getN`, `getT`, `countParameters`;
- `getSigma`, `getTheta`, `getCovarianceOfThetaT`;
- `getVMA_ts_sh_ho`, `getVMADerivatives_ts_sh_ho_dAL`;
- `stationarityTest`;
- `optimalNlagsByIC` — AIC/BIC/HQ lag selection.

Static helper:

- `readCSVheader`.

### `LSestimatesVAR`

Performs the actual reduced-form estimation.

Public methods:

- constructor `LSestimatesVAR(...)`;
- `getT`, `getThetaHat`, `getAL`, `getOmega`;
- `getBicAicHQic`;
- `getVMA_ts_sh_ho`, `getSigma`, `getVMADerivatives`.

Private/internal methods:

- `computeLSestimates` — OLS and residual covariance;
- `computeCovariance` — reduced-form covariance dispatch;
- `computeVMAandDerivatives`;
- `killiansBootstrap` — legacy bootstrap-after-bootstrap code based on Kilian (1998), not integrated into the current analytical inference path;
- `computeCovarianceOfTheta` — homoskedastic or heteroskedastic reduced-form covariance estimator.

## 9. `VecAR`: common reduced-form mathematics

`VecAR` is the abstract superclass shared by estimated and simulated reduced-form models.

Instance methods:

- `getIRFObjectiveFunctions` — point or cumulative VMA objectives according to configuration;
- `getIRFObjectiveFunctionsDerivatives` — corresponding derivatives;
- `getConfig`;
- `precomputeCache`;
- `ALSigmaFromThetaNandP` — maps the parameter vector back to \(A(L)\) and \(\Sigma\);
- `getMaxHorizons`.

Static methods:

- `thetaFromALSigma` — stacks `vec(A)` and `vech(Sigma)`;
- `getVMAfromAL` — VMA recursion;
- `getVMAderivatives` — analytical VMA derivatives;
- `getNPfromAL`;
- `stationarity` — companion-matrix stationarity diagnostic.

Abstract interface implemented by subclasses:

- `getVMA_ts_sh_ho`;
- `getVMADerivatives_ts_sh_ho_dAL`;
- `getSigma`;
- `getN`;
- `getTheta`;
- `getNames`;
- `getUnitsOfMeasurement`;
- `getCovarianceOfThetaT`.

## 10. Identification object: `IDassumptions`

`IDassumptions` stores the single-shock restriction scheme and converts it into the matrices required by the active-set algorithm.

Methods:

- constructor `IDassumptions(restMatShort,shockLabel)`;
- `getMaxTS`, `getMaxHorizon`;
- `tableForm`, `disp`;
- `getRestMat`;
- `countSignRestrictions`, `countZeroRestrictions`, `countActiveSets`;
- `getHorizonOfaRestriction`;
- `getTsOfaRestriction`;
- `aRestrictionIsCumulative`;
- `aRestrictionIsEquality`;
- `convertArestrictionToGeq`;
- `assertAssumptions`;
- `getLinearConstraintsAndDerivatives`;
- static `longForm`.

`getLinearConstraintsAndDerivatives` is especially important for the paper: it constructs the sample analogues of \(S(\mu)\), \(Z(\mu)\), and their derivatives with respect to the VAR coefficients.

## 11. Top-level facade: `SVAR`

`SVAR` joins the reduced form and identification scheme and loads two inference frameworks.

Construction and access methods:

- constructor `SVAR(VecARmodel,ID)`;
- `generateSamplesFromAsymptoticDistribution`;
- `loadMIframework`;
- `loadAnalyticFramwork`;
- `getConfig`, `getN`, `getT`, `getSigma`, `getTheta`;
- `getNamesOfTS`, `getUnitsOfMeasurement`, `getTSDescription`;
- `tableForm`, `disp`;
- `getShockLabel`, `getMaxHorizons`;
- `getCovarianceOfThetaT`, `getRestMat`;
- `getLinearConstraintsAndDerivatives`;
- `getIRFObjectiveFunctions`.

Inference/output methods:

- `enforceIRFRestrictions` — clips reported IRFs to explicitly imposed zero/sign restrictions;
- `IRFtwoSidedCS(level,type)` — facade for `Analytic` or `MI_Bonferroni` confidence sets;
- `estimatedIdentifiedSet` — analytical estimated identified set.

## 12. Analytical inference classes

### `SVARanalyticFramework`

Methods:

- constructor;
- `setOptimizationProblems` — lazy creation of active-set problems;
- `asymptoticStdDeviations`;
- `onesidedUpperIRFHat`;
- `onesidedLowerIRFHat`;
- `identifiedSet`;
- `onesidedUpperIRFCS`;
- `onesidedLowerIRFCS`;
- `twoSidedIRFCS`;
- `descriptionPointEstimates`, `descriptionStd`.

This is the preferred entry point for the Gafarov–Meier–Montiel Olea analytical method.

### `optimizationProblems`

Methods:

- constructor;
- `setSQRTSigma`, `getSQRTSigma`;
- `initializeSubproblems`;
- `getSigma`, `getCovarianceOfThetaT`;
- `getLinearConstraintsAndDerivatives`, `getObjectiveFunctions`, `getConfig`;
- `countSubproblems`, `getHorizons`, `getN`;
- `getMaxBounds`, `getMinBounds`;
- `countInequalityRestrictions`, `countEqualityRestrictions`;
- `getWorstCaseStdMat`;
- static `getAllsubsetsSize_k_from_N`.

### `subproblemsGivenActiveSet`

Methods:

- constructor;
- `setActiveConstraintsAndDerivatives`;
- `getProjectedSigma`;
- `computeKKTpointsAndValues`;
- `computeFeasibilityPenalty`;
- `computePenalizedMaximum`;
- `computePenalizedMinimum`;
- `asymptoticStandardDeviationForActiveSet`.

This class contains the closest code analogue of the mathematical derivations behind Theorems 1 and 2.

## 13. Moment-inequality / Bonferroni comparison path

The repository also provides a second inference option:

```matlab
model.IRFtwoSidedCS(level, 'MI_Bonferroni')
```

`SVAR.m` labels this as a Bonferroni-type method associated with Moon, Schorfheide, and Granziera. It is a **comparison method**, not the analytical method proposed by Gafarov, Meier, and Montiel Olea (2018).

### `SVARMomentInequalitiesFramework`

Methods:

- constructor;
- `setMomentInequalities`;
- `computeFirstColumnOfSigmaSqrtAtSphericalGridPoint`;
- `computeInequlaitySlackAtSphericalGridPoint`;
- `computeEqualityResidualAtSphericalGridPoint`;
- `computeObjectiveFunctionsAtSphericalGridPoint`;
- `twosidedIRFCSbonferroni`.

### `stochasticInequalities`

Main methods:

- constructor — generates asymptotic-normal reduced-form samples and random unit-sphere grid points;
- `getNsamples`, `getNgridPoints`;
- `computeBonferroniCS`;
- `descriptionBonferroniCS`;
- `generalizedMomentSelection`;
- `computeQuantileOfModifiedMethodOfMoments`;
- `testGridPoint`, `testAllPoints`;
- `computeConditionalIRFCS`;
- `computeCSObjectiveFunctionsAtGridPoint`;
- `resampleObjectiveFunctions`, `resampleResiduals`, `resampleSlacks`;
- `emptyIRF`;
- static `generateNdimUnitSphereGrid`;
- static `modifiedMethodOfMoments`;
- static `assertIfNoFeasiblePointsFound`.

The generalized moment selection threshold uses `SVARconfiguration.andrewsSoaresTunningSequence`.

## 14. Monte Carlo code

### `simulatedVecAR`

Represents a draw of reduced-form parameters from the estimated asymptotic Gaussian approximation.

Methods:

- constructor;
- `getN`;
- `generateThetaFromNormal`;
- `simulateVAR` — legacy residual-resampling VAR simulator;
- `getSigma`, `getAL_n_x_np`, `getCovarianceOfThetaT`, `getTheta`;
- `getVMA_ts_sh_ho`, `getVMADerivatives_ts_sh_ho_dAL`;
- `getNames`, `getUnitsOfMeasurement`.

### `mcExperiment`

Coverage simulation for the analytical intervals.

Methods:

- constructor;
- `getNumberOfSimulations`;
- `testCoverageOnesidedUpperIRFCSAnalytic`;
- `testCoverageOnesidedLowerIRFCSAnalytic`;
- `testCoverageTwoSidedIRFCSAnalytic`;
- `computeCoverageFrequency`;
- `MCwaitBar`.

The simulation design—drawing reduced-form estimates directly from an asymptotic normal distribution—is closely related to the Monte Carlo experiment in Section 4.3.3 of the published paper.

## 15. Output and utility classes

### `IRF`

Scalar-series IRF container with overloaded arithmetic/comparison operators.

Methods:

- constructor;
- `double`;
- `plus`, `mtimes`, `minus`, `le`, `ge`;
- `setValues`;
- `nNoncontemoraneousHorizons`;
- `setDescription`, `setDescriptionField`;
- `getLabelTS`;
- `plot`;
- `getMarker`.

### `IRFcollection`

Collection of `IRF` objects for multiple series and/or multiple lines.

Methods:

- constructor;
- `setValues`, `setDescription`, `setDescriptionField`;
- `matrixForm`, `arrayForm`, `double`, `disp`;
- `join`;
- `plotPanel`.

### `bootstrapSample`

Small wrapper around resampled arrays:

- constructor;
- `std`;
- `getN`;
- `standardized`;
- `studentized`;
- `quantile`.

### `converter`

Matrix-vectorization helpers:

- `getVecFromVech(n)` — duplication map satisfying `vec(Sigma)=D*vech(Sigma)`;
- `getVechFromVec(n)` — elimination map from `vec` to `vech`.

### `tensorOperations`

Static tensor helpers:

- `convWithVector`;
- `vectorDotTensor`;
- `rowForm`.

### `waitBarCustomized`

Progress-display utility:

- constructor;
- `setMessage`;
- `showProgress`;
- `elapsedTime`;
- `estimatedTime`.

## 16. Configuration

`SVARconfiguration.m` contains the default settings. Important fields include:

- `nLags` — VAR lag order;
- `nLagsMax` — maximum lag order for information-criterion search;
- `nNoncontemoraneousHorizons` — number of non-impact horizons;
- `scedasticity` — `'homo'` or `'hetero'`;
- `isCumulativeIRF` — `'yes'` or `'no'`;
- `largeNumber` — finite infeasibility penalty corresponding to the constant used in the active-set algorithm;
- `smallNumber` — feasibility tolerance;
- `nGridPoints`, `nBootstrapSamples`, `bonferroniStep1` — moment-inequality comparator settings;
- `masterSeed`, `MaxSimulations` — simulation settings.

The default configuration points to:

- `MSG/data.csv`
- `MSG/restMat.dat`

and applies an MSG-specific preprocessing function.

## 17. Empirical and comparison scripts

### `Example1.m`

Very early demonstration script for the default MSG data. It uses method names that no longer exist in the current `SVAR` facade and should be treated as stale legacy code.

### `Example2.m`

Later OOP-style script for a four-variable unconventional-monetary-policy model. Its architecture is useful as an example of constructing a custom `SVAR` object, but it is not currently a self-contained replication script:

- it loads `data/UMP`, which is not present in the current repository tree;
- it uses `level` without defining it locally;
- its fourth impact restriction is coded as negative, while Table 1 of the published paper imposes a **zero** response of the federal funds rate.

### `integreationTest1.m`

Legacy integration comparison against saved output from an older toolbox. Its four impact restrictions are `+,+,-,0`, which correspond to Table 1 of the published UMP application. It also depends on external files (`data/UMP` and `outputToolbox1/test_correct.mat`) that are not present in the current repository tree.

### `comparisonWithMSG.m`

Runs analytical and moment-inequality/Bonferroni intervals side by side for a legacy comparison design. It also depends on `data/UMP` and an externally supplied `level` variable.

### `MSG/`

Contains the default dataset and restriction file:

- `MSG/data.csv`: Output, Inflation, Interest Rate, Real Money.
- `MSG/restMat.dat`: default zero/sign restrictions used by the no-argument constructor.

The MSG dataset is a comparator/example and is not the data used for the UMP application in Section 5 of Gafarov, Meier, and Montiel Olea (2018).

## 18. Current repository audit notes

This repository is research code from the development period of the paper rather than a polished modern replication package. The following issues are visible from a static review of the current `master` branch and should be addressed before relying on it for new empirical work.

### 18.1 Likely bug in cumulative-restriction handling

`IDassumptions.aRestrictionIsCumulative` currently reads the horizon column and adds one instead of reading column 4 (`Cum`). This makes positive-horizon restriction construction inconsistent with the documented restriction matrix. Impact restrictions are not affected because the point and cumulative IRF coincide at horizon zero.

### 18.2 Heteroskedastic covariance path needs verification

In the heteroskedastic branch of `LSestimatesVAR.computeCovarianceOfTheta`, the code calls

```matlab
getVechFromVec(n)
```

without the class qualifier used elsewhere (`converter.getVechFromVec(n)`). The current repository does not define a separate top-level `getVechFromVec.m`, so the `'hetero'` option should be tested/fixed before use.

### 18.3 Legacy examples use stale API or missing inputs

- `Example1.m` calls old facade methods such as `onesidedLowerIRFHatAnalytic` and `onesidedUpperIRFCSAnalytic` that are no longer defined on `SVAR`.
- `Example2.m` and `comparisonWithMSG.m` rely on an undefined `level` variable.
- UMP data and saved comparison output referenced by the scripts are absent from the current repository tree.

### 18.4 `Example2.m` does not exactly match the published Table 1 restrictions

The published UMP identification scheme is:

```text
CPI   >= 0
IP    >= 0
2yTB  <= 0
FF     = 0
```

`integreationTest1.m` uses this pattern, but `Example2.m` currently codes the fourth restriction as `-1` rather than `0`.

### 18.5 Old manual overstates what is currently implemented

`manual/SVAR_toolbox_V2.tex` lists numerical `fmincon`, Bayesian grid-search, frequentist grid-search, robust-Bayesian, and projection methods. Searches of the current tracked MATLAB files do not reveal implementations of `fmincon`, a dedicated Bayesian sampler, or a projection routine. The manual should therefore be viewed as an old development plan rather than an accurate API description.

The current code visibly implements:

- the analytical delta-method procedure;
- a moment-inequality/Bonferroni comparison procedure;
- asymptotic-normal simulation/coverage utilities.

### 18.6 Robust-Bayesian credibility is a theorem, not a separate routine here

The published paper proves asymptotic robust-Bayesian credibility of the delta-method interval under additional assumptions. The current repository does not expose a separate robust-Bayesian posterior procedure in the main code path.

### 18.7 Legacy bootstrap code is not integrated

`LSestimatesVAR.killiansBootstrap` and `simulatedVecAR.simulateVAR` contain older bootstrap/simulation code and references to variables/functions from an earlier architecture. They are not used by the current analytical confidence-interval path.

### 18.8 Stationarity diagnostic contains an explicit FIXME

`VecAR.stationarity` calls `eigs(A,size(A,1)-2)` and contains the source comment `FIXME: WHY -2?`. The diagnostic should be rechecked before being treated as a formal stationarity test.

## 19. MATLAB dependencies

The exact minimum MATLAB release has not been reconstructed. The current code uses functions/features including:

- MATLAB class definitions (`classdef`);
- `table`;
- `lagmatrix`;
- `norminv`, `quantile`, `combnk`/combinatorial utilities;
- `eigs`;
- optional `parfor` / `gcp` for the moment-inequality grid calculation;
- plotting and `waitbar`.

Accordingly, Econometrics/Statistics functionality is needed for several routines, and the Parallel Computing Toolbox is optional for the parallel branch.

## 20. Interpretation of confidence sets

The analytical intervals are scalar-by-scalar IRF confidence sets. The published paper establishes **pointwise** frequentist validity under its assumptions; it does not claim uniform validity over the same broad class covered by projection procedures. The paper explicitly discusses projection methods as more conservative alternatives with different theoretical guarantees.

The estimated identified set and the confidence set should therefore not be conflated:

- `estimatedIdentifiedSet` estimates the range of structural responses compatible with the identifying restrictions at the estimated reduced form;
- `IRFtwoSidedCS(...,'Analytic')` expands those endpoint estimates to account for sampling uncertainty in the reduced-form parameters.

## 21. Suggested citation

If this code is used for the analytical set-identified SVAR procedure, cite:

```bibtex
@article{GafarovMeierMontielOlea2018,
  author  = {Gafarov, Bulat and Meier, Matthias and Montiel Olea, Jos\'e Luis},
  title   = {Delta-Method Inference for a Class of Set-Identified SVARs},
  journal = {Journal of Econometrics},
  year    = {2018},
  volume  = {203},
  number  = {2},
  pages   = {316--327},
  doi     = {10.1016/j.jeconom.2017.12.004}
}
```

## 22. Recommended modernization work

A clean next revision of the repository would ideally:

1. fix and unit-test cumulative restriction parsing;
2. fix and test the heteroskedastic covariance path;
3. replace the stale examples with one minimal runnable example and one exact UMP replication example;
4. add the UMP replication data or explicit download/preparation instructions if redistribution permits;
5. turn `integreationTest1.m` into automated regression tests with numerical tolerances;
6. test Theorem 1 endpoints against brute-force sphere sampling for small systems;
7. test analytical gradients against finite differences;
8. test equation (4.3) standard errors against direct numerical differentiation;
9. separate the GMM analytical method from MSG/GMS comparison code into clear subfolders;
10. archive or remove manual claims for methods not present in the tracked source;
11. add a license file and explicit MATLAB/toolbox version requirements;
12. add continuous integration using MATLAB’s unit-test framework.

---

### Short summary

The central research contribution implemented here is the **finite active-set analytical algorithm plus worst-case delta-method inference** from Gafarov, Meier, and Montiel Olea (2018). The core chain is:

```text
estimatedVecAR / LSestimatesVAR
        ↓
IDassumptions
        ↓
SVAR
        ↓
SVARanalyticFramework
        ↓
optimizationProblems
        ↓
subproblemsGivenActiveSet
        ↓
identified-set endpoints + worst-case delta-method standard errors
        ↓
scalar IRF confidence intervals
```

The moment-inequality code, simulation classes, old manual, and legacy scripts are auxiliary components rather than separate contributions of the 2018 analytical method.