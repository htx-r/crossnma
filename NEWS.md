## crossnma, version 1.3.0 (2024-11-28)

### Major changes

* Preferred citation for R package **crossnma**:
  [Hamza et al., 2024, BMC Med Res Meth](https://doi.org/10.1186/s12874-023-02130-0)

### Bug fixes

* crossnma.code():
  - fix bug for network meta-regression of only aggregated data

* league.crossnma():
  - fix bug for aggregated data only

### Internal changes

* R function pairwise() provided by R package **meta** not **netmeta**


## crossnma, version 1.2.0 (2023-09-18)

### Major changes

* Surface under the cumulative ranking (SUCRA) can be calculated

* Studies with zero or all events excluded from network meta-analysis
  with binary outcome

* More flexible network graphs

* Rename argument names containing *n.thin* to *thin* to match
  argument name in R packages **rjags** and **coda**

### User-visible changes

* crossnma.model():
  - new arguments 'sucra', 'small.values', 'cov1.value', 'cov2.value'
    and 'cov3.value' for SUCRAs
  - argument 'run.nrs.n.thin' renamed to 'run.nrs.thin'

* crossnma():
  - argument 'n.thin' renamed to 'thin'


## crossnma, version 1.1.0 (2023-05-08)

### Major changes

* More effect measures implemented, i.e., risk ratio (sm = "RR"), mean
  difference ("MD") and standardised mean difference ("SMD")

* Behaviour of print.crossnma() and print.summary.crossnma() switched
  (to be in line with other print and print.summary functions in R)

### Bug fixes

* crossnma.model():
  - fix an error in computing 'cov.ref' for network meta-regression
    with only IPDs or ADs (not both)

### User-visible changes

* crossnma.model():
  - new argument 'sm' to specify the effect measure
  - new argument 'se' to specify standard error for aggregated data
    with continuous outcome
  - replace list argument 'prior' with arguments 'prior.tau.trt',
    'prior.tau.reg0', 'prior.tau.regb', 'prior.tau.regw',
    'prior.tau.bias', 'prior.pi.low.rct', 'prior.pi.high.rct',
    'prior.pi.low.nrs', 'prior.pi.high.nrs'
  - replace list argument 'run.nrs' with arguments 'run.nrs.var.infl',
    'run.nrs.mean.shift', 'run.nrs.n.adapt', 'run.nrs.trt.effect',
    'run.nrs.n.iter', 'run.nrs.n.burnin', 'run.nrs.n.thin',
    'run.nrs.n.chains'

* crossnma.model(), crossnma(), summary.crossnma():
  - new argument 'level.ma' to specify the level for credible
    intervals

* crossnma(), crossnma.model(), heatplot.crossnma(),
  league.crossnma(), print.crossnma(), summary.crossnma()
  - new argument 'backtransf'

* heatplot.crossnma(), league.crossnma(), summary.crossnma()
  - argument 'order' replaced by 'seq'
  - argument 'exp' replaced by 'backtransf'

### Internal changes

* crossnma.model():
  - R function settings.meta() from **meta** package can be used to
    specify the default for arguments 'level.ma' and 'backtransf'

* print.crossnma(), print.summary.crossnma(), league.crossnma():
  - R function settings.meta() from **meta** package can be used to
    specify the default for argument 'digits'

* R function heatplot() renamed to heatplot.crossnma() which relies on
  generic function heatplot() from R package **netmeta**

* heatplot.crossnma():
  - R function settings.meta() from **meta** package can be used to
    specify the default for argument 'digits'

* New auxiliary functions from **meta** package:
  - chklevel(), deprecated2(), formatN(), formatCI(), rmSpace(),
    is.wholenumber()

* New auxiliary function crossnma.model2pairwise()


## crossnma, version 1.0.0 (2022-04-11)

* first submission to CRAN
