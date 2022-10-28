# BayesHistogram.jl
Optimal histogram binning based on piecewise constant model.

Paper: _Studies in Astronomical Time Series Analysis. VI. Bayesian Block Representations_ [https://arxiv.org/abs/1207.5578]

- [BayesHistogram.jl](#bayeshistogramjl)
  - [Introduction](#introduction)
  - [Installation](#installation)
  - [Simple usage](#simple-usage)
  - [Showcase](#showcase)
    - [it handles weighted data and errors correctly](#it-handles-weighted-data-and-errors-correctly)
    - [bins are determined automatically & optimally](#bins-are-determined-automatically--optimally)
    - [it routinely outperforms common binning rules](#it-routinely-outperforms-common-binning-rules)

## Introduction
Have you ever hated the default histogram binning rules in your favourite plotting/analysis library?

You can try to solve your problem by relying on BayesHistogram.jl! :)

This package constructs the bin sequence that maximises the probability of observing your data, under the assumption that it can be described by a histogram.

Or in other words, the package implements a complicated algorithm that returns the optimal histogram, respecting some simple constraints that can be customised.

For those who can't stand it any longer, skip directly to the usage example or, for more details, to the file `make_plot.jl` or if you want to know all the internal details `test/run_tests.jl`.

The optimal histogram is determined from four ingredients:
1) the likelihood corresponding to a given bin.
2) the a-priori probability of observing that bin.
3) the maximum resolution at which we want to separate the data (usually set as `Inf`).
4) the minimum possible number of observations we want for each bin (this is usually `1`).

How can we intervene on these?

1) is obviously not modifiable, for a long time now the Likelihood principle (or its natural Bayesian extension) has been part of the toolbox of all those who make use of statistics: it is a sound principle that we can always trust :)

2) on the other hand is modifiable, the package implements a wide choice of possibilities, the prior can be chosen from the following alternatives:
- `NoPrior`: for non-Bayesians, it always requires tuning (3) and (4).
- `BIC` (default): Bayesian information criterion, requires no parameters and is asymptotically consistent.
- `AIC`: Akaike information criterion: minimises prediction error, requires no parameters (in some cases adds too many bins, but this can be solved using (2) and (3)).
- `HQIC`: Hannan-Quinn criterion, has an intermediate behaviour between BIC and AIC, is close to being consistent, tries to minimise the prediction error.
- `Significance(p)`: Scargle criterion, a given bin is added if it has a statistical significance greater than `p`.
- `Geometric(gamma)`: by varying the parameter `gamma` we change the average number of bins we want to observe.
- `Pearson(p)`: this is useful when we want bins containing about `N*p` observations, where `N` is the total number of events.

3) If intervening on this parameter is necessary, simply add the keyword argument `resolution = ?` to the `bayesian_blocks` function, a typical value might be `100`.

4) Similar to (3), this can be configured via `min_counts = ?`.

Thank you for reading the (short?) introductory section.

## Installation
```julia
using Pkg
Pkg.add("BayesHistogram")
```

## Simple usage
```julia
using Plots, BayesHistogram
X = exp.(randn(5000)./3)
bl = bayesian_blocks(X)
# plot using "to_pdf"
support, density = to_pdf(bl)
plot(support, density)
# or using "edges" parameter and an histogramming procedure
histogram(X, bins=bl.edges, normalize = :pdf)
```

## Showcase
### it handles weighted data and errors correctly
![plot2.png](plot2.png "")
### bins are determined automatically & optimally 
![plot3.png](plot3.png "")
### it routinely outperforms common binning rules
![plot.png](plot.png "")

