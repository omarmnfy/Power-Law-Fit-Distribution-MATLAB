# Power-Law Fitting Distribution

# Power-Law Fitting Project

## Introduction

This project provides implementations of the `plfit` and `plpva` functions, which are essential for fitting power-law distributions to empirical data. Power-law distributions are significant in various fields, including physics, biology, economics, and social sciences, due to their ability to describe a wide range of natural and man-made phenomena.

## What is a Power-Law Distribution?

A power-law distribution is a type of probability distribution that has the form:

\[ P(x) \sim x^{-\alpha} \]

where \( \alpha \) is a positive constant known as the exponent or scaling parameter. In a power-law distribution, the probability of observing a value decreases polynomially with the size of the value.

### Continuous Case

For continuous data, a power-law distribution can be represented as:

\[ p(x) = C x^{-\alpha} \]

where \( C \) is a normalization constant.

### Discrete Case

For discrete data, a power-law distribution can be written as:

\[ p(x) = \frac{x^{-\alpha}}{\zeta(\alpha, x_{min})} \]

where \( \zeta(\alpha, x_{min}) \) is the generalized or Hurwitz zeta function.

## Importance of Power-Law Distributions

Power-law distributions are crucial because they often describe the distribution of various types of data, such as:

- Earthquake magnitudes
- Word frequencies in natural language
- Wealth distribution
- Sizes of cities and towns
- Internet traffic patterns

Understanding power-law distributions allows researchers and analysts to better comprehend the underlying mechanisms that generate such data and to make more accurate predictions and models.

## Project Components

### `plfit` Function

The `plfit` function fits a power-law distribution to empirical data using maximum likelihood estimation. It determines the scaling parameter \( \alpha \) and the lower bound \( x_{min} \) above which the power-law behavior holds.

#### Key Features:
- Estimates the scaling parameter \( \alpha \).
- Identifies the lower bound \( x_{min} \).
- Provides diagnostics for assessing the fit quality.

### `plpva` Function

The `plpva` function performs a statistical test to determine whether the power-law model is a good fit for the given data. It uses a p-value to quantify the goodness-of-fit by comparing the empirical data to synthetic datasets generated from the fitted power-law distribution.

#### Key Features:
- Computes the p-value for the power-law fit.
- Generates synthetic datasets for comparison.
- Assesses the statistical significance of the fit.

## Research Background

This project is based on the research paper "Power-Law Distributions in Empirical Data" by Aaron Clauset, Cosma Rohilla Shalizi, and M. E. J. Newman (2009). The paper presents a principled statistical framework for discerning and quantifying power-law behavior in empirical data. The approach combines maximum-likelihood fitting methods with goodness-of-fit tests based on the Kolmogorov–Smirnov (KS) statistic and likelihood ratios.

### Applying the Research

The `plfit` and `plpva` functions implemented in this project utilize the methodologies described in the paper:

1. **Maximum Likelihood Estimation (MLE)**: Used for fitting the power-law model to the data, providing robust parameter estimates.
2. **Goodness-of-Fit Tests**: Based on the KS statistic, these tests help determine the plausibility of the power-law model.

The MATLAB code provided in this project closely follows the algorithms and statistical techniques discussed in the paper, ensuring accurate and reliable power-law fits.

## Benefits of Power-Law Fitting

- **Identifying Patterns:** By fitting power-law distributions, we can identify patterns in data that are not apparent with other distributions.
- **Modeling Complex Systems:** Power-law models are instrumental in understanding and modeling complex systems with scale-invariant properties.
- **Predictive Analysis:** Accurate power-law fits allow for better predictive analysis in fields such as finance, risk assessment, and network analysis.
- **Informing Policy and Decision Making:** Insights from power-law fits can inform policy decisions in economics, urban planning, and disaster management.

## Usage

To use the `plfit` and `plpva` functions, follow these steps:

1. Clone this repository:

   ```bash
   git clone https://github.com/yourusername/power-law-fitting.git
   cd power-law-fitting
