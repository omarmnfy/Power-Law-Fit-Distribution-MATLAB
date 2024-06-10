# Power-Law Fitting Distribution

## Introduction

This project provides implementations of the `plfit` and `plpva` functions, which are essential tools for fitting power-law distributions to empirical data. Power-law distributions are significant in various fields, including physics, biology, economics, and social sciences, due to their ability to describe a wide range of natural and man-made phenomena. Understanding and accurately fitting these distributions allows researchers to better model complex systems and predict future events based on observed data.

The `plfit` function is designed to fit a power-law distribution to a given data set using maximum likelihood estimation (MLE). This method is known for providing robust parameter estimates, particularly in the presence of large fluctuations in the tail of the distribution, which are characteristic of power-law behaviors. The function determines the scaling parameter ${\alpha}$ and the lower bound $x_{min}$​, above which the power-law behavior holds. This is crucial for accurately characterizing the distribution and ensuring the reliability of the model.

The `plpva` function complements `plfit` by performing a statistical test to assess the goodness-of-fit of the power-law model. Using a p-value derived from the Kolmogorov–Smirnov (KS) statistic and likelihood ratios, `plpva` quantifies how well the power-law distribution matches the empirical data compared to synthetic datasets generated from the fitted model. This allows researchers to determine the plausibility of the power-law hypothesis and rule out alternative distributions.

This project builds upon the methodologies presented in the research paper [Power-Law Distributions in Empirical Data](https://arxiv.org/abs/0706.1062) by Aaron Clauset, Cosma Rohilla Shalizi, and M. E. J. Newman (2009). The paper provides a principled statistical framework for detecting and quantifying power-law behavior in empirical data, combining MLE with goodness-of-fit tests. Our MATLAB implementation closely follows the algorithms and statistical techniques discussed in the paper, ensuring accurate and reliable power-law fits.

The MATLAB code for this project has been further developed based on the work presented in another notable research titled [Molecular motors robustly drive active gels to a critically connected state.](https://www.nature.com/articles/nphys2715) This research used the power-law fit MATLAB code to produce significant figures, demonstrating the practical application of these functions in cutting-edge scientific studies. This project aims to document and expand upon this code, providing detailed explanations and enhanced usability for researchers and analysts.

By documenting and explaining the power-law fitting process in detail, this project serves as a comprehensive resource for those looking to understand and apply power-law models to their data. Whether for academic research, data analysis, or practical applications, the tools provided here offer robust solutions for fitting and validating power-law distributions.
## Power-Law Distribution

A power-law distribution is a type of probability distribution that has the form: $P(x) \sim x^{-\alpha}$ where ${\alpha}$ is a positive constant known as the exponent or scaling parameter. In a power-law distribution, the probability of observing a value decreases polynomially with the size of the value.

### Continuous vs. Discrete Power-Law Distributions

A power-law distribution can be characterized differently depending on whether the data is continuous or discrete.

#### Continuous Power-Law Distribution

In the continuous case, the variable $x$ can take any real value above a certain minimum $x_{min}$​. The probability density function (PDF) for a continuous power-law distribution is given by $p(x) = Cx^{-\alpha}$, where $\alpha$ is the scaling parameter and $C$ is a normalization constant that ensures the total probability integrates to 1 over the range $x \geq x_{min}$.

#### Discrete Power-Law Distribution

In the discrete case, the variable xxx takes discrete values, typically positive integers. The probability mass function (PMF) for a discrete power-law distribution is given by $p(x) = \frac{x^{-\alpha}}{\zeta(\alpha, x_{min})}$, where ${\zeta(\alpha, x_{min})}$ is the generalized or Hurwitz zeta function, ensuring the probabilities sum to 1 overall $x \geq x_{min}$.

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

The `plfit` function fits a power-law distribution to empirical data using maximum likelihood estimation. It determines the scaling parameter $\alpha$ and the lower bound $x_{min}$ above which the power-law behavior holds.

#### Key Features:

- Estimates the scaling parameter $\alpha$.
- Identifies the lower bound $x_{min}$.
- Provides diagnostics for assessing the fit quality.

### `plpva` Function

The `plpva` function performs a statistical test to determine whether the power-law model is a good fit for the given data. It uses a p-value to quantify the goodness-of-fit by comparing the empirical data to synthetic datasets generated from the fitted power-law distribution.

#### Key Features:

- Computes the p-value for the power-law fit.
- Generates synthetic datasets for comparison.
- Assesses the statistical significance of the fit.

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
   git clone https://github.com/omarmnfy/Power-Law-Fit-Distribution-MATLAB.git
   cd Power-Law-Fit-Distribution-MATLAB
```

