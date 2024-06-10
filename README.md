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

In the discrete case, the variable $x$ takes discrete values, typically positive integers. The probability mass function (PMF) for a discrete power-law distribution is given by $p(x) = \frac{x^{-\alpha}}{\zeta(\alpha, x_{min})}$, where ${\zeta(\alpha, x_{min})}$ is the generalized or Hurwitz zeta function, ensuring the probabilities sum to 1 overall $x \geq x_{min}$.

## Importance of Power-Law Distributions

Power-law distributions are crucial because they often describe the distribution of various types of data, such as:

1. Earthquake magnitudes.
2. Word frequencies in natural language.
3. Wealth distribution.
4. Sizes of cities and towns.
5. Internet traffic patterns.

Understanding power-law distributions allows researchers and analysts to better comprehend the underlying mechanisms that generate such data and to make more accurate predictions and models.

## Project Components

### `plfit` Function

The `plfit` function fits a power-law distribution to empirical data using maximum likelihood estimation. It determines the scaling parameter $\alpha$ and the lower bound $x_{min}$ above which the power-law behavior holds.

#### Key Features:

1. Estimates the scaling parameter $\alpha$.
2. Identifies the lower bound $x_{min}$.
3. Provides diagnostics for assessing the fit quality.

### Mathematical Calculations for `plfit` Function

#### Continuous Case (REAL)

In the continuous case, the power-law distribution is represented as $p(x) = Cx^{-\alpha}$, where $\alpha$ is the scaling parameter and $C$ is the normalization constant. The goal is to estimate $\alpha$ and $x_{min}$ using maximum likelihood estimation (MLE).

1. **Normalization Constant**: The normalization constant $C$ ensures that the total probability integrates to 1 over the range  $x \geq x_{min}$: $C = (\alpha - 1)x_{min}^{(\alpha-1)}$.
2. **Maximum Likelihood Estimation for $\alpha$**: Given a set of observed values $x_{i}$ such that $x_i \geq x_{min}$, the MLE for the scaling parameter $\alpha$ is $\hat{\alpha} = 1 + n \left[ \sum\limits_{i=1}^n \ln \left( \frac{x_i}{x_{\min}} \right) \right]^{-1}$, where $n$ is the number of observations with $x \geq x_{min}$.
3. **Estimating $x_{min}$**: 1. To find the optimal $x_{min}$, the function iteratively tests different values of $x_{min}$​ and selects the one that minimizes the Kolmogorov-Smirnov (KS) statistic, which measures the distance between the empirical distribution function and the fitted power-law model.

#### Discrete Case (INTS)

In the discrete case, the power-law distribution is represented as $p(x) = \frac{x^{-\alpha}}{\zeta(\alpha, x_{min})}$, where ${\zeta(\alpha, x_{min})}$ is the Hurwitz zeta function.

1. **Hurwitz Zeta Function**: The Hurwitz zeta function is defined as $\zeta(\alpha, x_{\min}) = \sum\limits_{n=0}^{\infty} (n + x_{\min})^{-\alpha}$.
2. **Maximum Likelihood Estimation for $\alpha$**: For discrete data, the MLE for $\alpha$ is found by solving the following equation numerically $\frac{\zeta'(\hat{\alpha}, x_{\min})}{\zeta(\hat{\alpha}, x_{\min})} = - \frac{1}{n} \sum\limits_{i=1}^n \ln x_i$, where $\zeta'(\hat{\alpha}, x_{\min})$ is the derivative of the Hurwitz zeta function with respect to $\alpha$.
3. **Estimating $x_{min}$**: Similar to the continuous case, the optimal $x_{min}$ is determined by iteratively testing different values and minimizing the KS statistic.

### `plpva` Function

The `plpva` function performs a statistical test to determine whether the power-law model is a good fit for the given data. It uses a p-value to quantify the goodness-of-fit by comparing the empirical data to synthetic datasets generated from the fitted power-law distribution.

#### Key Features:

1. Computes the p-value for the power-law fit.
2. Generates synthetic datasets for comparison.
3. Assesses the statistical significance of the fit.

### Mathematical Calculations for `plpva` Function

### Steps:
#### 1. Goodness-of-Fit (GoF) Statistic Calculation

The GoF statistic is the maximum absolute difference between the empirical cumulative distribution function (CDF) and the fitted CDF. This is similar to the Kolmogorov-Smirnov (K-S) statistic. The K-S statistic measures the maximum absolute difference between the empirical CDF and the theoretical CDF. This helps assess how well the power-law model fits the data.

The CDF of a random variable $X$ is a function that maps values to their cumulative probabilities. For a given value $x$, the CDF $F(x)$ gives the probability that the random variable $X$ is less than or equal to $x$.

#### Significance of the CDF

The CDF provides a complete description of the distribution of a random variable. It allows us to calculate the probability that the random variable falls within a certain range. It’s useful for comparing different distributions and for goodness-of-fit tests.
#### Continuous Case (REAL)

$\text{Goodness of Fit (GoF)} = \max \left| \frac{i-1}{n_z} - \left( 1 - \left( \frac{x_{\text{min}}}{z_i} \right)^{\alpha - 1} \right) \right|$, where $n_z$ is the number of data points, $z ≥ x_{min}$​, $z_i$ is the sorted data points $z$, and $\alpha = 1 + \frac{n_z}{\displaystyle \sum_{i=1}^{n_z} \ln\left(\frac{z_i}{x_{\min}}\right)}$.

#### Discrete Case (INTS)

For integer-valued data, $α$ is estimated using the method of maximum likelihood, where $L(\alpha) = -\alpha \displaystyle \sum_{i=1}^{n_z} \ln(z_i) - n_z \ln\left( \zeta(\alpha) - \displaystyle \sum_{x=1}^{x_{\min}-1} x^{-\alpha} \right)$.

Here, $α$ is the value that maximizes the log-likelihood $L(α)$. $\zeta(\alpha)$ is the Riemann zeta function, which is defined as $\zeta(\alpha) = \sum\limits_{k=1}^\infty k^{-\alpha}$.

$\text{Goodness of Fit (GoF)} = max$ $| CDF_{empirical}(z_i) - CDF_{fitted}(z_i) |$, where the fitted CDF for integer data is calculated using the discrete power-law distribution.

#### 2. Bootstrapping Procedure

1. **Generate Synthetic Data:** Generate synthetic datasets under the null hypothesis that the data follows a power-law distribution. For each bootstrap sample, generate synthetic data that follows the power-law distribution with parameters ( $\alpha$) and ($x_{min}$) estimated from the empirical data.

2. **Calculate GoF for Synthetic Data:** Calculate the GoF statistic for each synthetic dataset in the same way as for the empirical data.

3. **Compare GoF Statistics:** Compare the GoF statistic of the empirical data to the distribution of GoF statistics from the synthetic datasets.

#### 3. p-value Calculation

The p-value is the proportion of bootstrap GoF statistics that are greater than or equal to the empirical Goodness of Fit (GoF) statistic, such that $p = \frac{1}{B}$, where $B$ is the number of bootstrap samples, $\mathbb{I}(.)$ is the indicator function, which is 1 if the condition is true and 0 otherwise, $GoF$ is the GoF statistic for the $i$-th bootstrap sample, and $GoF_{empirical}$ is the GoF statistic for the empirical data.

### Determination of p-value Factors

1. Estimating the scaling parameter ($α$) and the lower bound ($x_{min}$​) using the `plfit` function.
2. Calculating the GoF statistic for the empirical data.
3. Generating multiple synthetic datasets (typically 1000) that follow the power-law distribution with the estimated parameters.
4. Calculating the GoF statistic for each synthetic dataset.
5. Determining the proportion of synthetic GoF statistics that are greater than or equal to the empirical GoF statistic, which gives the p-value.

### Applying the Research

The `plfit` and `plpva` functions implemented in this project utilize the methodologies described in the power-fit distribution original paper; Maximum Likelihood Estimation (MLE) and Goodness-of-Fit Tests.

1. **Maximum Likelihood Estimation (MLE)**: Used for fitting the power-law model to the data, providing robust parameter estimates.
2. **Goodness-of-Fit Tests**: Based on the KS statistic, these tests help determine the plausibility of the power-law model.

The MATLAB code provided in this project closely follows the algorithms and statistical techniques discussed in the paper, ensuring accurate and reliable power-law fits.

## Benefits of Power-Law Fitting

1. **Identifying Patterns:** By fitting power-law distributions, we can identify patterns in data that are not apparent with other distributions.
2. **Modeling Complex Systems:** Power-law models are instrumental in understanding and modeling complex systems with scale-invariant properties.
3. **Predictive Analysis:** Accurate power-law fits allow for better predictive analysis in fields such as finance, risk assessment, and network analysis.
4. **Informing Policy and Decision Making:** Insights from power-law fits can inform policy decisions in economics, urban planning, and disaster management.

## Usage

To use the `plfit` and `plpva` functions, clone this repository.

   ```bash
   git clone https://github.com/omarmnfy/Power-Law-Fit-Distribution-MATLAB.git
   cd Power-Law-Fit-Distribution-MATLAB
```
