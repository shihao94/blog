---
layout: post
title: Understanding NONMEM estimation methods
date: Aug-19-2025
categories: pharmacometrics
tags: pharmacometrics
permalink: 
banner: ./pics/ODE.png
math: true
---

A couple of days ago, I encountered [a nice blog](https://www.vrognas.com/docs/tools-of-the-trade/nonmem/nm-est) on practical tips of using different NONMEM estimation methods. 
It reminds me of my last summer's attempt to understand the basic principles and differences of these estimation methods. It puzzled me that we tend to use FOCE by default. 
So when should we use other estimation methods like SAEM, and how? Also, there seems to different SAEM methods on different software. How to understand them?

Here, I try to summarize the basic theory of these estimation methods, mainly for my own learning. But I also welcome other opinions on this topic and I will be happy 
to learn new updates in this field. As I work on this summary, I had a lot of interactions with Chatgpt.

# 0. maximum likelihood estimates for single measured data

<img src="{{ '/pics/linear_regression.png' | relative_url }}" width="300" height="300" />

- one subject only has one measurement; Normally-distributed independent errors:
  $Y_i = N(f(x_i, \theta), \sigma^2)$, i = 1, 2, …n,
  where $y_i$ is the dependent variable, $x_i$ is independent variable, $\theta$ is model parameters, $\sigma^2$ is residual variance. $\epsilon_j=y_i-f(x_i, \theta)$ follows a normal distribution and are statistically independent.
- Different approaches exist to get “best” estimates of $\theta$: ordinary least squares, weight least squares, maximal likelihood...
- The likelihood of the data given the parameters is $L(Y \vert \theta)$, which is the y axis value of distribution. In this case, the value of $\theta$ to maximize $L(Y \vert \theta)$ is knows as the maximum likelihood estimate of $\theta$. 

<img src="{{ '/pics/ML.png' | relative_url }}" width="400" height="200" />

- Based on probability density function of normal distribution, we have: for each observation $L(y_i \vert \theta)=\frac{1}{\sigma \sqrt{2\pi}}e^{\frac{(y_i-f(x_i, \theta))^2}{-2\sigma^2}}$. Then, for all the observations: $L(Y \vert \theta)=\prod_{i=1}^{n}\frac{1}{\sigma \sqrt{2\pi}}e^{\frac{(y_i-f(x_i, \theta))^2}{-2\sigma^2}}$.
- OFV = $-2lnL(Y \vert \theta)=\sum_{i=1}^{n}ln2\pi+\sum_{i=1}^{n}ln\sigma^2+\sum_{i=1}^{n}\frac{(y_i-f(x_i, \theta))^2}{\sigma^2}$, the first term is a constant, so OFV can be defined as: $nln\sigma^2+\sum_{i=1}^{n}\frac{(y_i-f(x_i, \theta))^2}{\sigma^2}$. To find the maximal likelihood is to find $\theta$ that minimizes OFV.

# 1. Mixed effects modeling: to get maximum likelihood estimates for repeated measured data

- repeated measured data: Multiple observations on one subject. Observations on the same subject are correlated due to the subject-specific $\theta_i$. $y_i \sim N (f(t_i, \theta), \sigma^2)$
- Mixed-effects models:
    
    Stage 1: within the subject
    
     $Y_{i} = h_i(\theta_i)+e_i$, i = 1,…,N
    
    $e_i \sim N(0, G_i(h_i(\theta_i), \sigma^2))$
    
    $Y_{i}$ is the i-th subject’s observation, $h_i(\theta_i)$ is IPRED, $\sigma^2$ represent residual model parameters
    
    Stage 2: inter-subject variability
    
    $\theta_i \sim N(\mu, \Sigma)$ or $LN(\mu, \Sigma)$
    
**Our goal is to find the best estimate of $\mu$, $\Sigma$, $\sigma^2$ to fit the data**  
- The likelihood of the data $Y$ given the parameters is $L(Y \vert \theta)$. The value of $\theta$ to maximize $L(Y \vert \theta)$ is known as the maximum likelihood estimate of $\theta$
- $\theta_i$ is latent (missing), but we assume $\theta$ follows a log-normal or normal distribution, so we can consider the conditional likelihood at every possible value of $\theta_i$ and get a weighted average (expectation). 
Overall conditional data likelihood from all observations (assuming conditionally independent, i.e. the conditional likelihood for an individual observation does not depend on previous observations) is the multiplication of likelihood of all observations. $L(\mu, \Sigma, \sigma^2)=\prod_{i=1}^{N}\int l_i(Y_i \vert \theta, \sigma^2)p(\theta \vert \mu,\Sigma)d\theta$,
where $l_i(Y_i \vert \theta, \sigma^2)$ defines the fit to the data, and $p(\theta \vert \mu,\Sigma)$ is the prior on $\theta$.  
- The integration is computationally difficult to solve, and various methods try to solve this problem
    1. Directly maximize
    2. Approximate likelihood (FO, FOCE, Laplace)
    3. EM algorithm: iterative solution to 2 simpler problems, plus sampling-based methods. This produces the Exact Maximum Likelihood Estimate
    
# 2. Approximate likelihood (FO, FOCE, Laplace)

## an example to understand FO method

One-compartment PK model with IV bolus dosing, and clearance (CL) as the model parameter of interest

predicted concentration: $C_i(t) = \frac{D}{V} \cdot e^{- \frac{\text{CL}_i}{V} t}$; assume $\text{CL}_i = \text{CL} \cdot e^{\eta{i}}$, where $\eta_i \sim \mathcal{N}(0, \omega^2)$

Observed data: $y_{ij} = C_i(t_{ij}) + \epsilon_{ij}, \quad \epsilon_{ij} \sim \mathcal{N}(0, \sigma^2)$

Problem: the function $C_i(t)$ is **nonlinear in** $\eta_i$, so the likelihood:

$L(\theta) = \prod_i \int \prod_j p(y_{ij} \mid \eta_i) p(\eta_i) d\eta_i$ 

is not analytically solvable.

To make this tractable, **FO linearizes** the model $f(t, \eta_i) = \frac{D}{V} e^{- \frac{\text{CL} \cdot e^{\eta_i}}{V} t}$ **around** $\eta_i = 0$

Let: $f(t, \eta_i) \approx f(t, 0) + \left. \frac{\partial f(t, \eta)}{\partial \eta} \right \vert_{\eta=0} \cdot \eta_i$

### Step 1: Evaluate $f(t, 0)$

$f(t, 0) = \frac{D}{V} e^{- \frac{\text{CL}}{V} t}$

### Step 2: Compute derivative w.r.t. $\eta$

Recall $\text{CL}_i = \text{CL} \cdot e^{\eta_i}$, so:

$\frac{df}{d\eta} = \frac{d}{d\eta} \left( \frac{D}{V} \cdot e^{- \frac{\text{CL} e^{\eta}}{V} t} \right)
= \frac{D}{V} \cdot e^{- \frac{\text{CL} e^{\eta}}{V} t} \cdot \left( -\frac{\text{CL} t}{V} \cdot e^{\eta} \right)$

Evaluated at $\eta = 0$:

$\left. \frac{df}{d\eta} \right\vert_{\eta=0} = \frac{D}{V} \cdot e^{- \frac{\text{CL}}{V} t} \cdot \left( -\frac{\text{CL} t}{V} \right)$

### Step 3: Linearized function

So, $f(t, \eta_i) \approx \frac{D}{V} \cdot e^{- \frac{\text{CL}}{V} t} -\left( \frac{D}{V} \cdot e^{- \frac{\text{CL}}{V} t} \cdot \frac{\text{CL} t}{V} \right) \cdot \eta_i$ is now **linear in** $\eta_i$, so:

You can compute the integral over $\eta_i$ analytically

## how about FOCE method?

The difference between FO and FOCE is that FO uses $\eta_i = 0$ as the linerization point, while FOCE $\eta_i = \hat{\eta}_i$ ($\hat{\eta}_i$ is the individual-specific conditional mode)

**First-Order (FO):** $f(t, \eta_i) \approx f(t, 0) + J_i(0) \cdot \eta_i$

**FOCE:**

1. For each individual, plug in current population parameters $\theta$, $\Omega$ and optimize $\eta_i$ by maximizing likelihood (This step requires an inner optimization loop per subject, hence “conditional” estimation)
    
    $\hat{\eta}_i = \arg \max[\log p(y_i \mid \eta_i, \theta) + \log p(\eta_i)]$
    
2. Linearize the model around $\hat{\eta}_i$ to optimize $\theta$, $\Omega$, $\sigma^2$ (population-level optimization)
    
     $f(t, \eta_i) \approx f(t, \hat{\eta}_i) + J_i(\hat{\eta}_i) \cdot (\eta_i - \hat{\eta}_i)$
     
## when to consider Laplace method

**The Laplace method** is the same as FOCE except that it uses a second derivative assessment of the variance of the joint density with respect to the ETAs. More accurate, more computational cost.

# 3. Expectation maximization (EM) algorithm

1. General idea of EM algorithm
    
    Complete data log-likelihood: $\log p(y, \eta \mid \theta, \omega, \sigma^2) = \sum_{i=1}^{N} \left[ \log p(y_i \mid \eta_i, \theta, \sigma^2) + \log p(\eta_i \mid \omega^2) \right]$
    
    1. E-step: For the E-step, we compute the **expected value** of this log-likelihood under the distribution of $\eta_i$ given data and current parameters (optimize a **tractable expected surrogate** instead)
        
        $Q(\theta, \omega^2, \sigma^2 \mid \theta^{(k)}, \omega^{2(k)}, \sigma^{2(k)}) = \mathbb{E}_{\eta \mid y, \theta^{(k)}, \omega^{2(k)}, \sigma^{2(k)}} \left[ \log p(y, \eta \mid \theta, \omega, \sigma^2) \right]$
        
    2. M-step: Given the approximate $Q$, maximize w.r.t. $\theta, \omega^2, \sigma^2$
    3. Repeat until convergence.
2. Monte Carlo approach to implement the E step
    - replace this **intractable expectation** Q with a **Monte Carlo average** based on M samples of $\eta_i$ from the posterior at the current parameter estimates. Then plug it to approximate true expectation. As the sample size N gets large, we can better approximate normal distribution, with mean of true expectation, and variances of variances of the original function/N.
        
        Sample $\eta_i^{(m)} \sim p(\eta_i \mid y_i, \theta^{(k)})$ and compute $Q(\theta) \approx \frac{1}{M} \sum \log p(y, \eta^{(m)} \mid \theta)$
        
        more samples improve convergence stability
        
    - importance sampling: use a new distribution q(x), so that the variances gets smaller. To take samples in the important region.
3. Difference between IMP EM and Complex nonlineaer models and/or many parameters are summarized in the table below

    | Feature | **Monte Carlo Importance Sampling EM (IMP)** | **MCMC SAEM** |
    | --- | --- | --- |
    | **NONMEM Option** | `METHOD=IMP` | `METHOD=SAEM` |
    | **Sampling method** | **Importance Sampling**: weighted samples from a proposal distribution | **MCMC (typically Metropolis-Hastings)**: samples from posterior using Markov chain |
    | **E-step Approximation** | Importance-weighted average of complete-data log-likelihood | Running average over MCMC samples (stochastic approximation) |
    | **Initialization** | Usually requires **well-tuned proposal distribution** | Robust to poor initial values (thanks to burn-in) |
    | **Convergence behavior** | More **deterministic**, but sensitive to sampling efficiency | **Stochastic**, but more robust for complex posteriors |
    | **Accuracy** | Can be **very accurate**, but depends on proposal quality | **Highly accurate** even in highly nonlinear models |
    | **Performance on complex models** | Less robust for strongly nonlinear or multimodal posteriors | Better for **highly nonlinear**, multimodal models |
    | **Computational cost** | Moderate | High (many iterations needed) |
    | **Variability in estimates** | Lower (deterministic samples) | Higher (due to stochastic updates) |
    | **Supports likelihood estimation (e.g., OFV)** | Yes | Yes (via Monte Carlo integration) |

# 4. Summary table: comparison of different estimation algorithms
**from top down: more accurate, better to deal with nonliearity, more computational expensive, greater incidence of success**

| Algorithm | Developed time | general idea | applications to models | applications to data | pros | cons |
| --- | --- | --- | --- | --- | --- | --- |
| FO | 1970s | 1st order partial derivative |  |  |  |  |
| Iterative two stage (ITS). Approximate EM method MCMC Bayesian analysis | 1984 | approximate, deterministic EM method | Rapid, exploratory method | Rich data |  |  |
| FOCE | 1992 | 1st order partial derivative | moderately nonlinear models; may fail for highly nonlinear models | Rich and semi-rich continuous data | Requires less fuss with adjusting options than SAEM or IMP | Convergence problem; Does not handle full OMEGA blocks as easily as EM |
| Laplace | 1992 | 2nd order partial derivative | Highly nonlinear models (e.g., turnover, Emax, logit/probit) | Large IIV (CV > 100%) |  |  |
| Monte Carlo importance sampling (IMP) expectation maximization (EM) | 2000s | EM+Monte Carlo sampling | Complex nonlineaer models and/or many parameters  | Sparse (fewer data points per subject than ETAs to be estimated) or rich data | Can handle full OMEGA blocks | less precise than FOCE because the results have stochastic variability |
| Markov Chain Monte Carlo (MCMC) stochastic approximation expectation maximization (SAEM) | 2000s | EM+Mone Carlo sampling | Complex nonlineaer models and/or many parameters  | Categorical data; Very sparse, sparse, or rich data | Can handle full OMEGA blocks | less precise than FOCE because the results have stochastic variability |
| Maximal likelihood expectation maximization (MLEM) in Adapt 5, **another way of IMP** | 1997 | EM+Monte Carlo sampling |  |  |  |  |