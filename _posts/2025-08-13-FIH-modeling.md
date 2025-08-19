---
layout: post
title: PK Modeling Journey for a First-in-Human (FIH) Study
date: Aug-13-2025
categories: pharmacometrics
tags: pharmacometrics
permalink: 
banner: ./pics/FIH_trial.png
math: true
---

FIH programs usually generate single ascending dose (SAD) and multiple ascending dose (MAD) data. 
One early and important question is to model PK using the observed data to interpolate reliable PK predictions for other dosing regimens. 
This post lays out a simple decision path and the modeling approach. 

For an oral drug, one typical scenario is to have the SAD cohorts with a fasted condition and MAD cohorts with a fed condition. Under fed conditions, food effects can come into play.

First let us define dose-normalized exposure metrics:

  $dnAUC_{inf} = AUC_{inf} / Dose$
  
  $dnC_{max} = C_{max} / Dose$
  
  If we see similar dose-normalized exposure metrics across different dose levels, then it is suggesting linear PK; increasing trend with dose ≈ supraproportional; decreasing with dose ≈ subproportional.

## 1. Start with SAD (fasted) PK modeling: check dose proportionality using $dnAUC_{inf}$ and $dnC_{max}$

If approximately linear in SAD, a fixed F1=1 can be used to estimate apparent clearance ($CL/F$) and apparent volume terms for a compartmental model. First-order absorption model can be tried first as the initial model. 

If subproportional in SAD (i.e., lower $dnC_{max}$), we can suspect it is absorption-limited, especially for BCS class II drug. 

If supraproportional in SAD, maybe the clearance is saturated at higher doses, leading to higher dose-normalized exposure. 

Flag which direction the deviation goes and carry the hypothesis into MAD.

## 2. Move to MAD (fed): check steady-state proportionality using $dnAUC_{\tau, SS}$

If linear across MAD cohorts, a constant F1=$\theta_1$ can be estimated. Here, $\theta_1$ represents the relative bioavailabilty at fed compared to the fasted condition.

If $dnAUC_{tau, SS}$ increases with dose, two common mechanisms need to explored:

* Possibility 1 — Nonlinear clearance (capacity-limited elimination)
* Possibility 2 — Increased F1 under fed conditions (food effect)

These can confound each other. Here’s how to separate them.

## 3. Exploration & modeling implementation

### Possibility 1: Nonlinear clearance (capacity-limited elimination)

* **What you should see:** dose-dependent increase in $t_{₁/₂}$ (because effective clearance decreases with concentration). $t_{1/2}=\frac{\ln 2 \cdot V}{CL}$
* **Modeling approch:** use Michaelis-Menten elimination alone **or** linear clearance plus Michaelis-Menten elimination

### Possibility 2: Food effect on F1

* **What you should see:** $t_{₁/₂}$ similar across MAD cohorts; across MAD cohorts: higher exposure after the first dose with higher doses (use $AUC_{last}$ or $AUC_{\infty}$), because food effect can show right after the first dose; $AUC_{\infty,Day1}\approx AUC_{\tau,ss}$
* **Modeling approch:**

  * Treat SAD (fasted) as reference.
  * For MAD (fed), we can model dose-dependent food effect (e.g., Emax form). 

    \[F_{MAD}=FE \times \Bigl(1 + \frac{F_{\max}\cdot Dose}{Dose+FD_{50}}\Bigr)\]

    * **FE:** baseline food effect vs fasted at a minimal dose.
    * **$F_{max}$:** maximum fractional increase in F1 due to dose under fed conditions.
    * **$FD_{50}$:** dose producing half of $F_{max}$ under fed conditions (same units as Dose).
    * In this equation, we fix $F_{SAD}$=1 as a constant. This assumption can be adjusted based on SAD PK linearity. 

## 4. Math/PK signatures

 For linear PK, single-dose $AUC_{\infty}$ and $AUC_{\tau,ss}=\frac{F\cdot \text{Dose}}{CL}$ equals the same (see Blood Levels of Drug at the Equilibrium State after Multiple Dosing, 1965):

  \[AUC_{\infty,\text{single dose}}=\frac{F\cdot \text{Dose}}{CL}\]

  ⇒ At a given dose, **$AUC_{\infty,\text{single dose}}\approx AUC_{\tau,ss}$** if kinetics are linear and F, CL don’t change across occasions.


## 5. Sanity checks

* %AUC extrapolation for single-dose $AUC\_{\infty}$: keep ≤ 20% or flag as unreliable.
* Steady state achieved? Confirm ≥ 4–5 terminal half-lives before relying on $AUC_{\tau,ss}$.
* The absorption or disposition mechanisms included in the model should be checked against results from preclinical study (e.g., BCS class, ADME pathways) to have a consistent mechanistic understanding.



