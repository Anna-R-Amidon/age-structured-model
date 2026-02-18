# age-structured-model
Density Dependent Age-0 Survival of Sablefish (Anoplopoma fimbria): Age-Structured Modeling Assessment of Commercial Harvest Feasibility

This repository demonstrates parameter estimation for an age-structured population dynamics model under commercial harvest pressure. 
The project estimates early-life survival and density dependence using maximum likelihood methods and evaluates model fit to observed abundance-at-age data.

This repo is intended as a **portfolio example** showcasing:
- Statistical modeling of population dynamics  
- Likelihood-based parameter estimation  
- Reproducible analysis workflows  
- Scientific reporting and visualization


 ---

## Problem Statement

A managed fish population has experienced periodic commercial harvest followed by closures. Abundance-at-age (ages 1–4+) and age-specific harvest rates were 
monitored after the first closure. Survival and fecundity are known for ages 1–4, but early-life survival (age-0) and density dependence are unknown.

Age-0 survival is modeled as:

<sub>0</sub> = α · e<sup>−βᵈ N<sub>t</sub></sup>

where:
- **α** = baseline age-0 survival  
- **βᵈ** = strength of density dependence  
- **Nₜ** = total population abundance  

The objective is to estimate **α** and **βᵈ**, quantify uncertainty, and assess model fit.

---

## Methods 

**Methods**
- Age-structured population model  
- Maximum likelihood estimation  
- Optimization using mle2, SANN 
- Model validation via observed vs. predicted abundance-at-age  

**Tools**
- R  
- Reproducible scripts and version-controlled workflow  
- 'bbmle', 'tidyverse', 'ggplot2'
