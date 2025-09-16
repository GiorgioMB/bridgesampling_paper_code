# Estimating Marginal Likelihoods Using Bridgesampling
This repository contains the code and data for a [paper](https://arxiv.org/abs/2508.14487) proposing diagnostic methods for bridgesampling. The experiments are divided into two main parts: analysis using posterior distributions from the `posteriorDB` project and a toy example exploring the effects of increasing the number of covariates in a generalized linear model (GLM).

## Repository Structure
- `posteriordb/`: Contains the main experiments that apply bridgesampling to various posterior distributions curated from the `posteriorDB`. This folder includes both R scripts and Stan models that are essential for replicating the findings and further exploration.
- `toy_example/`: Includes experiments designed to study the impact of an increasing number of covariates in a GLM on the efficiency and accuracy of the bridgesampling method. This section is helpful for understanding scalability and performance in simpler, controlled scenarios.

### Contents

- **R files**: Scripts for setting up the statistical models, executing the bridgesampling algorithm, and processing the results.
- **Stan files**: Stan models used to define the Bayesian models whose posteriors are being explored.
- **Data files**: Data files used in the experiments, stored in csv format.

## Getting Started

### Prerequisites

Ensure you have the following software and libraries installed:
- R (version 4.0 or higher)
- RStan (version 2.21 or higher)
- bridgesampling package in R
- CmdstanR (bridgesampling version)

### Installation

Clone the repository to your local machine using:

```bash
git clone https://github.com/GiorgioMB/bridgesampling_paper_code.git
cd bridgesampling_paper_code
```

#### Usage
Navigate to either of the experiment directories and run the R scripts provided. For example:
```bash
cd posteriordb
Rscript low_dim_gauss_mix.R
```
This will execute the analysis using the predefined model and data from `posteriorDB`.

## Contacts
For any queries related to the repository, please contact:

ext-giorgio.micaletto@aalto.fi
