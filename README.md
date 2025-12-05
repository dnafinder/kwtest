[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/kwtest&file=kwtest.m)

# kwtest

## 📘 Overview
kwtest is a MATLAB function implementing the Kruskal-Wallis one-way analysis of variance by ranks, a nonparametric alternative to classical one-way ANOVA.

The Kruskal-Wallis test evaluates whether three or more independent groups come from populations with equal medians. It does not assume normality and is based on ranks of the pooled observations.

Unlike MATLAB’s built-in kruskalwallis, which reports a chi-square approximation only, kwtest also provides additional approximations of the Kruskal-Wallis statistic to give a more informative view of significance:

- Chi-square approximation (most conservative)
- F approximation (less conservative)
- Beta approximation
- Gamma approximation

This can be useful when you want to understand how sensitive your inference might be under different classical approximations of the same statistic.

## ✨ Features
- Input format suitable for quick analysis: N-by-2 matrix [data, group]
- Works with three or more groups
- Supports non-consecutive group labels (positive integers)
- Computes group summaries:
  - sample size
  - median
  - sum of ranks
  - mean rank
- Tie correction for H statistic
- Four distributional approximations with p-values
- Script-friendly structured output
- Optional Display flag to control printing

## 📥 Installation
1. Download or clone the repository:
   https://github.com/dnafinder/kwtest

2. Add the folder to your MATLAB path:
      addpath('path_to_kwtest')

3. Verify installation:
      which kwtest

## ⚙️ Requirements
- MATLAB (any recent version)
- Statistics and Machine Learning Toolbox (required for:
  - tiedrank
  - chi2cdf
  - fcdf
  - betacdf
  - gamcdf)

## 📈 Usage

Basic usage:

    STATS = kwtest(x);

Silent mode (no console output):

    STATS = kwtest(x, 'Display', false);

## 🔢 Inputs

kwtest(x)  
kwtest(x, 'Display', DISPLAY)

- x
  - N-by-2 numeric matrix
  - x(:,1): data
  - x(:,2): group labels
  - group labels must be positive whole numbers
  - labels do not need to be consecutive

- 'Display'
  - true/false (default true)
  - controls printing of the results to the Command Window

## 📤 Outputs
kwtest returns a structure STATS containing:

- STATS.N  
  total number of observations

- STATS.k  
  number of groups

- STATS.groups  
  unique group labels in ascending order

- STATS.GroupTable  
  a table with:
  Group, Samples, Median, Ranks_sum, Mean_rank

- STATS.Hbiased  
  Kruskal-Wallis statistic before tie correction

- STATS.CF  
  tie correction factor (CF = 1 when no ties)

- STATS.H  
  tie-corrected Kruskal-Wallis H statistic

Approximations:

- STATS.Chi_square
  - chi2
  - df
  - pvalue

- STATS.F
  - F
  - dfn
  - dfd
  - pvalue

- STATS.Beta
  - m
  - s2
  - eta
  - B
  - alpha
  - beta
  - pvalue

- STATS.Gamma
  - m
  - s2
  - G
  - alpha
  - beta
  - pvalue

## 🧠 Interpretation
- Small p-values suggest that at least one group differs in median from the others.
- The chi-square approximation is typically the most conservative.
- The F approximation is typically less conservative.
- Beta and Gamma approximations provide alternative classical views of the same test statistic.

If your chosen approximation yields significance, you may consider a multiple-comparison post-hoc procedure to identify which groups differ.

## 📌 Example

Example dataset from the function help:

    x = [7.79 9.16 7.64 10.28 9.12 9.24 8.40 8.60 8.04 8.45 9.51 8.15 7.69 ...
         8.84 9.92 7.20 9.25 9.45 9.14 9.99 9.21 9.06 8.65 10.70 10.24 8.62 ...
         9.94 10.55 10.13 9.78 9.01; ...
         repmat(1,1,13) repmat(2,1,9) repmat(3,1,9)]';

    STATS = kwtest(x);

## 🧾 Citation
If you use kwtest in research, analysis, or publications, please cite:

Cardillo G. (2009). KWTEST: Kruskal-Wallis non parametric test for ANalysis Of VAriance.  
Available at: https://github.com/dnafinder/kwtest

## 👤 Author
Giuseppe Cardillo  
Email: giuseppe.cardillo.75@gmail.com  
GitHub: https://github.com/dnafinder

## 📄 License
The code is provided as-is, without any explicit warranty.  
Please refer to the repository for licensing details if a LICENSE file is present.
