function STATS = kwtest(x,varargin)
%KWTEST Kruskal-Wallis test for the non parametric ANOVA
% In statistics, the Kruskal–Wallis one-way analysis of variance by ranks
% (named after William Kruskal and W. Allen Wallis) is a non-parametric
% method for testing equality of population medians among groups. It is
% identical to a one-way analysis of variance with the data replaced by
% their ranks. It is an extension of the Mann–Whitney U test to 3 or more
% groups.
%
% The MATLAB function KRUSKALWALLIS only uses the chi-square distribution,
% which is the most conservative approximation of the Kruskal-Wallis
% statistic. This function also computes the F, Beta and Gamma
% approximations, providing a more informative output.
%
% Syntax:
%     STATS = kwtest(x)
%     STATS = kwtest(x,'Display',DISPLAY)
%
% Inputs:
%     X - N-by-2 data matrix:
%         X(:,1) = data
%         X(:,2) = group labels (whole numbers)
%
% Name-Value Pair:
%     'Display' - true/false (default true). If true, results are printed.
%
% Outputs:
%     - Statistics of each group (samples, median, sum of ranks and mean rank)
%     - Correction factor for ties and H statistics
%     - Chi-square approximation
%     - F approximation
%     - Beta approximation
%     - Gamma approximation
%
% Example:
% x=[7.79 9.16 7.64 10.28 9.12 9.24 8.40 8.60 8.04 8.45 9.51 8.15 7.69 ...
% 8.84 9.92 7.20 9.25 9.45 9.14 9.99 9.21 9.06 8.65 10.70 10.24 8.62 ...
% 9.94 10.55 10.13 9.78 9.01; ...
% repmat(1,1,13) repmat(2,1,9) repmat(3,1,9)]'
%
%     STATS = kwtest(x)
%
% Created by Giuseppe Cardillo
% Email: giuseppe.cardillo.75@gmail.com
% GitHub: https://github.com/dnafinder
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2009). KWTEST: Kruskal-Wallis non parametric test for ANalysis
% Of VAriance.
% Available at: https://github.com/dnafinder/kwtest

% Input Error handling
p = inputParser;
addRequired(p,'x',@(y) validateattributes(y,{'numeric'}, ...
    {'real','finite','nonnan','nonempty','ncols',2}));
addParameter(p,'Display',true, ...
    @(d) (islogical(d) && isscalar(d)) || ...
         (isnumeric(d) && isscalar(d) && (d==0 || d==1)));
parse(p,x,varargin{:});
x = p.Results.x;
Display = logical(p.Results.Display);
clear p

% Validate group labels
assert(all(x(:,2) == fix(x(:,2))), ...
    'kwtest:InvalidGroupLabels', ...
    'All elements of column 2 of input matrix must be whole numbers.');
assert(all(x(:,2) > 0), ...
    'kwtest:InvalidGroupLabels', ...
    'Group labels must be positive whole numbers.');

tr = repmat('-',1,80);

% Total observations
N = size(x,1);

% Rank all observations (pooled)
[R,T] = tiedrank(x(:,1));

% Groups
groups = x(:,2);
ug = unique(groups);
k = numel(ug);

% Compute group statistics
KW = zeros(4,k); % rows: Samples, Median, Ranks_sum, Mean_rank

for i = 1:k
    gval = ug(i);
    idx = (groups == gval);
    xi = x(idx,1);
    ri = R(idx);

    KW(1,i) = sum(idx);          % Samples
    KW(2,i) = median(xi,1);      % Median
    KW(3,i) = sum(ri);           % Ranks sum
    KW(4,i) = KW(3,i) / KW(1,i); % Mean rank
end

GroupTable = table(ug, KW(1,:)', KW(2,:)', KW(3,:)', KW(4,:)', ...
    'VariableNames',{'Group','Samples','Median','Ranks_sum','Mean_rank'});

% Mean rank under H0
Rbar = (N+1)/2;

% Variability measure
D = sum(KW(1,:) .* (KW(4,:) - Rbar).^2);

% Uncorrected H
Hbiased = 12 * D / N / (N+1);

% Tie correction
if T == 0
    CF = 1;
    H = Hbiased;
else
    CF = 1 - 2*T/N/(N^2-1);
    H = Hbiased / CF;
end

% Approximations
df = k - 1;

% Chi-square approximation (most conservative)
P1 = 1 - chi2cdf(H, df);

% F approximation (less conservative)
dfd = N - df; % kept for backward compatibility with historical output style
F = ((dfd + 1) * H) / ((k - 1) * (N - 1 - H));
P2 = 1 - fcdf(F, df, N - k - 1);

% Beta/Gamma parameters
m = k - 1; % expected value
s2 = 2*df ...
    - 2*(3*k^2 - 6*k + N*(2*k^2 - 6*k + 1)) / (5*N*(N+1)) ...
    - (6/5 * sum(1 ./ KW(1,:)));
eta = (N^3 - sum(KW(1,:).^3)) / N / (N+1);

% Beta approximation
B = H / eta;
alpha1 = m * ((m*(eta-m) - s2) / (eta*s2));
beta1  = alpha1 * ((eta - m) / m);
P3 = 1 - betacdf(B, alpha1, beta1);

% Gamma approximation
alpha2 = m^2 / s2;
beta2  = s2 / m;
P4 = 1 - gamcdf(H, alpha2, beta2);

% Display results
if Display
    disp('KRUSKAL-WALLIS TEST')
    disp(tr)
    disp(GroupTable)

    fprintf('Correction factor for ties: %0.4f\tH: %0.4f\n', CF, H)
    disp(tr); disp(' ')

    fprintf('Chi-square approximation (the most conservative)\n')
    disp(tr)
    disp(table(H, df, P1, ...
        'VariableNames',{'Chi_square','df','p_value'}))

    fprintf('F-statistic approximation (the less conservative)\n')
    disp(tr)
    disp(table(F, df, dfd, P2, ...
        'VariableNames',{'F','df_num','df_denom','p_value'}))

    fprintf('Beta distribution approximation\n')
    disp(tr)
    disp(table(B, alpha1, beta1, P3, ...
        'VariableNames',{'B','alpha','beta','p_value'}))

    fprintf('Gamma distribution approximation\n')
    disp(tr)
    disp(table(H, alpha2, beta2, P4, ...
        'VariableNames',{'G','alpha','beta','p_value'}))
end

% Build output structure
STATS = struct();

STATS.N = N;
STATS.k = k;
STATS.groups = ug;
STATS.GroupTable = GroupTable;

STATS.Hbiased = Hbiased;
STATS.CF = CF;
STATS.H = H;

STATS.Chi_square.chi2 = H;
STATS.Chi_square.df = df;
STATS.Chi_square.pvalue = P1;

STATS.F.F = F;
STATS.F.dfn = df;
STATS.F.dfd = dfd;
STATS.F.pvalue = P2;

STATS.Beta.m = m;
STATS.Beta.s2 = s2;
STATS.Beta.eta = eta;
STATS.Beta.B = B;
STATS.Beta.alpha = alpha1;
STATS.Beta.beta = beta1;
STATS.Beta.pvalue = P3;

STATS.Gamma.m = m;
STATS.Gamma.s2 = s2;
STATS.Gamma.G = H;
STATS.Gamma.alpha = alpha2;
STATS.Gamma.beta = beta2;
STATS.Gamma.pvalue = P4;
end
