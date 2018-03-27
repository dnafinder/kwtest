# kwtest
Kruskal-Wallis test for the non parametric ANOVA<br/>
In statistics, the Kruskal–Wallis one-way analysis of variance by ranks (named
after William Kruskal and W. Allen Wallis) is a non-parametric method for
testing equality of population medians among groups. It is identical to a
one-way analysis of variance with the data replaced by their ranks. It is an
extension of the Mann–Whitney U test to 3 or more groups. Since it is a
non-parametric method, the Kruskal–Wallis test does not assume a normal
population, unlike the analogous one-way analysis of variance. However, the
test does assume an identically-shaped and scaled distribution for each group,
except for any difference in medians. The exact distribution of the
Kruskal-Wallis statistics is very time, space and memory expensive, so it is
approximated using other distribution.
The MatLab function KRUSKALWALLIS only uses the chi-square distribution that
is the most conservative (and this means that it accepts the H0 hypothesis more
than you want). This function computes also the F, Beta and Gamma
approximations (the F distribution is the less conservative and this means that
it refuses the H0 more than you want), giving a more informative output (if you
want, see http://www.jmu.edu/assessment/JPM%20AERA%20SP%2008.pdf). If you
believe that the p-value you choose is smaller than your cut-off (usually
0.05), you can use my function Dunn-Sidak to isolate the differences among
groups (http://www.mathworks.com/matlabcentral/fileexchange/12827).

 Syntax: 	STATS=kwtest(X)
      
     Inputs:
           X - data matrix (Size of matrix must be n-by-2; data=column 1, group=column 2). 

     Outputs:
           - Statistics of each group (samples, median, sum of ranks and mean rank)
           - Correction factor for ties and H statistics
           - Chi-square approximation
           - F approximation
           - Beta approximation
           - Gamma approximation
        If STATS nargout was specified the results will be stored in the STATS
        struct.

To cite this file, this would be an appropriate format:
Cardillo G. (2009). KWTEST: Kruskal-Wallis non parametric test for ANalysis Of VAriance
http://www.mathworks.com/matlabcentral/fileexchange/25860
