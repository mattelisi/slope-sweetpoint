These functions implement an adaptive maximum-likelihood procedure that optimise stimulus placement for measurement of the slope of the psychometric function, using a minimum-variance criterion. That is, at each step these functions allow to compute the position of the two sweet points (see Green, 1990, J. Acoust.Soc. Am.) that minimise the expected variance in the slope estimate.

The current implementation is suited for standard 2AFC task, use a cumulative Gaussian psychometric function, and assumes that the probability of lapsing is equal for small vs large stimulus values. (In the computations the probability of lapsing is considered a free parameter.)

Things to do:
- make it Bayesian? i.e. use a prior and place some trials also at the sweetpoint for the threshold (instead of starting with few random trials). It may not be ideal for some applications because if the prior is mis-specified there is little or no gain at all in efficiency with respect to random. 
