# Unifying probabilistic models for time-frequency analysis

https://arxiv.org/abs/1811.02489

In our paper we show equivalence between probabilistic time-frequency models (e.g. the probabilistic phase vocoder) and Spectral Mixture Gaussian processes. Therefore this code serves 3 novel purposes:

- Providing an easy way to construct more complex probabilistic time-frequency models by swapping in different kernel functions.

- Converting Spectral Mixture GPs to state space form so we can apply Kalman smoothing for efficient inference that scales linearly in the number of time steps.

- Hyperparameter tuning in spectral mixture GPs via a maximum likelihood approach in the frequency domain (Bayesian spectrum analysis).




matlab/ folder contains the code and example scripts.


matlab/experiments/ folder allows you to rerun the missing data synthesis experiments from the paper and produce the plots.


matlab/prob_filterbank folder contains Richard Turner's standard probabilistic time-frequency analysis code.



#### Reference:
```
@article{wilkinson2018unifying,
    author = {Wilkinson, William J. and Andersen, Michael Riis and Reiss, Joshua D. and Stowell, Dan and Solin, Arno},
    title = {Unifying probabilistic models for time-frequency analysis},
    journal = {ArXiv},
    volume = {1811.02489},
    year = 2018,
    }
```
