# numapprox-is
An importance sampling framework for reliable and efficient inference in Bayesian models 
that require numerical approximations, such as solutions of implicitly defined functions. 
Ordinary differential equation (ODE) models are one example.

## Experiments

Code for reproducing the experiments of the [paper](https://onlinelibrary.wiley.com/doi/full/10.1002/sta4.614) is in the `experiments` subdirectory. Running
them requires the [odemodeling](https://github.com/jtimonen/odemodeling) R package,
which was developed for these experiments. Experiments were run using version 0.2.0 of *odemodeling*. 

## References

Timonen, J., Siccha, N., Bales, B., Lähdesmäki, H., & Vehtari, A. (2023). **An importance sampling approach for reliable and efficient inference in Bayesian ordinary differential equation models**. *Stat*, 12(1), e614. https://doi.org/10.1002/sta4.614

