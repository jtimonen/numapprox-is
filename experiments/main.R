
# Requirements
library(cmdstanr)
library(posterior)
library(bayesplot)
library(checkmate)
library(loo)
library(stats)
library(outbreaks)
library(scales)
library(ggplot2)

# Options
stan_opts <- list(
  seed = 233,  # random seed for Stan
  sig_figs = 18 # number of significant figures to store in floats
)

# R functions
source("functions.R")
source("setup_sir.R")

code_parts <- setup_stancode_sir(solver="rk45")
sir_models <- create_cmdstan_models(code_parts)



prior_fit <- sir_models$prior$sample(
  sig_figs = stan_opts$sig_figs, 
  seed = stan_opts$seed,
  refresh = 1000
)
print(prior_fit)
prior_draws <- prior_fit$draws(c("beta", "gamma", "phi_inv"))


These draws can be passed to the second Stan model
```{r print_stan_model_sim}
sir_models$sim
```
that can be called in the `generate_quantities` mode to 
solve the system using the prior parameter draws. The solution provided by a numerical method is always an approximation to the
true solution. For the ODE solvers in Stan, the control arguments `rel_tol` 
and `abs_tol` (relative and absolute tolerance) have an effect on the accuracy
of this solution. In addition,
`max_num_steps` affects whether a solution can be obtained at all. Since we are
now doing the solves in `generated quantities`, we must set `max_num_steps`
large enough so that the solver never fails. However, this is not a problem at
this point since don't need gradients of the ODE solutions so solving the
system is fairly easy.

For now, we set some values for the solver arguments and plot the obtained
solutions.

```{r sir_sim_model, message=FALSE}
solver_args_default <- list(
  rel_tol = 1e-6,
  abs_tol = 1e-6,
  max_num_steps = 1e6
)
# Simulate solutions using prior draws
prior_sim <- simulate(sir_models$sim, prior_draws, data_list, 
                      solver_args_default, stan_opts)

# Plot solutions
title <- "Solutions using prior param draws"
plot_sir_example_solutions(prior_sim, data_list, thin=10, main=title)
```

## Analyzing the numerical method
How can we know that the approximate solution given by the numerical solver
is accurate enough? We approach this question by repeating the prior ODE
solving with different values of `rel_tol` and `abs_tol`.

```{r sir_test_prior_tols_time, results = FALSE}
TOLS <- 10^seq(-12,-4) # put less values here to knit the case study faster
atols <- TOLS
rtols <- TOLS
sims <- simulate_many(sir_models$sim, prior_draws, data_list, stan_opts,
                      atols, rtols, solver_args_default$max_num_steps)
```

We plot the mean and maximum absolute error compared to the solution given
by the most strict tolerances.

```{r sir_sol_errors, fig.width=8.2, fig.height=4.5}
mean_abs_sol_error <- compute_sol_errors(sims$sims, "mean")
max_abs_sol_error <- compute_sol_errors(sims$sims, "max")
plot_sim_errors(atols, rtols, mean_abs_sol_error)
plot_sim_errors(atols, rtols, max_abs_sol_error)
```

We see that lower tolerances generally produce more accurate results. Finally
we also plot the error in log likelihood, compared to the
one obtained with the strictest tolerances.

```{r sir_loglik_errors, fig.width=8.2, fig.height=4.5}
mean_abs_loglik_error <- compute_loglik_errors(sims$log_liks, "mean")
max_abs_loglik_error <- compute_loglik_errors(sims$log_liks, "max")
plot_sim_errors(atols, rtols, mean_abs_loglik_error, log=FALSE)
plot_sim_errors(atols, rtols, max_abs_loglik_error, log=FALSE)
```

Note that here `abs_tol` dominates because we are adding it to the solution
before likelihood computation to avoid negative values for the mean parameter,
which are not allowed. Negative solutions are a numerical artefact, because
they should not be possible with our ODE system with positive parameters.


## Computational challenges

Why not then always set the tolerances to very low values? At this point
we plot the simulation runtime for each pair of tolerances in the
previous experiment.

```{r print_res, fig.width=8.2, fig.height=4.5}
plot_sim_times(atols, rtols, sims$times)
```

We see that really small values for the tolerances start to require much
more computation time. This problem will become more pronounced when
we start to actually do posterior sampling.

When performing posterior inference for the parameters $\theta$
(and possible other model parameters, note that also $\textbf{y}_0$ can be a
parameter) using MCMC methods, the system needs to be solved numerically on each log posterior probability evaluation. And the number of these evaluations is
much higher than the number of obtained (effective) draws. 
Furthermore, when using a gradient
based MCMCM algorithm like Stan does, also _sensitivities_, i.e. 
derivatives of the solution w.r.t. parameters, are needed. This is done
in Stan by numerically solving an extended ODE system consisting of the
sensitivities and the original dynamical variables. Solving this larger
system takes more computation and also the solved sensitivities need to
achieve the set tolerances. So in order to be able to perform the model inference
in a reasonable time, we don't want to set the tolerances unnecessarily low.

As said before, the solution provided by a numerical method
is always an approximation to the true solution. 
This means that the computed posterior probability density and its
gradient are also approximations and the whole MCMC inference is biased
to some degree. So we are performing Bayesian inference for the parameters, the error in the solutions propagates also to the statistical properties of the posterior draws, in ways that cannot be predicted beforehand.

How could this error be studied, and how can be sure that our inference
results are correct? If the model parameters are fixed,
we can verify that the solution at all points is
appropriately close to a more accurate reference solution. 
There is a vast amount of classical literature that studies this error, and the ODE solvers built into
Stan are these classical methods that obtain estimate their error and adapt
their step size so that a given tolerance is achieved. However, we are doing statistics, and so we need to
know that the solution is accurate enough across all relevant parts of
parameter space. Additionally, it is not known beforehand where the 
"relevant parts of parameter space" are!
  
  The problem of validating the use of a numerical method for a Bayesian model
is therefore significantly more complicated than in the classical numerical
analysis world. One idea is to run sampling repeatedly, gradually increasing
the solver accuracy, until the posterior inference
results don't change anymore. But this might not
be computationally very attractive, as posterior sampling until $S$ parameter
draws are obtained is
computationally much more demanding than just solving the ODE system $S$
times, for the reasons stated above.

These problems  occur also
with any other type of model that requires numerically solving an implicitly defined function or variable, but in this case study we focus on ODEs.

The point of this case study is to show how to obtain reliabile inferences
efficiently. This involves the use of Pareto-Smoothed Importance Sampling
(PSIS) [@yao2018; @vehtari2019].

# Workflow

Let $M$ be the model for which we would like to perform inference, but which we
cannot evaluate since the likelihood is defined implicitly through an
ODE system that is not analytically tractable. MCMC inference for $M$
can be seen actually as inference for another model $M_{high}$, which
is the same model as $M$ but using a numerical solver, and can therefore be
evaluated.

Our workflow addresses the problem of defining the high-precision numerical
method in $M_{high}$ so that $M_{high}$ can trusted to have essentially the
same posterior as $M$. We define a way to perform inference for $M_{high}$ 
without needing to compute gradients or HMC trajectories for it. This involves
another model $M_{low}$, which is again the same model, except that $M_{low}$
uses a cheaper and less accurate numerical methods (or just looser tolerances
and/or coarser discretization) to compute the required ODE solutions,
and is therefore potentially faster to fit. The posterior densities are denoted
$p_{low}$ and $p_{high}$, respectively.

To understand how PSIS comes into play, we must first discuss importance
sampling. If we want to compute expectations with the high precision model, we
can take draws from the low precision models and reweight these according to the
importance weights $\frac{p_{high}}{p_{low}}$. If these models are too
different, then the reweighting will produce noisy estimates that are not
useful. PSIS and particularly the Pareto $k$-diagnostic (denoted $\hat{k}$), is
the tool that tells us when we can or cannot rely on the importance weights. If
$\hat{k} < 0.5$ we are safe to do the importance sampling, if $\hat{k} < 0.7$
the importance sampling will start to converge more slowly, and if
$\hat{k} > 0.7$ the importance sampling estimates are unreliable. For simplicity
we will only consider the $\hat{k} < 0.5$ threshold.

Ideally, $M_{high}$ would involve a numerical method that we can trust
completely in all parts of the parameter space so that, as long as
$\hat{k} < 0.5$, importance weights can be used to reweight the low precision
approximation $p_{low}$ to the high precision approximation $p_{high}$. We can
think of $M_{high}$ as a reference model, because it is the baseline to which
we compare. It is difficult in practice to know if a given model is a good
reference model in all parts of parameter space, due to curse of dimensionality
and the fact that analysed system can have different properties in different
parts of the parameter space. For example, ODEs can qualitatively change their
behaviour as a function of parameters (bifurcation), or become stiff or
chaotic in some parameter regions. Accuracy can be checked at a given set of
parameters fairly easily, but not over a high dimensional parameter space.
Under these conditions it is necessary to compromise to develop a reference
model that works only over a range of parameter space, but even then it is hard
to know *a priori* what range that is.

We propose the following workflow:

::: {#workflow .box}
1. Generate draws from $p_{low}$.
2. Tune the numerical method in $M_{high}$ so that it is reliable at
these draws. All application specific knowledge and classical numerical
analysis can be used here.
3. Compute importance weights $\frac{p_{high}}{p_{low}}$ and the $\hat{k}$ 
diagnostic. If $\hat{k} > 0.5$, raise precision of
the numerical method in $M_{low}$ and go back to step 1.
4. Resample draws of $M_{low}$, using importance weights to get draws from
$M_{high}$, and therefore essentially draws from $M$.
:::

The next two sections of this case study outline how to apply this workflow to
do fast but reliable inference for an ODE model using a built-in Stan solvers.
The importance sampling diagnostics are handled with the
[`loo`](https://mc-stan.org/loo) package and the resampling is handled with the
[`posterior`](https://github.com/stan-dev/posterior) package.


## SIR Example: going posterior
We consider again the SIR example and the influenza data presented earlier.
The parameters $\beta$ and $\gamma$
will be estimated from time series observations of the number of infected
people (I). For the purposes of this case study, the goal is to use a very low
precision ODE solver to do inference and then check it afterwards against
a high precision solver. This is useful in practice if sampling with the
high precision solver itself would take an inordinate amount of time.

We will now need the third Stan model that we built earlier
```{r print_posterior_model}
sir_models$posterior
```
for posterior sampling.

### Step 1. Generate draws from $p_{low}$.

We start with `abs_tol=rel_tol=1e-4`. Generally, it has been observed to be
a bad idea to use
tolerances larger than this, because the gradients start to be badly wrong
and everything seems to break down so that sampling actually can become
slower than with stricter tolerances.

```{r sample_posterior, message=FALSE}
solver_args_sample <- list(
  abs_tol = 1e-4,
  rel_tol = 1e-4,
  max_num_steps=1e6
)
fit <- sample_posterior(sir_models$posterior, data_list, 
                                  solver_args_sample, stan_opts, refresh=1000,
                                  init=0)
```

```{r print_post}
print(fit$time())
print(fit)
post_draws <- fit$draws(c("beta", "gamma", "phi_inv"))
```

These draws are again passed again to `sir_models$sim` and we vizualize the
simulated model dynamics.

```{r sir_sim_posterior, message=FALSE}
# Simulate and plot generated solutions
post_sim <- simulate(sir_models$sim, post_draws, data_list,
                          solver_args_sample, stan_opts)

title <- "Solutions using posterior param draws"
plot_sir_example_solutions(post_sim, data_list, thin=10, main=title)
```

### Step 2. Tune the numerical method in $M_{high}$ so that it is reliable at these draws.

We write a function that does this tuning for us, by gradually strictening
the tolerances until the solutions effectively don't change anymore 
(maximum absolute error changes less than `tuning_tol`).

```{r sim_posterior_high, message=FALSE}
TOLS <- 10^seq(-4,-12) # could be just halving for example
tuning_tol <- 0.0001
tuning <- tune_solver(TOLS, sir_models$sim, post_draws, data_list, stan_opts, 
                      solver_args_sample$max_num_steps, tuning_tol)

# Print the final tolerance
print(tuning$last_tol)
```


### Step 3. Compute importance weights and Pareto-$k$

```{r importance_weights}
is <- use_psis(post_sim, tuning$last_sim)
print(is)
print(is$diagnostics$pareto_k)
```

### Step 4. Importance resampling

```{r resampling}
post_draws_resampled <- posterior::resample_draws(post_draws,
                                                  weights=is$log_weights)

# Now we could use the resampled draws for whatever...
```


### Comparing runtimes

The total time needed to run the workflow comes from the inital posterior
sampling and the time required for running `tune_solver()` (other computations
                                                            are negligible). Had we needed to return to Step 1, the total time would have
of course increased.

```{r timing}
workflow_time <- fit$time()$total + sum(tuning$times)
cat("Total workflow time was", workflow_time, "seconds.\n", sep=" ")
```

Let's see how long it would have taken to do the original posterior
sampling with the tolerances that we found necessary for the reference method.


```{r reference_fit, message=FALSE}
solver_args_refsample <- list(
  abs_tol = tuning$last_tol,
  rel_tol = tuning$last_tol,
  max_num_steps=solver_args_sample$max_num_steps
)
fit_ref <- sample_posterior(sir_models$posterior, data_list, 
                                  solver_args_refsample, stan_opts, refresh=1000,
                                  init=0)
print(fit_ref$time()$total)
```
