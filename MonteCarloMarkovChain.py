import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal, norm


def prior(mean, covariance, data):
    if mean[0] < 1:
        return 0
    else:
        return 1


def acceptance_rule(log_posterior_new, log_posterior_old):
    if log_posterior_new > log_posterior_old:
        return True
    else:
        acceptance_ratio = min(1, (log_posterior_old / log_posterior_new))
        u = np.random.uniform()
        if u < 0.1: #acceptance_ratio:
            return True
        else:
            return False


def log_likelihood(mean, covariance, data):
    return np.sum(np.log(multivariate_normal(mean, cov=covariance).pdf(data)))


def mean_transition_model(mean, covariance):
    return multivariate_normal.rvs(mean=mean, cov=covariance)


def covariance_transition_model(covariance, data):
    cnt = 0
    accepted = False
    while not accepted:
        try:
            cov_perturbation = multivariate_normal.rvs(mean=covariance.reshape(-1), cov=[0, 0, 0, 0]).reshape((2, 2))
            multivariate_normal([1, 1], cov=cov_perturbation).pdf(data)
            accepted = True
        except:
            print(f"trying to find a covariance, attempt no: {cnt}", flush=True, end='\r')
            cnt += 1
            pass

    return cov_perturbation