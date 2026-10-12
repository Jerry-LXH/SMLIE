"""Small log-domain HMM, A[i,j] = p(z_t=j | z_{t-1}=i)."""
import numpy as np
from itertools import product
from math import lgamma

def lse(x, axis=None):
    return np.logaddexp.reduce(x, axis=axis)

def poisson_log_emissions(y, rates):
    y=np.asarray(y); rates=np.asarray(rates)
    return y[:,None]*np.log(rates)[None,:]-rates[None,:]-np.array([lgamma(float(v)+1) for v in y])[:,None]

def forward_backward(log_b, pi, A):
    log_b=np.asarray(log_b); T,K=log_b.shape
    if T<1 or np.shape(A)!=(K,K) or np.shape(pi)!=(K,):
        raise ValueError('Expected nonempty [T,K] emissions, [K] pi and [K,K] A')
    if np.any(np.asarray(A)<0) or np.any(np.asarray(pi)<0) or not np.allclose(np.sum(A,axis=1),1) or not np.isclose(np.sum(pi),1):
        raise ValueError('pi and A must be probability distributions')
    with np.errstate(divide='ignore'):
        lp=np.log(pi); la=np.log(A)
    alpha=np.empty((T,K)); alpha[0]=lp+log_b[0]
    for t in range(1,T): alpha[t]=log_b[t]+lse(alpha[t-1,:,None]+la,axis=0)
    ll=lse(alpha[-1])
    if not np.isfinite(ll): raise ValueError('Observation has zero or nonfinite probability')
    beta=np.zeros((T,K))
    for t in range(T-2,-1,-1): beta[t]=lse(la+log_b[t+1][None,:]+beta[t+1][None,:],axis=1)
    gamma=np.exp(alpha+beta-ll)
    xi=np.array([np.exp(alpha[t,:,None]+la+log_b[t+1][None,:]+beta[t+1][None,:]-ll) for t in range(T-1)]).reshape(T-1,K,K)
    filtering=np.exp(alpha-lse(alpha,axis=1)[:,None])
    return float(ll),gamma,xi,filtering

def brute_force(log_b,pi,A):
    T,K=log_b.shape; paths=np.array(list(product(range(K),repeat=T)))
    with np.errstate(divide='ignore'):
        scores=np.array([np.log(pi[z[0]])+sum(log_b[t,z[t]] for t in range(T))+sum(np.log(A[z[t-1],z[t]]) for t in range(1,T)) for z in paths])
    ll=lse(scores); w=np.exp(scores-ll)
    g=np.array([[w[paths[:,t]==k].sum() for k in range(K)] for t in range(T)])
    return float(ll),g,paths,w
