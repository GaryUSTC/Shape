import numpy as np
from numpy import linalg as LA
from numpy.linalg import inv
import random
import matplotlib.pyplot as plt


# helper function
def decompose_kernel(M):
    [D, V] = LA.eig(M)
    return M, np.real(V), np.real(D)

def elem_sympoly(k_lambda,k):
    N = len(k_lambda)
    E = np.zeros((k+1, N+1))
    E[0,:] = 1
    for l in xrange(1,k+1):
        for n in xrange(1, N+1):
            E[l,n] = E[l, n-1] + k_lambda[n-1]*E[l-1, n-1]
    return E

def sample_k(k_lambda, k):
    E = elem_sympoly(k_lambda, k)
    i = len(k_lambda) - 1
    remaining = k-1
    S = np.zeros((k,1))
    # Note: here >= to fit matlab
    while remaining >= 0:
        if i == remaining:
            marg = 1
        else:
            marg = k_lambda[i]*E[remaining, i] / E[remaining+1, i+1]
        if random.random() < marg:
            S[remaining] = i
            remaining = remaining - 1

        # print remaining, i
        i=i-1
    return S

def sample_dpp(L, k):
    # Perform DPP sampling
    M, V, D = decompose_kernel(L)
    # DPP and k_DPP
    if k != 0:
        # do k-DPP
        v = sample_k(D, k)
    else:
        # do DPP
        D = np.divide(D,1.0+D)
        v = np.where(np.random.random(len(D)) <= D)[0]

    k = len(v)
    V = V[:,v.astype(int)]
    V = np.squeeze(V)
    if len(V.shape) == 1:
        V = V.reshape(len(V), 1)
    Y = np.zeros(k)

    for i in np.arange(k-1,-1,-1):
        P = np.sum(V*V, 1)
        P = P / np.sum(P)
        Y[i] = np.where(random.random() <= np.cumsum(P))[0][0]

        # correct
        j = np.nonzero(V[int(Y[i]),:])[0][0]
        Vj = V[:, j]
        V = V[:, [i for i in range(0,len(V[0])) if i != j]]
        # update V
        V = V - np.multiply.outer(Vj, V[int(Y[i]),:]/Vj[int(Y[i])])
        for a in range(0, i-1):
            for b in range(0, a-1):
                V[:,a] = V[:,a] - np.transpose(V[:,a])*V[:,b]*V[:,b]
            V[:,a] = V[:,a] / np.linalg.norm(V[:,a])
    return np.sort(Y)

def plot(x_reshape, y_reshape, dpp_sample, ind_sample):
    # plot
    plt.subplot(1,2,1)
    plt.plot(x_reshape[dpp_sample], y_reshape[dpp_sample], '.')
    for i in range(len(x_reshape)):
        if i in dpp_sample:
            plt.text(x_reshape[i], y_reshape[i], str(i), color='red', fontsize=0)
        else:
            pass
    plt.title('DPP')

    plt.subplot(1,2,2)
    plt.plot(x_reshape[ind_sample], y_reshape[ind_sample], '.')
    for i in range(len(x_reshape)):
        if i in ind_sample:
            plt.text(x_reshape[i], y_reshape[i], str(i), color='red', fontsize=0)
        else:
            pass
    plt.title('Random')
    plt.show()

def standard_experiment(L, k, x_reshape, y_reshape):
    # checkpoint; eigenvectors different from matlab, not sure why?

    dpp_sample = sample_dpp(L, k).astype(int)
    ind_sample = np.random.choice(int(len(L)),k)

    return dpp_sample

    # plot(x_reshape, y_reshape, dpp_sample, ind_sample)

def DPP(n,k):
    sigma = 1
    [x, y] = np.meshgrid(np.arange(1,float(n+1))/n, np.arange(1,float(n+1))/n)

    x_reshape = x.flatten()
    y_reshape = y.flatten()

    L = np.exp(-(np.square(np.subtract.outer(x_reshape,np.transpose(x_reshape))) + np.square(np.subtract.outer(y_reshape,np.transpose(y_reshape)))) / (sigma*sigma))
    # standard_experiment(L, k, x_reshape, y_reshape)
    dpp_sample = sample_dpp(L, k).astype(int)

    points=[]
    for i in dpp_sample:
        points.append(x_reshape[i])
        points.append(y_reshape[i])
    return points
