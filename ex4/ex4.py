# Bioinfotools exercise 4
# Put your code instead of the 'pass' statements
import math
import random
import numpy as np
import scipy as sp
import random


# helpers
def find_k_nearest_indexes(X, u, k):
    distances = np.linalg.norm(X - u, axis=1)
    idx = distances.argsort()[:k]
    return idx.tolist()


def knn(X, y, u, K):
    k_nearest_idx = find_k_nearest_indexes(X, u, K)
    return round(1 / K * y[k_nearest_idx, :].sum())


def randonly_divide_idx(X, L):
    n = X.shape[0]
    lst = np.array(range(n))
    np.random.shuffle(lst)
    groups = np.array_split(lst, L)
    return [group.tolist() for group in groups]


def L_fold_knn(L, K, X, y):
    # print("X: ", X)
    n = X.shape[0]
    groups_idx = randonly_divide_idx(X, L)
    cnt = 0
    for i in range(len(groups_idx)):
        current_group_idx = groups_idx[i]
        other_groups = groups_idx[:i] + groups_idx[i + 1:]
        other_groups_idx = [item for sublist in other_groups for item in sublist]
        new_X = X[other_groups_idx, :]
        new_y = y[other_groups_idx, :]
        # check if knn result is like the real result in y for each item in current group
        for idx in current_group_idx:
            res = knn(new_X, new_y, X[idx], K)
            # TODO fix this
            if res == y[idx]:
                cnt += 1
    return cnt / n


def lin_reg(X, y):
    pass


# you may change the function signature and add parameters (with or without default values)
def my_procedue():
    pass


if __name__ == '__main__':
    # you can write whatever you want here
    pass
