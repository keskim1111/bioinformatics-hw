# Bioinfotools exercise 4
# Put your code instead of the 'pass' statements
import math
import random
import numpy as np
import scipy as sp
from scipy.linalg import lstsq
import random
from collections import defaultdict



# helpers
def find_k_nearest_indexes(X, u, k):
    distances = np.linalg.norm(X - u, axis=1)
    idx = distances.argsort()[:k]
    return idx.tolist()


def knn(X, y, u, K):
    # Expecting y to be with shape == (n,1)
    k_nearest_idx = find_k_nearest_indexes(X, u, K)
    return round(1 / K * y[k_nearest_idx, :].sum())


def randonly_divide_idx(X, L):
    n = X.shape[0]
    lst = np.array(range(n))
    np.random.shuffle(lst)
    groups = np.array_split(lst, L)
    return [group.tolist() for group in groups]


def L_fold_knn(L, K, X, y):
    # Expecting y to be with shape == (n,1)
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
            if res == y[idx]:
                cnt += 1
    return cnt / n


def calc_ssx(b,X,y,is_ssr=False, is_sst=False):
    curr_sum = 0
    for i in range(X.shape[1]):
        y_avg = y.mean()
        y_original = y[i]
        y_predicted = np.dot(X, b)
        if is_ssr:
            return sum(pow( y_predicted - y_avg , 2))
        elif is_sst:
            return sum(pow(y - y_avg , 2))
    return curr_sum


def lin_reg(X, y):
    # Expecting y to be with shape == (n,1)
    ones = np.ones(shape=(X.shape[0], 1), dtype=np.int64)
    X = np.c_[X, ones]
    # X = np.matrix(np.array(ones,X)).transpose()
    b = lstsq(X, y)[0]
    ssr = calc_ssx(b,X,y,is_ssr=True)
    sst = calc_ssx(b,X,y,is_sst=True)
    return b , ssr/sst

# you may change the function signature and add parameters (with or without default values)
def my_procedue(ge_path, age_path, bmi_path, gender_path, gene_names_path):
    X = load_txt_as_numpy(ge_path)
    age = load_txt_as_numpy(age_path).reshape((-1, 1))
    bmi = load_txt_as_numpy(bmi_path).reshape((-1, 1))
    gender = load_txt_as_numpy(gender_path).reshape((-1, 1))
    ABG = np.c_[age,bmi,gender]
    n = X.shape[0]
    gene_to_num_dict = defaultdict(int)
    line_num_to_gene_name_dict  = line_num_to_gene_name(gene_names_path)

    for i in range(100):
        # Preparing data
        indexes = random.sample(range(n), int(n / 2))
        new_X = X[indexes]
        new_ABG = ABG[indexes]
        # print(f"----Iteration number {i}-----\n")
        best_gene, max_gene_index, max_r_power_2 = one_iteration(new_X, new_ABG, line_num_to_gene_name_dict)
        # print(f"Winner gene with highst R^2 is {best_gene}\n {max_r_power_2}")
        gene_to_num_dict[best_gene]+=1
    return sorted(((v, k) for k, v in gene_to_num_dict.items()), reverse=True)[:1]

def load_txt_as_numpy(path):
    return np.loadtxt(path)

def one_iteration(X, ABG, line_num_to_gene_name_dict):
    max_r_power_2 = float("-inf")
    max_gene_index = 0
    for j in range(X.shape[1]): #for col in X - for gene
        X_j = X[:,j]
        b, r_power_2 = lin_reg(ABG ,X_j)
        if r_power_2 > max_r_power_2:
            max_r_power_2 = r_power_2
            max_gene_index = j
    best_gene = line_num_to_gene_name_dict[max_gene_index]
    return best_gene, max_gene_index, max_r_power_2

def line_num_to_gene_name(file_name):
    d = {}
    curr_line = 0
    with open(file_name, 'r') as read_obj:
        for line_read in read_obj:
            gene_read = "".join(line_read.split())
            d[curr_line] = gene_read
            curr_line += 1
    return d

if __name__ == '__main__':
    print(my_procedue("ge.txt","age.txt","bmi.txt","gender.txt","gene_names.txt"))
    pass
