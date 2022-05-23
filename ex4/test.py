import pickle
from datetime import datetime
from pprint import pprint

import numpy as np
from ex4 import *
import os
import matplotlib.pyplot as plt
from collections import defaultdict
import itertools


def _pickle(fp, object="", is_load=False, is_dump=False):
    if is_dump:
        with open(fp, "wb") as f:
            pickle.dump(object, f)
    elif is_load:
        with open(fp, "rb") as f:
            return pickle.load(f)
    return None


def write_to_file(file, content, is_log=False):
    now = datetime.now()
    if is_log:
        current_time = now.strftime("%H:%M:%S")
        line = f"[{current_time}]: {content}\n"
        print(line)
    else:
        line = f"{content}"
    with open(file, "a") as f:
        f.write(line)
    return file


def test_a():
    X = np.array([
        [1, 1, 1],
        [4, 4, 4]])
    y = np.array([
        [1],
        [0],
    ])
    u = np.array([3, 3, 3])
    u2 = np.array([1, 0, 0])
    K = 1
    assert knn(X, y, u, K) == 0
    assert knn(X, y, u2, K) == 1


def test_b():
    X = np.array([
        [0, 0, 0],
        [1, 1, 1],
        [2, 2, 2],
        [3, 3, 3],
        [4, 4, 4],
        [5, 5, 5],

    ])
    y = np.array([
        [0],
        [0],
        [0],
        [1],
        [1],
        [1]
    ])
    print(y.shape)
    K = 2
    L = 3
    print(L_fold_knn(L, K, X, y))


def load_txt_as_numpy(path):
    return np.loadtxt(path)


def test_d(X,y):
    L = 10
    x_lst = []
    y_lst = []
    for k in range(1, 11):
        res = L_fold_knn(L, k, X, y)
        x_lst.append(k)
        y_lst.append(res)
    return x_lst, y_lst


def plot_fig_d(x, y, run_kind):
    plt.plot(x, y,label=run_kind)
    plt.title(f"Results d")
    plt.xlabel("K value")
    plt.ylabel("Accuracy")
    plt.legend(loc="upper left")
    plt.savefig(f"Results.png")


def runs_d():
    ### Using all the features in the data
    X = load_txt_as_numpy("ge.txt")
    y = load_txt_as_numpy("caco.txt")
    y = y.reshape((-1, 1))
    x_res, y_res = test_d(X,y)
    plot_fig_d(x_res, y_res, "only ge.txt")

    ### Using only the subset of genes that you found in section (c)
    paths = ["KEGG_WNT_SIGNALING_PATHWAY", "KEGG_LYSOSOME","KEGG_PPAR_SIGNALING_PATHWAY"]
    ## geting the relavent cols for sets chosen
    set_name_to_genes_list = create_set_name_to_genes_list_dict()
    gene_names_file_name_to_line_dict = gene_names_file_name_to_line_nums()
    name_to_idx = create_set_name_to_idx_gene_listdict(set_name_to_genes_list, gene_names_file_name_to_line_dict)
    relavent_ids_lists = [name_to_idx[path] for path in paths]
    relavent_cols_3_sets = []
    for listt in relavent_ids_lists:
        for item in listt:
            relavent_cols_3_sets.append(item)
    B = load_txt_as_numpy("ge.txt")
    y = load_txt_as_numpy("caco.txt")
    y = y.reshape((-1, 1))
    B = B[:,relavent_cols_3_sets]
    x_res, y_res = test_d(B,y)
    plot_fig_d(x_res, y_res, "3 pathways")

    ### Using only the subset of genes and age, gender and BMI as additional features
    C = load_txt_as_numpy("ge.txt")
    C = C[:,relavent_cols_3_sets]
    y = load_txt_as_numpy("caco.txt").reshape((-1, 1))
    age = load_txt_as_numpy("age.txt").reshape((-1, 1))
    bmi = load_txt_as_numpy("bmi.txt").reshape((-1, 1))
    gender = load_txt_as_numpy("gender.txt").reshape((-1, 1))
    C = np.c_[C, age,bmi,gender]
    x_res, y_res = test_d(C,y)
    plot_fig_d(x_res, y_res, "3 pathways +  3 files")
################################################################################################## C #############################################################################################
def create_set_name_to_genes_list_dict():
    f = open("kegg.MSigDB.symbols.gmt", "r")
    set_name_to_genes_list_dict = {}
    for line in f:
        line_splited = line.split("\t")
        name = "".join(line_splited[0].split())
        set_name_to_genes_list_dict[name] = line_splited[2:]
    return set_name_to_genes_list_dict
def gene_names_file_name_to_line_nums():
    file_name = "gene_names.txt"
    d = {}
    curr_line = 0
    with open(file_name, 'r') as read_obj:
        for line_read in read_obj:
            gene_read = "".join(line_read.split())
            d[gene_read] = curr_line
            curr_line += 1
    return d

def create_set_name_to_idx_gene_listdict(set_name_to_genes_list_dict, names_file_name_to_line_nums):
    set_name_to_idx_gene_list = defaultdict(list)
    set_genes_percent_dict = defaultdict(int)
    for gene_set in set_name_to_genes_list_dict.keys():
        for gene in set_name_to_genes_list_dict[gene_set]:
            if gene in names_file_name_to_line_nums:
                set_name_to_idx_gene_list[gene_set].append(names_file_name_to_line_nums[gene])
    return set_name_to_idx_gene_list
def check_if_string_in_file(file_name, string_to_search):
    """ Check if any line in the file contains given string """
    # Open the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            if string_to_search in line:
                return True
    return False


def search_for_3_best_sets(k):
    X = load_txt_as_numpy("ge.txt")
    y = load_txt_as_numpy("caco.txt")
    y = y.reshape((-1, 1))
    set_name_to_genes_list = create_set_name_to_genes_list_dict()
    gene_names_file_name_to_line_dict = gene_names_file_name_to_line_nums()
    name_to_idx = create_set_name_to_idx_gene_listdict(set_name_to_genes_list, gene_names_file_name_to_line_dict)
    set_name_to_accuracy = {}
    for set,idx in name_to_idx.items():
        # take only certain cols
        res = L_fold_knn(L=10,K=k,X=X[:,idx],y=y)
        set_name_to_accuracy[set] = res
    return set_name_to_accuracy

def test_e():
    X = np.array([
        [0],
        [1],
        [2],

    ])
    y = np.array([
        [0],
        [1],
        [2],
    ])

    X1 = np.array([
        [0, 0, 0],
        [1, 1, 1],
        [2, 2, 2],
        [3, 3, 3],
        [4, 4, 4],
        [5, 5, 5],

    ])
    y1 = np.array([
        [0],
        [0],
        [0],
        [1],
        [1],
        [1]
    ])

    print(lin_reg(X,y))

def line_num_to_gene_name():
    file_name = "gene_names.txt"
    d = {}
    curr_line = 0
    with open(file_name, 'r') as read_obj:
        for line_read in read_obj:
            gene_read = "".join(line_read.split())
            d[curr_line] = gene_read
            curr_line += 1
    return d

def test_f():
    X = load_txt_as_numpy("ge.txt")
    age = load_txt_as_numpy("age.txt").reshape((-1, 1))
    bmi = load_txt_as_numpy("bmi.txt").reshape((-1, 1))
    gender = load_txt_as_numpy("gender.txt").reshape((-1, 1))
    ABG = np.c_[age,bmi,gender]
    max_r_power_2 = float("-inf")
    max_gene_index = 0
    for j in range(X.shape[1]): #for col in X - for gene
        X_j = X[:,j]
        b, r_power_2 = lin_reg(ABG,X_j)
        if r_power_2 > max_r_power_2:
            max_r_power_2 = r_power_2
            max_gene_index = j
    d  = line_num_to_gene_name()
    best_gene = d[max_gene_index]
    return best_gene, max_gene_index, max_r_power_2


if __name__ == '__main__':
    # you can write whatever you want here
    # test_a()
    # test_b()
    # runs_d()
    # for k in range(1,10):
    #     d=  search_for_3_best_sets(k)
    #     print(sorted(((v, k) for k, v in d.items()), reverse=True)[:4])
    # pprint(_pickle("genes_missing_set.dat", is_load=True))
    print(test_f())
