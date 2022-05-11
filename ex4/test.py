import pickle
from datetime import datetime
from pprint import pprint

import numpy as np
from ex4 import *
import os
import matplotlib.pyplot as plt
from collections import defaultdict


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
    X = np.array([[1, 1, 1], [4, 4, 4]])
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


def test_d(X):
    y = load_txt_as_numpy("caco.txt")
    y = y.reshape((-1, 1))
    L = 10
    x_lst = []
    y_lst = []
    for k in range(1, 11):
        res = L_fold_knn(L, k, X, y)
        x_lst.append(k)
        y_lst.append(res)
    return x_lst, y_lst


def plot_fig_d(x, y, run_kind):
    plt.scatter(x, y)
    plt.title(f"Accuracy per k value {run_kind}")
    plt.xlabel("K value")
    plt.ylabel("Accuracy")
    plt.savefig(f"{run_kind}.png")


def runs_d():
    ### Using all the features in the data
    X = load_txt_as_numpy("ge.txt")
    x, y = test_d(X)
    plot_fig_d(x, y, "Using all the features in the data")
    ### Using only the subset of genes that you found in section (c)
    ### Using only the subset of genes and age, gender and BMI as additional features


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
    genes_from_keg = set()
    genes_missing = set()
    for gene_set in set_name_to_genes_list_dict.keys():
        for gene in set_name_to_genes_list_dict[gene_set]:
            genes_from_keg.add(gene)
            if gene in names_file_name_to_line_nums:
                set_name_to_idx_gene_list[gene_set].append(names_file_name_to_line_nums[gene])
            else:
                genes_missing.add(gene)
    print("genes_from_keg: \n", genes_from_keg)
    print("genes_missing_set: \n", genes_missing)
    print(f"per of genes not apearing is {len(genes_missing) / len(genes_from_keg)}")
    _pickle("genes_missing_set.dat", genes_missing, is_dump=True)
    write_to_file("genes_missing_set.txt", genes_missing)
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


def search_for_3_best_sets():
    set_name_to_genes_list = create_set_name_to_genes_list_dict()
    gene_names_file_name_to_line_dict = gene_names_file_name_to_line_nums()
    name_to_idx = create_set_name_to_idx_gene_listdict(set_name_to_genes_list, gene_names_file_name_to_line_dict)


if __name__ == '__main__':
    # you can write whatever you want here
    # test_a()
    # test_b()
    runs_d()
    # search_for_3_best_sets()
    # pprint(_pickle("genes_missing_set.dat", is_load=True))
