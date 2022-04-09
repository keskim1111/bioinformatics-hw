# Bioinfotools exercise 3
# Put your code instead of the 'pass' statements
import math
import random
from Bio.Align import substitution_matrices
from Bio.SubsMat import MatrixInfo
from Bio import pairwise2


# Helpers for upgma
def find_smallest_distance_cell(D):
    smallest_val = float("inf")
    smallest_row = -2
    smallest_col = -2
    for i in range(len(D)):
        for j in range(len(D[i])):
            val = D[i][j]
            if val < smallest_val:
                smallest_val = val
                smallest_row, smallest_col = i, j
    return smallest_row, smallest_col


def update_row_until_first(d, first, second):
    res = []
    for i in range(0, first):
        res.append(d[first][i] / 2 + d[second][i] / 2)
    return res


def update_col(d, first, second):
    for i in range(first + 1, second):
        d[i][first] = (d[i][first] + d[second][i]) / 2


def update_row_from_second(d, first, second):
    for i in range(second + 1, len(d)):
        d[i][first] = (d[i][first] + d[i][second]) / 2
        del d[i][second]


def update_d(d, i, j):
    if i < j:
        join_i_j_in_d(d, i, j)
    else:
        join_i_j_in_d(d, j, i)


def join_i_j_in_d(d, first, second):
    d[first] = update_row_until_first(d, first, second)
    update_col(d, first, second)
    update_row_from_second(d, first, second)
    del d[second]


def update_names(names, i, j):
    if i < j:
        names[i] = f"({names[i]},{names[j]})"
        del names[j]
    else:  # j < i
        names[j] = f"({names[i]},{names[j]})"
        del names[i]


# Helpers for globalpw_dist”
def calc_score(S, T, sub_matrix):
    alignment_score = pairwise2.align.globalds(S, T, sub_matrix, -5, -5)[0].score
    return alignment_score


def create_d_matrix(similarity_matrix, n):
    max_val = similarity_matrix.find_maximal_value()
    distance_matrix = matrix(n, n)
    for i in range(len(similarity_matrix.matrix) - 1):
        for j in range(len(similarity_matrix.matrix[i]) - 1):
            similarity_value = similarity_matrix.get(i, j)
            distance_matrix.set(max_val - similarity_value + 1, i, j)
    return distance_matrix


def create_similarity_matrix(seq_lst, sub_matrix, n):
    similarity = matrix(n, n)
    for i in range(n):
        for j in range(n):
            S = seq_lst[i]
            T = seq_lst[j]
            similarity.set(calc_score(S, T, sub_matrix), i, j)
    return similarity


# Helpers for k_mer

def find_k_mers(S, k):
    res = set()
    optional_kmers_indexes = len(S) - k + 1
    for i in range(optional_kmers_indexes):
        res.add(S[i:i + k])
    return res


def num_of_k_mer_in_s(S, k_mer):
    count = 0
    k = len(k_mer)
    optional_kmers_indexes = len(S) - k + 1
    for i in range(optional_kmers_indexes):
        w = S[i:i + k]
        if w == k_mer:
            count += 1
    return count


def seq_lst_to_kmers_dict(seq_lst, k):
    d = dict()
    for seq in seq_lst:
        d[seq] = find_k_mers(seq, k)
    return d


def create_k_mer_dist_matrix(seq_lst, seq_to_kmers):
    n = len(seq_lst)
    d = matrix(n, n)
    for i in range(n):
        for j in range(n):
            distance = calcule_k_mer_dist(seq_lst[i], seq_lst[j], seq_to_kmers)
            d.set(distance, i, j)
    return d


def calcule_k_mer_dist(seq1, seq2, seq_to_kmers_dict):
    res = 0
    k_mers_union = seq_to_kmers_dict[seq1].union(seq_to_kmers_dict[seq2])
    for k_mer in k_mers_union:
        c_res = num_of_k_mer_in_s(seq1, k_mer) - num_of_k_mer_in_s(seq2, k_mer)
        res += pow(c_res, 2)
    return math.sqrt(res)


class matrix:
    def __init__(self, n_rows, n_cols):
        assert n_cols > 0 and n_cols > 0
        self.matrix = [[0] * (n_cols + 1) for i in range(n_rows + 1)]

    def get(self, i, j):
        return self.matrix[i + 1][j + 1]

    def set(self, x, i, j):
        self.matrix[i + 1][j + 1] = x

    def get_readable_matrix_string(self, matrix):
        strings = []
        for row in matrix:
            strings.append(str(row))
        return '\n'.join(strings)

    def __str__(self):
        return self.get_readable_matrix_string(self.matrix)

    def find_maximal_value(self):
        maximal_value = float("-inf")
        for i in range(len(self.matrix) - 1):
            for j in range(len(self.matrix[i]) - 1):
                maximal_value = max(maximal_value, self.get(i, j))
        return maximal_value


# TODO find out D how it looks
def upgma(D, seq_names_lst):
    d = D
    names = seq_names_lst
    while len(names) > 1:
        i, j = find_smallest_distance_cell(D)
        update_d(d, i, j)
        update_names(names, i, j)
    return names[0]


def globalpw_dist(seq_lst):
    n = len(seq_lst)
    sub_matrix = MatrixInfo.blosum62
    similarity_matrix = create_similarity_matrix(seq_lst, sub_matrix, n)
    distance_matrix = create_d_matrix(similarity_matrix, n)
    return distance_matrix


def kmer_dist(seq_lst, k=3):
    seq_to_kmers = seq_lst_to_kmers_dict(seq_lst, k)
    k_mer_distance_matrix = create_k_mer_dist_matrix(seq_lst, seq_to_kmers)
    return k_mer_distance_matrix


def eval_dist(seq_lst, msa_aln_path, dist_func=globalpw_dist):
    """

    :param seq_lst: list
        list of n sequences S1, S2, ..., Sn
    :param msa_aln_path: str
        ClustalW FASTA alignment file
    :param dist_func: a distance function name
    :return: dict
        keys are tuples representing each of the splits in the UPGMA tree T
        values are the counters of the splits in the random permutations
    """
    pass


if __name__ == '__main__':
    MAFFT_EXE_PATH = r"WHERE_YOU_UNZIPPED_MAFFT\usr\bin\mafft"  # depends on your operating system
    seqs_path = r"sequences.fasta"
    msa_aln_path = r"sequences.aln.fasta"
    # you can write whatever you want here
    pass
