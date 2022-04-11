# Bioinfotools exercise 3
# Put your code instead of the 'pass' statements
import math
import random
from Bio.Align import substitution_matrices
from Bio.SubsMat import MatrixInfo
from Bio import pairwise2, SeqIO
import copy


# Helpers for upgma
def find_smallest_distance_cell(distance_matrix):
    try:
        smallest_val = float("inf")
        smallest_row = -2
        smallest_col = -2
        for i in range(distance_matrix.rows):
            for j in range(i):
                val = distance_matrix.get(i, j)
                if val < smallest_val:
                    smallest_val = val
                    smallest_row, smallest_col = i, j
        return smallest_row, smallest_col
    except IndexError:
        print(f"indexes:{i},{j}")
        raise IndexError


def update_row_until_first(distance_matrix, first, second, names):
    # first item in each row is 0
    res = [0]
    for i in range(0, first):
        res.append(
            (names[first][1] * distance_matrix.get(first, i) + names[second][1] * distance_matrix.get(second, i)) / (
                    names[first][1] + names[second][1]))
    zeros = distance_matrix.rows - len(res)
    for j in range(zeros):
        res.append(0)
    distance_matrix.set_row(first, res)


def update_col(distance_matrix, first, second, names):
    for i in range(first + 1, second):
        distance_matrix.set(
            (names[first][1] * distance_matrix.get(i, first) + names[second][1] * distance_matrix.get(second, i)) / (
                    names[first][1] + names[second][1]), i, first)


def update_row_from_second(distance_matrix, first, second, names):
    for i in range(second + 1, distance_matrix.rows):
        distance_matrix.set(
            (distance_matrix.get(i, first) * names[first][1] + names[second][1] * distance_matrix.get(i, second)) / (
                    names[first][1] + names[second][1]), i, first)
        distance_matrix.delete_index(i, second)
    # removes a col
    distance_matrix.cols -= 1


def update_distance_matrix(distance_matrix, i, j, names):
    if j < i:
        i, j = j, i
    join_i_j_in_d(distance_matrix, i, j, names)


def join_i_j_in_d(distance_matrix, first, second, names):
    update_row_until_first(distance_matrix, first, second, names)
    update_col(distance_matrix, first, second, names)
    update_row_from_second(distance_matrix, first, second, names)
    distance_matrix.delete_row(second)
    # removed a row
    distance_matrix.rows -= 1


def update_names(names, i, j):
    if j < i:
        i, j = j, i
    names[i] = ((names[i][0], names[j][0]), names[i][1] + names[j][1])
    del names[j]


# Helpers for globalpw_dist”
def calc_score(S, T, sub_matrix):
    alignment_score = pairwise2.align.globalds(S, T, sub_matrix, -5, -5)[0].score
    return alignment_score


def create_d_matrix(similarity_matrix, n):
    max_val = similarity_matrix.find_maximal_value()
    distance_matrix = matrix(n, n)
    for i in range(similarity_matrix.rows):
        for j in range(similarity_matrix.cols):
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
            distance = calculate_k_mer_dist(seq_lst[i], seq_lst[j], seq_to_kmers)
            d.set(distance, i, j)
    return d


def calculate_k_mer_dist(seq1, seq2, seq_to_kmers_dict):
    res = 0
    k_mers_union = seq_to_kmers_dict[seq1].union(seq_to_kmers_dict[seq2])
    for k_mer in k_mers_union:
        c_res = num_of_k_mer_in_s(seq1, k_mer) - num_of_k_mer_in_s(seq2, k_mer)
        res += pow(c_res, 2)
    return math.sqrt(res)


# Helpers for eval_dist

def tree_splits_names(tree, lst):
    lst1 = []
    lst2 = []
    if isinstance(tree[0], tuple):
        lst1 = tree_splits_names(tree[0], lst)
    if isinstance(tree[1], tuple):
        lst2 = tree_splits_names(tree[1], lst)
    return lst1 + lst2 + [tree]

def read_fasta(path):
    seq_lst = []
    name_lst = []
    record_iterator = SeqIO.parse(path, "fasta")
    # infinite loop
    while True:
        try:
            # get the next item
            element = next(record_iterator)
            seq_lst.append(str(element.seq))
            name_lst.append(str(element.name))
        except StopIteration:
            # if StopIteration is raised, break from loop
            break
    return seq_lst, name_lst


class matrix:
    def __init__(self, n_rows, n_cols):
        assert n_cols > 0 and n_cols > 0
        self.matrix = [[0] * (n_cols + 1) for i in range(n_rows + 1)]
        self.rows = n_rows
        self.cols = n_cols

    def get(self, i, j):
        return self.matrix[i + 1][j + 1]

    def set(self, x, i, j):
        self.matrix[i + 1][j + 1] = x

    def delete_index(self, i, j):
        del self.matrix[i + 1][j + 1]

    def delete_row(self, row):
        del self.matrix[row + 1]

    def set_row(self, row_index, new_row):
        self.matrix[row_index + 1] = new_row

    def get_readable_matrix_string(self, matrix):
        strings = []
        for row in matrix:
            strings.append(str(row))
        return '\n'.join(strings)

    def __str__(self):
        return self.get_readable_matrix_string(self.matrix)

    def set_whole_matrix(self, input_matrix):
        new_rows_len = len(input_matrix)
        new_cols_len = len(input_matrix[0])
        if new_rows_len == self.rows and new_cols_len == self.cols:
            for i in range(self.rows):
                for j in range(self.cols):
                    val = input_matrix[i][j]
                    self.set(val, i, j)

    def find_maximal_value(self):
        maximal_value = float("-inf")
        for i in range(len(self.matrix) - 1):
            for j in range(len(self.matrix[i]) - 1):
                maximal_value = max(maximal_value, self.get(i, j))
        return maximal_value


# TODO find out D how it looks
def upgma(D, seq_names_lst):
    distance_matrix = copy.deepcopy(D)
    names = [(name, 1) for name in seq_names_lst]
    while len(names) > 1:
        i, j = find_smallest_distance_cell(distance_matrix)
        update_distance_matrix(distance_matrix, i, j, names)
        update_names(names, i, j)
    return names[0][0]


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
    distanse_matrix = dist_func(seq_lst)
    seq_lst_msa, name_lst = read_fasta(msa_aln_path)
    result_tree = upgma(distanse_matrix, name_lst)
    print(result_tree)
    print(type(result_tree))
    keys = tree_splits_names(result_tree, [])
    print(keys)

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
