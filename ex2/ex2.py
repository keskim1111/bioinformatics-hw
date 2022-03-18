# Bioinfotools exercise 2
# Put your code instead of the 'pass' statements
import random
import re
import time

# Auxiliary function to read a FASTA file. Returns a string

def read_fasta(filename):
    f = open(filename)
    header = f.readline()
    x = f.read()
    x = x.replace("\n", "")
    return x


# Part A
# Scoring regime
# TODO return to before
def score(a, b):
    if (a == b):
        return 3
    elif ((a == "-" and b != "-") or (a != "-" and b == "-")):
        return -2
    else:
        return -3


def traceback(S, T, index, B):
    allign_s = ""
    allign_t = ""
    cur = index
    while cur != (None, None):
        i, j = cur
        in_cur = B[cur[0]][cur[1]]
        diagonal = (i - 1, j - 1)
        up = (i - 1, j)
        left = (i, j - 1)
        if (in_cur == diagonal):
            allign_s = S[i - 1] + allign_s
            allign_t = T[j - 1] + allign_t
        elif (in_cur == up):
            allign_s = S[i - 1] + allign_s
            allign_t = "-" + allign_t
        elif (in_cur == left):
            allign_s = "-" + allign_s
            allign_t = T[j - 1] + allign_t
        cur = in_cur
    return allign_s, allign_t


def fill_matrix(S, T, score):
    # create matrix
    A = [] # will save the scores
    B = [] # will save the path
    for i in range(len(S) + 1):
        A.append([0] * (len(T) + 1))
        B.append([(None, None) for i in range(len(T) + 1)])
    best = 0
    optimal_location = (0, 0)
    for i in range(1, len(S) + 1):
        for j in range(1, len(T) + 1):
            my_dict = {
                (i, j - 1): A[i][j - 1] + score("-", T[j - 1]),
                (i - 1, j): A[i - 1][j] + score(S[i - 1], "-"),
                (i - 1, j - 1): A[i - 1][j - 1] + score(S[i - 1], T[j - 1]),
            }
            A[i][j] = max(
                0, max((my_dict.values()))
            )
            if A[i][j] >= best:
                best = A[i][j]
                optimal_location = (i, j)
            if A[i][j] == 0:
                continue
            B[i][j] = max(my_dict, key=my_dict.get)
    return best, optimal_location, A, B


def local_pwalignment(S, T, score=score):
    best, optloc, A, B = fill_matrix(S, T, score)
    s, t = traceback(S, T, optloc, B)
    return best, s, t


# Part B
def find_strs(S, s, r):
    reg_rule = f"({s}){{{r},}}"
    res = re.findall(reg_rule, S)
    return len(res)


def find_strs3(S, r):
    letters = ['A', 'C', 'G', 'T']
    options = ["".join([x, y, z]) for x in letters for y in letters for z in letters]
    count = 0
    for option in options:
        if option not in ['AAA', 'TTT', 'GGG', 'CCC']:
            count += find_strs(S, option, r)
    return count


def permutation_test(S, r):
    count = 0
    samples = 100
    for i in range(100):
        random_s = ''.join(random.sample(S, len(S)))
        if find_strs3(S, r) < find_strs3(random_s, r):
            count += 1

    return count / samples < 0.05


if __name__ == '__main__':
    ## Part a
    foxp1_human_path = 'Foxp1_Homo_sapiens.fasta'
    foxp1_mus_path = 'Foxp1_Mus_musculus.fasta'
    foxp1_bos_path = 'Foxp1_Bos_taurus.fasta'
    foxp1_human = read_fasta(foxp1_human_path)
    foxp1_mus = read_fasta(foxp1_mus_path)
    foxp1_bos = read_fasta(foxp1_bos_path)
    print(local_pwalignment("GGTTGACTA", "TGTTACGG", score=score)[0] == 13)
    # Part b
    print(find_strs("AAGAGAGTTAGAGTCAGC", "AG", 2) == 2)
    print(find_strs("AAGAGAGTTAGAGTCAGC", "AG", 3) == 1)
    print(find_strs3("AAAGGAGGTGTTCGGTCGTCGTC", 2) == 4)
    print(find_strs3("AAAGGAGGTGTTCGGTCGTCGTC", 3) == 1)
    genome_path = "genome.fasta"
    pass
