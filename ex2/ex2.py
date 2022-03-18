# Bioinfotools exercise 2
# Put your code instead of the 'pass' statements
import random
import re
import time

letters = ['A', 'C', 'G', 'T']
options = ["".join([x, y, z]) for x in letters for y in letters for z in letters]
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


#
def find_aligned(optloc, B, S, T):
    i,j = optloc
    print(f"optloc is {i},{j}")
    curr_pos = optloc
    align_S = ""
    align_T = ""
    i, j = curr_pos
    align_S = S[i] + align_S
    align_T = T[j] + align_T
    while curr_pos != (None, None):
        print(align_S, align_T)
        next_pos = B[curr_pos[0]][curr_pos[1]]
        if next_pos ==(None, None):
            break
        i,j = curr_pos
        i_next, j_next = next_pos
        if i_next - i == 1 and j_next - j == 1:
            align_S = S[i_next ] + align_S
            align_T = T[j_next ] + align_T
        elif i_next - i == 1:
            align_S = S[i_next ] + align_S
            align_T = "-" + align_T
        else:# j_old - j == 1
            align_S = "-" + align_S
            align_T = T[j_next ] + align_T
    align_S = S[i_next ] + align_S
    align_T = T[j_next ] + align_T
    return align_S, align_T


def traceback(S, T, index, B):
    allign_s = ""
    allign_t = ""
    cur = index
    while (cur != (None,None)):
        i, j = cur
        in_cur = B[cur[0]][cur[1]]
        diagonal = (i-1, j-1)
        up = (i-1, j)
        left = (i, j-1)
        if (in_cur == diagonal):
            allign_s = S[i - 1] + allign_s
            allign_t = T[j - 1] + allign_t
        elif (in_cur == up):
            allign_s = S[i - 1] + allign_s
            allign_t = "-" + allign_t
        elif (in_cur == left ):
            allign_s = "-" + allign_s
            allign_t = T[j - 1] + allign_t
        cur = in_cur
    return (allign_s, allign_t)

def fill_matrix(S, T, score=score):
    # create matrix
    A = []
    B = []
    for i in range(len(S)+1):
        A.append([0] * (len(T) + 1))
        B.append( [(None,None) for i in range(len(T) + 1)])

    best = 0
    optloc = (0, 0)
    for i in range(1, len(S)+1):
        for j in range(1, len(T)+1):
            my_dict = {
                (i, j-1): A[i][j - 1] + score("-", T[j-1]),
                (i-1, j): A[i - 1][j] + score(S[i-1], "-"),
                (i-1, j-1): A[i - 1][j - 1] + score(S[i-1], T[j-1]),
                       }
            A[i][j] = max(
                0, max((my_dict.values()))
            )

            if A[i][j] >= best:
                best = A[i][j]
                optloc = (i, j)
            if A[i][j] == 0:
                continue
            B[i][j] = max(my_dict, key=my_dict.get)
    return best, optloc, A, B


def local_pwalignment(S, T, score=score):
    best, optloc, A, B = fill_matrix(S, T)
    s, t = traceback(S, T, optloc, B)
    return best, s, t


# Part B
def find_strs(S, s, r):
    reg_rule = f"({s}){{{r},}}"
    res = re.findall(reg_rule, S)
    return len(res)


def find_strs3(S, r):
    count = 0
    for option in options:
        if option not in ['AAA', 'TTT', 'GGG', 'CCC']:
            count += find_strs(S, option, r)
    return count

# TODO remove prints
def permutation_test(S, r):
    start1 = time.time()
    count = 0
    samples = 100
    for i in range(100):
        print(f"starting the {i}th test...")
        start = time.time()
        random_S = ''.join(random.sample(S, len(S)))
        if find_strs3(S, r) < find_strs3(random_S, r):
            count+=1
            print(f"the answer for the {i}th test is True ")
        else:
            print(f"the answer for the {i}th test is False ")
        end = time.time()
        print(f"the {i}th test took {end-start} seconds")
        print(f"{100-i} tests left\n\n")
    end1 = time.time()
    print(count)
    print(f"the whole program took {end1-start1}")

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