from ex2 import *

expected_cases = [
    ("GGTTGACTA", "TGTTACGG", 13, 'GTTGAC', 'GTT-AC'),
                  ("AA", "AA", 6, 'AA', 'AA'),
                  ("AA", "AAA", 6, 'AA', 'AA'),
                  ]
files = [
         "Foxp1_Bos_taurus",
         "Foxp1_Homo_sapiens",
         "Foxp1_Mus_musculus"
]
pairs = [
    ["Foxp1_Homo_sapiens", "Foxp1_Bos_taurus"],
    ["Foxp1_Mus_musculus", "Foxp1_Bos_taurus"],
    ["Foxp1_Homo_sapiens","Foxp1_Mus_musculus"]
]
find_strs_cases = [
    [("AAGAGAGTTAGAGTCAGC", "AG", 2), 2],
    [("AAGAGAGTTAGAGTCAGC", "AG", 3), 1],
    [("AAA", "A", 3), 1],
    [("AGAA", "AG", 1), 1],
    [("AGAAGAG", "AG", 2), 1],
    [("AGAAGAG", "AG", 3), 0],
]

find_strs3_cases = [
    [("AAAGGAGGTGTTCGGTCGTCGTC", 2), 4],
    [("AAAGGAGGTGTTCGGTCGTCGTC", 3), 1],
    [("AAAAAA", 3), 0],
    [("ATGATG", 2), 1],
    [("ATGATG", 3), 0],
]

def read_file(path_file):
    f = open(path_file, "r")
    return f.read()

def write_to_file(path,output):
    f = open(path, "a")
    f.write(str(output))
    f.close()

def check_first_function():
    try:
        for case in expected_cases:
            S, T, expected_score, expected_s, expected_t = case
            res_score, s_res, t_res = local_pwalignment(S, T, score=score)
            assert res_score == expected_score
            assert s_res == expected_s
            assert t_res == expected_t
    except AssertionError:
        print("a. local_pwalignment is NOT OK")
        raise
    print("a. local_pwalignment is  OK")

def compare_dicts(x, y):
    shared_items = {k: x[k] for k in x if k in y and x[k] == y[k]}
    return len(shared_items) == len(x.keys())


def compare_fastas():
    try:
        mydict = dict()
        for file in files:
            mydict[file] = read_fasta(f"{file}.fasta")
        for file1, file2 in pairs:
                score, s, t = local_pwalignment(mydict[file1], mydict[file2])
                name = f"{file1}_vs_{file2}"
                score_k, s_k, t_k = (read_file(f"Kim_files_d\{name}.fasta")).split(",")
                assert int(score_k) == int(score)
                assert s_k == s
                assert t_k == t
    except AssertionError:
        print(f"d. In comparison with kim's files there is a problem with {name}")
        raise
    print("d. In comparison with kim's files all is good")


def create_fastas():
    mydict = dict()
    res_dict = dict()
    for file in files:
        mydict[file] = read_fasta(f"{file}.fasta")
    for file1, file2 in pairs:
            res = local_pwalignment(mydict[file1], mydict[file2])
            res_score, res_s, res_t = res
            name = f"{file1}_vs_{file2}"
            print(f"the score for {name} is {res_score}")
            res_dict[name] = f">{file1}\n{res_s}\n>{file2}\n{res_t}"
    for key, value in res_dict.items():
        write_to_file(f"{key}.fasta",value)
    return res_dict

def largest_indel():
    mydict = dict()
    for file in files:
        mydict[file] = read_fasta(f"{file}.fasta")
    for file1, file2 in pairs:
            res = local_pwalignment(mydict[file1], mydict[file2])
            res_score, res_s, res_t = res
            indel = max(find_largest_indel(res_s),find_largest_indel(res_t))
            print(f"the indel length for {file1} and {file2} is {indel}")


def find_largest_indel(S):
    best = 0
    curr = 0
    for i in range(len(S)):
        if S[i] == '-':
            curr+=1
            continue
        else:
            best = max(best,curr)
            curr = 0
    return best



def create_fastas_for_comparison():
    mydict = dict()
    res_dict = dict()
    for file in files:
        mydict[file] = read_fasta(f"{file}.fasta")
    for file1, file2 in pairs:
            res = local_pwalignment(mydict[file1], mydict[file2])
            res_score, res_s, res_t = res
            name = f"{file1}_vs_{file2}"
            print(f"the score for {name} is {res_score}")
            res_dict[name] = f"{res_score},{res_s},{res_t}"
    for key, value in res_dict.items():
        write_to_file(f"{key}.fasta",value)
    return res_dict
# B
def find_strs_test():
    try:
        for case in find_strs_cases:
            S,s,r = case[0]
            expectes_result = case[1]
            assert find_strs(S,s,r) == expectes_result
    except AssertionError:
        print("g. find_strs_test is NOT OK")
        raise
    print("g. find_strs_test is OK")


def find_strs_test3():
    try:
        for case in find_strs3_cases:
            S, r = case[0]
            expectes_result = case[1]
            assert find_strs3(S,r) == expectes_result
    except AssertionError:
        print("h. find_strs_test3 is NOT OK")
        raise
    print("h. find_strs_test3 is OK")

def find_strs_test3_on_fasta():
    fasta = read_fasta('genome.fasta')
    result = find_strs3(fasta ,3)
    print(f"h. The result of running find_strs3 on genome.fasta is {result}")


def permutation_test_t():
    S =  read_fasta('genome.fasta')
    r = 3
    res = permutation_test(S,r)
    print(f"i. the answer is {res}")

if __name__ == '__main__':
    # A
    check_first_function()
    # create_fastas()
    # create_fastas_for_comparison()
    largest_indel()
    # compare_fastas()
    # # B
    find_strs_test()
    find_strs_test3()
    find_strs_test3_on_fasta()
    # permutation_test_t()
