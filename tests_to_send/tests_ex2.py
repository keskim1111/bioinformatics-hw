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
    ["Foxp1_Bos_taurus","Foxp1_Homo_sapiens"],
    ["Foxp1_Bos_taurus","Foxp1_Mus_musculus"],
    ["Foxp1_Mus_musculus", "Foxp1_Homo_sapiens"]
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
directory = "Kim_files_d"

success_msg = "is OK according to the limited tests in this file "
fail_msg = "is NOT OK according to the limited tests in this file "

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
        print(f"a. local_pwalignment {fail_msg}")
        raise
    print(f"a. local_pwalignment {success_msg}")

def compare_fastas():
    try:
        mydict = dict()
        for file in files:
            mydict[file] = read_fasta(f"{file}.fasta")
        for file1, file2 in pairs:
                score, s, t = local_pwalignment(mydict[file1], mydict[file2])
                name = f"{file1}_vs_{file2}"
                score_k, s_k, t_k = (read_file(f"{directory}\{name}.fasta")).split(",")
                assert int(score_k) == int(score)
                assert s_k == s
                assert t_k == t
    except AssertionError:
        print(f"d. In comparison with kim's files there is a problem with {name}")
        raise
    print(f"d. In comparison with kim's files the test {success_msg}")

# B
def find_strs_test():
    try:
        for case in find_strs_cases:
            S,s,r = case[0]
            expectes_result = case[1]
            assert find_strs(S,s,r) == expectes_result
    except AssertionError:
        print(f"g. find_strs_test{fail_msg}")
        raise
    print(f"g. find_strs_test {success_msg}")


def find_strs_test3():
    try:
        for case in find_strs3_cases:
            S, r = case[0]
            expectes_result = case[1]
            assert find_strs3(S,r) == expectes_result
    except AssertionError:
        print(f"h. find_strs_test3 {fail_msg}")
        raise
    print(f"h. find_strs_test3 {success_msg}")

def find_strs_test3_on_fasta():
    fasta = read_fasta('../genome.fasta')
    result = find_strs3(fasta ,3)
    print(f"h. The result of running find_strs3 on genome.fasta is {result}")



if __name__ == '__main__':
    # A
    check_first_function()
    compare_fastas()
    # # B
    find_strs_test()
    find_strs_test3()
    find_strs_test3_on_fasta()
