from ex3 import upgma, globalpw_dist, find_k_mers, num_of_k_mer_in_s, kmer_dist


def alpha_labels(start, end):
    labels = []
    for i in range(ord(start), ord(end) + 1):
        labels.append(chr(i))
    return labels


# Test table data and corresponding labels
M_labels = alpha_labels("A", "G")  # A through G
M = [
    [0],
    [19],
    [27, 31],
    [8, 18, 26],
    [33, 36, 41, 31],
    [18, 1, 32, 17, 35],
    [13, 13, 29, 14, 28, 12]
]


# should output: '((((A,D),((B,F),G)),C),E)'
# print(upgma(M, M_labels))
# seq_lst = ['ACCGT', "ACG"]
# print(globalpw_dist(seq_lst))


def compare_sets(set1, set2):
    if len(set1) == len(set2):
        if set1 == set2:
            return True
        else:
            return False
    else:
        return False


def test_kmer():
    inputs = [("GATGAT", 3), ("GATGAT", 2), ("GATGAT", 5)]
    results = [{'ATG', 'GAT', 'TGA'}, {'AT', 'GA', 'TG'}, {'ATGAT', 'GATGA'}]
    for i in range(len(results)):
        word, k = inputs[i]
        assert find_k_mers(word, k) == results[i]


def test_num_of_k_mer_in_s():
    inputs = [
        ("GATGAT", 'ATG'),
        ("GATGAT", 'GAT'),
        ("GATGAT", 'TGA'),
        ("GATGAT", 'AT'),
        ("GATGAT", 'GATGA'),
        ("CACAC", "CAC")
    ]
    results = [
        1, 2, 1, 2, 1, 2

    ]
    for i in range(len(results)):
        word, k_mer = inputs[i]
        assert num_of_k_mer_in_s(word, k_mer) == results[i]


# TODO write test
def kmer_dist_test():
    lst = ["GATGAT", "CACAC"]
    print(kmer_dist(lst))


def closed_tests():
    test_kmer()
    test_num_of_k_mer_in_s()


if __name__ == '__main__':
    closed_tests()
    kmer_dist_test()
    pass
