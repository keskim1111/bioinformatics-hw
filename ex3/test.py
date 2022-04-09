from ex3 import upgma, globalpw_dist, find_k_mers, num_of_k_mer_in_s, kmer_dist, matrix


def alpha_labels(start, end):
    labels = []
    for i in range(ord(start), ord(end) + 1):
        labels.append(chr(i))
    return labels


def test_upgma():
    # Test table data and corresponding labels
    M_labels = alpha_labels("A", "G")  # A through G
    dist_matrix = matrix(7, 6)
    M = [
        [0, 0, 0, 0, 0, 0],
        [19, 0, 0, 0, 0, 0],
        [27, 31, 0, 0, 0, 0],
        [8, 18, 26, 0, 0, 0],
        [33, 36, 41, 31, 0, 0],
        [18, 1, 32, 17, 35, 0],
        [13, 13, 29, 14, 28, 12]
    ]
    dist_matrix.set_whole_matrix(M)

    # should output: '((((A,D),((B,F),G)),C),E)'
    print(upgma(dist_matrix, M_labels))


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
    kmer_dist_test()


if __name__ == '__main__':
    # closed_tests()
    test_upgma()
    pass
