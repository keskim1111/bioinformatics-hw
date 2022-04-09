import functools
import time
from datetime import datetime


from ex3 import upgma, globalpw_dist, find_k_mers, num_of_k_mer_in_s, kmer_dist, matrix
from Bio import SeqIO

def timeit(func):
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        elapsed_time = time.time() - start_time
        print('function [{}] finished in {} ms'.format(
            func.__name__, int(elapsed_time * 1_000)))
        return result
    return new_func

dist_kmer = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 0.0, 33.15116890850155, 19.28730152198591, 43.64630568559039, 39.01281840626232, 35.04283093587046, 45.8257569495584, 44.31703961232068, 42.5440947723653, 38.98717737923585, 44.28317965096906, 44.48595283907045],
[0, 33.15116890850155, 0.0, 33.24154027718932, 44.69899327725402, 41.19465984809196, 28.442925306655784, 47.275786614291256, 46.357307945997036, 43.89760813529594, 40.779897008207364, 45.60701700396552, 46.32493928760188],
[0, 19.28730152198591, 33.24154027718932, 0.0, 44.07947368106838, 39.57271787481876, 35.58089374931439, 46.9041575982343, 44.54211490264017, 42.731721238442994, 39.44616584663204, 44.50842616853577, 44.91102314577124],
[0, 43.64630568559039, 44.69899327725402, 44.07947368106838, 0.0, 44.170125650715555, 45.81484475582123, 49.48737212663449, 47.02127178203499, 44.21538193886829, 43.0, 46.58325879540846, 44.11349000022555],
[0, 39.01281840626232, 41.19465984809196, 39.57271787481876, 44.170125650715555, 0.0, 42.28474902373195, 48.062459362791664, 44.11349000022555, 42.35563716909474, 39.03844259188627, 43.78355855797927, 44.21538193886829],
[0, 35.04283093587046, 28.442925306655784, 35.58089374931439, 45.81484475582123, 42.28474902373195, 0.0, 48.641546028061235, 47.60252094164762, 45.36518488885502, 42.14261501141095, 46.67976006793523, 47.57099956906519],
[0, 45.8257569495584, 47.275786614291256, 46.9041575982343, 49.48737212663449, 48.062459362791664, 48.641546028061235, 0.0, 50.23942674832188, 48.55924216871593, 46.75467891024384, 50.82322303829225, 50.46781152378217],
[0, 44.31703961232068, 46.357307945997036, 44.54211490264017, 47.02127178203499, 44.11349000022555, 47.60252094164762, 50.23942674832188, 0.0, 46.4327470649756, 43.3358973600409, 41.291645644125154, 47.738873049120045],
[0, 42.5440947723653, 43.89760813529594, 42.731721238442994, 44.21538193886829, 42.35563716909474, 45.36518488885502, 48.55924216871593, 46.4327470649756, 0.0, 42.40283009422838, 45.967379738244816, 40.95119045888654],
[0, 38.98717737923585, 40.779897008207364, 39.44616584663204, 43.0, 39.03844259188627, 42.14261501141095, 46.75467891024384, 43.3358973600409, 42.40283009422838, 0.0, 41.773197148410844, 44.30575583375144],
[0, 44.28317965096906, 45.60701700396552, 44.50842616853577, 46.58325879540846, 43.78355855797927, 46.67976006793523, 50.82322303829225, 41.291645644125154, 45.967379738244816, 41.773197148410844, 0.0, 47.20169488482379],
[0, 44.48595283907045, 46.32493928760188, 44.91102314577124, 44.11349000022555, 44.21538193886829, 47.57099956906519, 50.46781152378217, 47.738873049120045, 40.95119045888654, 44.30575583375144, 47.20169488482379, 0.0]]

dist_global =[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[0, 1609.0, 3104.0, 1959.0, 6048.0, 4600.0, 3500.0, 5461.0, 6027.0, 5320.0, 4681.0, 5799.0, 5574.0],
[0, 3104.0, 1370.0, 3016.0, 5914.0, 4795.0, 2436.0, 5050.0, 5800.0, 5060.0, 4950.0, 5554.0, 5307.0],
[0, 1959.0, 3016.0, 1542.0, 6002.0, 4630.0, 3410.0, 5344.0, 5980.0, 5239.0, 4752.0, 5720.0, 5451.0],
[0, 6048.0, 5914.0, 6002.0, 1158.0, 6168.0, 5970.0, 6409.0, 6182.0, 5996.0, 6257.0, 6105.0, 5953.0],
[0, 4600.0, 4795.0, 4630.0, 6168.0, 1671.0, 5073.0, 6004.0, 5928.0, 5539.0, 4889.0, 5832.0, 5719.0],
[0, 3500.0, 2436.0, 3410.0, 5970.0, 5073.0, 932.0, 4795.0, 5483.0, 5050.0, 5349.0, 5323.0, 5160.0],
[0, 5461.0, 5050.0, 5344.0, 6409.0, 6004.0, 4795.0, 1.0, 5389.0, 5637.0, 6258.0, 5433.0, 5585.0],
[0, 6027.0, 5800.0, 5980.0, 6182.0, 5928.0, 5483.0, 5389.0, 505.0, 5778.0, 5800.0, 3400.0, 5641.0],
[0, 5320.0, 5060.0, 5239.0, 5996.0, 5539.0, 5050.0, 5637.0, 5778.0, 991.0, 5672.0, 5586.0, 3738.0],
[0, 4681.0, 4950.0, 4752.0, 6257.0, 4889.0, 5349.0, 6258.0, 5800.0, 5672.0, 1800.0, 5511.0, 5912.0],
[0, 5799.0, 5554.0, 5720.0, 6105.0, 5832.0, 5323.0, 5433.0, 3400.0, 5586.0, 5511.0, 755.0, 5536.0],
[0, 5574.0, 5307.0, 5451.0, 5953.0, 5719.0, 5160.0, 5585.0, 5641.0, 3738.0, 5912.0, 5536.0, 834.0]]
def alpha_labels(start, end):
    labels = []
    for i in range(ord(start), ord(end) + 1):
        labels.append(chr(i))
    return labels


# from https://github.com/lex8erna/UPGMApy/blob/master/UPGMA.py
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

    assert upgma(dist_matrix, M_labels) == '((((A,D),((B,F),G)),C),E)'


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
    test_upgma()


######################## F
def read_fasta():
    seq_lst = []
    name_lst = []
    record_iterator = SeqIO.parse("sequences.fasta", "fasta")
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

@timeit
def run_kmer(seq_lst, name_lst):
    distance_matrix_kmer = kmer_dist(seq_lst, 3)
    # distance_matrix_kmer = matrix(len(dist_kmer), len(dist_kmer[0]))
    # distance_matrix_kmer.set_whole_matrix(dist_kmer)
    print("distance_matrix_kmer")
    print(distance_matrix_kmer)
    print("*****************kmer************************")
    result_of_upgma_kmer = upgma(distance_matrix_kmer, name_lst)
    content_kmer = f"result_of_upgma_kmer:\n {result_of_upgma_kmer}"
    write_to_file("results/result_kmer", content_kmer)

@timeit
def run_global(seq_lst, name_lst):
    distance_matrix_globalpw = globalpw_dist(seq_lst)
    # distance_matrix_globalpw = matrix(len(dist_global), len(dist_global[0]))
    # distance_matrix_globalpw.set_whole_matrix(dist_kmer)
    print("distance_matrix_globalpw")
    print(distance_matrix_globalpw)
    print("*****************globalpw************************")
    result_of_upgma_globalpw = upgma(distance_matrix_globalpw, name_lst)
    content_global = f"result_of_upgma_globalpw:\n {result_of_upgma_globalpw}\n"
    write_to_file("results/result_global", content_global)
@timeit
def f_question():
    seq_lst, name_lst = read_fasta()
    run_kmer(seq_lst, name_lst)
    run_global(seq_lst, name_lst)

def current_time():
    now = datetime.now()
    dt_string = now.strftime("%d-%m-%Y--%H-%M-%S")
    return dt_string

def write_to_file(file, content):
    with open(f"{file}_{current_time()}.txt", "w") as f:
        f.write(content)
    return file

####################### end F

if __name__ == '__main__':
    # closed_tests()
    f_question()
    pass
