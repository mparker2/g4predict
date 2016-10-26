'''
Methods for filtering a set of overlapping bed records.

author: Matthew Parker
'''

from operator import itemgetter
from bisect import bisect_left


def check_overlapping(cluster_range, record_start):
    '''
    check if record y overlaps with cluster range x
    '''
    if record_start < cluster_range[1]:
        return True
    else:
        return False


def cluster_overlapping(bed):
    '''
    generator yielding new clusters, treats each strand separately
    '''
    cluster = {'+': [], '-': []}
    cluster_range = {'+': [0, 0], '-': [0, 0]}
    cluster_chrom = {'+': None, '-': None}
    while True:
        try:
            record = next(bed).split()
        except StopIteration:
            if cluster['+']:
                yield cluster['+']
            if cluster['-']:
                yield cluster['-']
            break

        record[1] = int(record[1])
        record[2] = int(record[2])
        record[4] = float(record[4])
        s = record[5]

        # check cluster not empty:
        if not cluster[s]:
            cluster[s].append(record)
            cluster_range[s] = [record[1], record[2]]
            cluster_chrom[s] = record[0]

        # check we are still on the same chromosome
        elif cluster_chrom[s] != record[0]:
            yield cluster[s]
            cluster[s] = [record, ]
            cluster_range[s] = [record[1], record[2]]
            cluster_chrom[s] = record[0]

        # check if record overlaps with cluster
        elif check_overlapping(cluster_range[s], record[1]):
            cluster[s].append(record)

            # if neccessary, increase cluster range
            if record[2] > cluster_range[s][1]:
                cluster_range[s][1] = record[2]

        # if record does not overlap with cluster, yield the cluster
        # and start a new one
        else:
            yield cluster[s]
            cluster[s] = [record, ]
            cluster_range[s] = [record[1], record[2]]


def join_records(cluster):
    '''
    list of lists to list of tab delim strings
    '''
    return ['\t'.join([str(f) for f in record]) for record in cluster]


def filter_overlapping(cluster):
    '''
    find the non-overlapping records in a cluster which yield
    the best total score.
    '''
    # if cluster is only one record, return it
    if len(cluster) == 1:
        return join_records(cluster)
    # if cluster is only two records, return higher scoring
    if len(cluster) == 2:
        return join_records([max(cluster, key=itemgetter(4)), ])
    # cluster is sorted by stop-values
    end_sorted_cluster = sorted(cluster, key=itemgetter(2))
    end_sorted_vals = [x[2] for x in end_sorted_cluster]

    # vector of maximum scores, one longer than cluster size
    max_of_cluster = [[0, 0, 0], ] * (len(cluster) + 1)

    # calculate max score from cluster:
    # for each record in the cluster:
    for i, record in enumerate(end_sorted_cluster):

        # calculate the score from prev records if not including this record
        not_incl = max(max_of_cluster[i][:2])

        # find the closest non-overlapping record left of the current one
        closest = bisect_left(end_sorted_vals, record[1])
        if closest == 0:
            closest = -1

        # score including this record is total from closest plus current
        incl = max(max_of_cluster[closest][:2]) + record[4]

        # store scores and closest record in matrix
        max_of_cluster[i+1] = [not_incl, incl, closest]
    # backtrack through matrix to get high scoring records:
    i = len(cluster)
    incl_records = []
    while i != -1:
        # if not_incl score is less than incl
        if max_of_cluster[i][1] >= max_of_cluster[i][0]:

            # include cluster[i]
            incl_records.append(end_sorted_cluster[i-1])

            # jump to closest
            i = max_of_cluster[i][2]
        else:

            # try next record
            i -= 1

    return join_records(incl_records)


def merge_overlapping(cluster):
    '''
    flatten all the records in a cluster into one bed6 formatted record.
    score is set as the total number of records in the cluster (may give
    some indication of G4 variants at the locus, though the regex does
    not capture all potential variants, just one per left mapping pos).
    '''
    cluster_min = min(cluster, key=itemgetter(1))[1]
    cluster_max = max(cluster, key=itemgetter(2))[2]
    score = len(cluster)
    return join_records([[cluster[0][0], cluster_min, cluster_max,
                          'PG4_cluster', score, cluster[0][5]]])


def apply_filter_method(file_handle, filter_method):
    '''
    Cluster overlapping records then apply merge or filter methods.
    '''
    for cluster in cluster_overlapping(file_handle):
        for record in filter_method(cluster):
            yield record
