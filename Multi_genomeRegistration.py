
import argparse
import csv
import pandas as pd
import numpy
from tqdm import tqdm_notebook
import os
import gzip
import time

#
# parser = argparse.ArgumentParser()
# #parser.add_argument('-Let', '--Letter',help='chromosome used to do search',required=True)
# parser.add_argument('-strain', '--strain', help='strain used to do search', required=True)
# args = parser.parse_args()
# #Letter = args.Letter
# strain = args.strain
#
# print strain
# #print Letter
#
# #letter = Letter


import CCGG_SH as CCGG
import csv
import json
import numpy
import time
import pandas as pd
from tqdm import tqdm_notebook
from collections import defaultdict
import matplotlib.pyplot as plot



# In[2]:


import json
def load_anchorcandidate_dict(filename):
    '''load the anchor candidate dictionary
    Output {anchor_name:sequence}'''
    with open(filename) as json_file: # anchor : sequence
        anchor_Dict = json.load(json_file)
    return anchor_Dict


# In[3]:


def create_kmer_profile(anchor_Dict):
    """Input: Anchor Candidates Dictionary (anchor_name: sequence)
    Output: initialize empty kmer dictionary for recording mapping positions of both forward and reverse complement
    (sequence: [] )
    """
    KmerProfile = defaultdict(lambda: defaultdict(list))
    for anchor,seq in anchor_Dict.iteritems():
        KmerProfile[seq]
        rev = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(seq)])
        KmerProfile[rev]
    # print len(KmerProfile)
    return KmerProfile


# In[4]:


def mapping(filename, chromo, Position_Dict, k):

    #filename = '/csbiodataxw/KeaneGenomes/genomes/DBA2/Chr%s.seq' % str(chromo)
    with open(filename, "r") as f:
        seq = f.read().upper()
    print 'chromosome', chromo
    if chromo == 'X':
        chromo = 20

    for i in tqdm_notebook(range(1, len(seq) - k + 1)):
        kmer = seq[i:i+45]
        if 'N' in kmer:
            continue

        c = Position_Dict.get(kmer, -1)
        if c < 0:
            continue
        Position_Dict[kmer]['chr' + str(chromo)] += [i]
    return Position_Dict


# In[5]:


def get_mapping_position(anchor_Dict, s, k = 45):
    letter2strain = {'C':'129S1', 'A':'AJ', 'B': 'B6', 'F': 'CAST', "I":'DBA2', 'D': 'NOD', 'E':'NZO','G':'PWK', 'H': "WSB"}
    strain = letter2strain[s]
    kmer_d = create_kmer_profile(anchor_Dict) #initialize the profile

    chrlist = range(1,20) + ['X']
    # mapping
    for chromo in chrlist:
        filename = '/csbiodataxw/KeaneGenomes/genomes/%s/Chr%s.seq' % (strain, str(chromo))
        kmer_d = mapping(filename, chromo, kmer_d, k)
    return kmer_d


# In[6]:


#filter anchors
def filter_anchors(anchor_Dict,Kmer_Dict):
    '''Given kmer mapping position, return the mapping position of each anchor'''
    Duplication = {}
    Deletion = []
    Inversion = {}
    wrong_anchor = []

    stime = time.time()
    for anchor, kmer in anchor_Dict.iteritems():
        count = 0
        if len(Kmer_Dict[kmer]) < 1:
            wrong_anchor.append(anchor) # anchor forward sequence do not occur in linear assembly

            # check reverse complement occurrence
            krev = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(kmer)])
            for chromo, poslist in Kmer_Dict[krev].iteritems():
                count += len(poslist)
            # classify
            if count >0:
                Inversion[anchor] = Kmer_Dict[krev]
            else:
                Deletion.append(anchor)

        else:
            # check forward occurrence
            for chromo, poslist in Kmer_Dict[kmer].iteritems():
                count += len(poslist)
            # check reverse complement occurrence
            krev = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(kmer)])
            for chromo, poslist in Kmer_Dict[krev].iteritems():
                count += len(poslist)

            # classify
            if count > 1:
                wrong_anchor.append(anchor)
                Duplication[anchor] = Kmer_Dict[kmer]
                Duplication[anchor].update(Kmer_Dict[krev])
    print "Using %s sec" % str(time.time() - stime)
    return Duplication, Inversion, Deletion, wrong_anchor


# In[7]:


def write_file(strain, Duplication, Inversion, Deletion):
    fileName = "./SV_anchor_info/%s_Deletion_anchor.csv" % strain
    with open(fileName, "w") as csvfile:
        csvwriter = csv.writer(csvfile,  delimiter=',')
        csvwriter.writerow(["Deletion"])
        for anchor in sorted(Deletion):
            csvwriter.writerow([anchor])

    fileName = "./SV_anchor_info/%s_Duplication_anchor.json" % strain
    with open(fileName, 'w') as fp:
        json.dump(Duplication, fp)

    fileName = './SV_anchor_info/%s_Inversion_anchor.json' % strain
    with open(fileName, 'w') as fp:
        json.dump(Inversion, fp)


# In[8]:


def kmer_write_file(strain, Duplication, Inversion, Deletion):
    fileName = "./SV_anchor_info/%s_Deletion_kmer.csv" % strain
    with open(fileName, "w") as csvfile:
        csvwriter = csv.writer(csvfile,  delimiter=',')
        csvwriter.writerow(["Deletion"])
        for anchor in sorted(Deletion):
            csvwriter.writerow([anchor])

    fileName = "./SV_anchor_info/%s_Duplication_kmer.json" % strain
    with open(fileName, 'w') as fp:
        json.dump(Duplication, fp)

    fileName = './SV_anchor_info/%s_Inversion_kmer.json' % strain
    with open(fileName, 'w') as fp:
        json.dump(Inversion, fp)


# In[9]:


def get_anchorposlist(anchor_Dict, wrong_anchor, Kmer_Dict):
    translocation = {}
    Original_Anchors = anchor_Dict.keys()
    remaining_anchor = sorted(set(Original_Anchors) - set(wrong_anchor))
    PositionDict = defaultdict(list)
    AnchorDict = defaultdict(list)

    for anchor in tqdm_notebook(remaining_anchor):
        seq = anchor_Dict[anchor]
        l = Kmer_Dict[seq].items()
        chromosome, poslist =  l[0]
        assert len(poslist) == 1
        chromo = "chr" + str(int(anchor[1:3]))

        if chromo != chromosome:
            translocation[anchor] = l
        else:
            PositionDict[chromosome] += poslist
            AnchorDict[chromosome] += [anchor]
    return translocation, PositionDict, AnchorDict


# In[10]:


def binary_search(arr, val, l, r):
    if l == r:
        if arr[l] > val:
            return l
        else:
            return l+1
    if l > r:
        return l

    mid = (l+r)/2
    if arr[mid] < val:
        return binary_search(arr, val, mid+1, r)
    elif arr[mid] > val:
        return binary_search(arr, val, l, mid-1)
    else:
        return mid

# NlogN
from tqdm import tqdm_notebook
def efficientDeletionSort(array):
    subsequence_end = [0] # index of the element
    predecessors = [-1] # predecessor index
    for i in tqdm_notebook(range(1,len(array))):
        arr = array[i]
        # can do binary search instead, just skip to make it faster
        if arr > array[subsequence_end[-1]]:
            predecessors += [subsequence_end[-1]]
            subsequence_end += [i]
        else:
            # preform binary search
            minimum_end = [array[j] for j in subsequence_end] # element in current subsequence

            insert_point = binary_search(minimum_end, arr, 0, len(minimum_end)-1)
            if insert_point > 0:
                predecessors += [subsequence_end[insert_point-1]]
            else:
                predecessors += [-1]

            if insert_point > len(subsequence_end)-1: # arr == array[subsequence_end[-1]]
                subsequence_end += [i]
            elif arr < array[subsequence_end[insert_point]]:
                subsequence_end[insert_point] = i

    # backtrack
    pre = subsequence_end[-1]
    listIndex = []
    while pre != -1:
        listIndex.append(pre)
        pre = predecessors[pre]
    listIndex.reverse()
    longest_subsequence = [array[i] for i in listIndex]
    return listIndex, longest_subsequence


# In[11]:


def filter_misordered_anchor(PositionDict, AnchorDict):
    Misordered = defaultdict(list)
    Current_anchor = defaultdict(list)
    for chromo in range(1,21):
        chromosome = 'chr' + str(chromo)
        print chromosome
        poslist = PositionDict[chromosome]
        anchorlist = AnchorDict[chromosome]
        listIndex, longest_subsequence = efficientDeletionSort(poslist)
        for i in tqdm_notebook(listIndex):
            Current_anchor[chromosome] += [(anchorlist[i], poslist[i])]

        # misindex
        misindex = list(set(range(len(poslist))) - set(listIndex))
        print len(misindex)
        for index in misindex:
            Misordered[chromosome] += [(anchorlist[index], poslist[index])]
    return Misordered, Current_anchor


# In[12]:


def main(k,strain, acfilename, genomepath):
    chrlist = range(1,20) + ['X']
    anchor_Dict = load_anchorcandidate_dict(acfilename)
    #k = 45
    # original = 0
    # final = 0
    # Wrong_anchor = set()
    demote_anchor = 0
    kmer_d = create_kmer_profile(anchor_Dict)
    print strain, kmer_d.items()[0]
    # mapping
    for chromo in chrlist:
        filename = genomepath + '/%s/Chr%s.seq' % (strain, str(chromo))
        kmer_d = mapping(filename, chromo, kmer_d, k)

    # filter
    Duplication, Inversion, Deletion, wrong_anchor = filter_anchors(anchor_Dict,kmer_d)
    print len(Duplication), len(Inversion), len(Deletion), len(wrong_anchor)
    demote_anchor += len(Duplication)+ len(Inversion)+ len(Deletion)

    write_file(strain, Duplication, Inversion, Deletion)
    translocation, PositionDict, AnchorDict = get_anchorposlist(anchor_Dict, wrong_anchor, kmer_d)
    demote_anchor += len(translocation)
    t = 0
    for c in PositionDict.keys():
        t += len(PositionDict[c])
    assert demote_anchor + t == len(anchor_Dict)
    print len(translocation), len(PositionDict)
    fileName = './SV_anchor_info/%s_Translocation_anchor.json' % strain
    with open(fileName, 'w') as fp:
        json.dump(translocation, fp)

    # monotonicity
    Misordered, Current_anchor = filter_misordered_anchor(PositionDict, AnchorDict)
    m = 0
    for c in Misordered.keys():
        m += len(Misordered[c])
    demote_anchor += m
    r = 0
    for c in Current_anchor.keys():
        r += len(Current_anchor[c])

    assert demote_anchor + r == len(anchor_Dict)
    print len(Misordered), len(Current_anchor)
    fileName = './SV_anchor_info/%s_Misordered_anchor.json' % strain
    with open(fileName, 'w') as fp:
        json.dump(Misordered, fp)
    #final anchor
    fileName = './anchor_info/%s_Anchor_pos.json' % strain
    with open(fileName, 'w') as fp:
        json.dump(Current_anchor, fp)

    del kmer_d
    del Duplication
    del Inversion
    del Deletion
    del wrong_anchor
    del Misordered
    del Current_anchor
    del translocation
    del PositionDict
    del AnchorDict


# In[ ]:

# Strainlist = ['AJ', 'CAST', 'DBA2', 'NOD', 'NZO',  'PWK',  'WSB']
# for strain in Strainlist:
#     main(45, strain)


# In[12]:


def SV_search_main(k, strain):
    chrlist = range(1,20) + ['X']
    anchor_Dict = load_anchorcandidate_dict()
    #k = 45
    # original = 0
    # final = 0
    # Wrong_anchor = set()
    demote_anchor = 0
    kmer_d = create_kmer_profile(anchor_Dict)
    print strain, kmer_d.items()[0]
    # mapping
    for chromo in chrlist:
        filename = '/csbiodataxw/KeaneGenomes/genomes/%s/Chr%s.seq' % (strain, str(chromo))
        kmer_d = mapping(filename, chromo, kmer_d, k)

    # filter
    Duplication, Inversion, Deletion, wrong_anchor = filter_anchors(anchor_Dict,kmer_d)
    print len(Duplication), len(Inversion), len(Deletion), len(wrong_anchor)
    demote_anchor += len(Duplication)+ len(Inversion)+ len(Deletion)

    fileName = "./SV_anchor_info/%s_Duplication_anchor.json" % strain
    with open(fileName, 'w') as fp:
        json.dump(Duplication, fp)

    # write_file(strain, Duplication, Inversion, Deletion)
    # translocation, PositionDict, AnchorDict = get_anchorposlist(anchor_Dict, wrong_anchor, kmer_d)
    # demote_anchor += len(translocation)
    # t = 0
    # for c in PositionDict.keys():
    #     t += len(PositionDict[c])
    # assert demote_anchor + t == len(anchor_Dict)
    # print len(translocation), len(PositionDict)
    # fileName = './SV_anchor_info/%s_Translocation_anchor.json' % strain
    # with open(fileName, 'w') as fp:
    #     json.dump(translocation, fp)
    #
    # monotonicity
    # Misordered, Current_anchor = filter_misordered_anchor(PositionDict, AnchorDict)
    # m = 0
    # for c in Misordered.keys():
    #     m += len(Misordered[c])
    # demote_anchor += m
    # r = 0
    # for c in Current_anchor.keys():
    #     r += len(Current_anchor[c])
    #
    # assert demote_anchor + r == len(anchor_Dict)
    # print len(Misordered), len(Current_anchor)
    # fileName = './SV_anchor_info/%s_Misordered_anchor.json' % strain
    # with open(fileName, 'w') as fp:
    #     json.dump(Misordered, fp)
    # #final anchor
    # fileName = './anchor_info/%s_Anchor_pos.json' % strain
    # with open(fileName, 'w') as fp:
    #     json.dump(Current_anchor, fp)
    #
    # del kmer_d
    # del Duplication
    # del Inversion
    # del Deletion
    # del wrong_anchor
    # del Misordered
    # del Current_anchor
    # del translocation
    # del PositionDict
    # del AnchorDict
main(k,strain, acfilename, genomepath)
