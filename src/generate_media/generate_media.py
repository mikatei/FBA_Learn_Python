#This file contains a function that generates media based on the input.

#from MyOtherModule.parsing import convert_d2_list_to_tab_separated_string
from parsing_functions import tsv_to_d2_list
from parsing import string_to_file, convert_d2_list_to_tab_separated_string 
import logging
import math
import os
'''
logger = logging.getLogger('debuger')
logger.setLevel(logging.DEBUG)
'''
logging.basicConfig(level=logging.DEBUG)




#The following takes a list of all compounds of interest and returns the permutations of media with and without that compound.
def media_permutations(cmpnd_list_d2):
    if len(cmpnd_list_d2) > 12:
        raise Exception("too many compounds - max compounds should be 12. Currently, the total number of media is " + str(len(cmpnd_list_d2)))
    else:
        #medial_list_d3 holds the different cmpnd_list_d2s which each contain information for a different media file.
        media_list_d3 = []
        current_compound_list = []
        recursion_on_cmpnds(cmpnd_list_d2, media_list_d3, current_compound_list, 0)
        return media_list_d3


def recursion_on_cmpnds(cmpnd_list_d2, media_list_d3, current_compound_list, i):
    if i == len(cmpnd_list_d2):
        media_list_d3.append(current_compound_list)
    else:
        crnt_cmpnd = cmpnd_list_d2[i]
        new_compound_list_with = current_compound_list[:] + [crnt_cmpnd]
        recursion_on_cmpnds(cmpnd_list_d2,media_list_d3,new_compound_list_with, i+ 1)
        new_compound_list_without = current_compound_list[:]
        recursion_on_cmpnds(cmpnd_list_d2,media_list_d3,new_compound_list_without, i + 1)

def only_cmpnd_ids(cmpnd_list_d2):
    cmpnd_ids_d1 = []
    for cmpnd in cmpnd_list_d2:
        cmpnd_ids_d1.append(cmpnd[0])
    return cmpnd_ids_d1

def list_of_media_ids(media_list_d3):
    media_ids_d2 = []
    for cmpnd_list_d2 in media_list_d3:
        media_ids_d2.append(only_cmpnd_ids(cmpnd_list_d2)
    return media_ids_d2

#The input to this function is a list: [[cpdid1,name1,formula1, minFlux1,maxFlux1,concentration1],[cpdid2,name2,...],...]
def make_media_tsv_str(cmpnd_list_d2):

    header_list = ['compounds','name','formula','minFlux','maxFlux','concentration']
    final_compound_list_d2 = [header_list] + cmpnd_list_d2
    media_tsv_string = convert_d2_list_to_tab_separated_string(final_compound_list_d2)

    return media_tsv_string


def media_list_to_dir_files(media_list_d3, dir_path): 
    
    for i in range(len(media_list_d3)):
        crnt_media_d2 = media_list_d3[i]
        crnt_media_string = make_media_tsv_str(crnt_media_d2)
        new_file_path = os.path.join(dir_path,'media' + str(i) + '.tsv')
        string_to_file(crnt_media_string,new_file_path)



def test():
    '''
    cmpnd_list_d2 = [['a',1],['b',2],['c',3]]
    media_list_d3 = media_permutations(cmpnd_list_d2)
    print(media_list_d3)
    '''
    dir_path = '/Users/omreeg/KBase/apps/omreegalozMediaPermutations/data/tmpmedia' 
    filename = '/Users/omreeg/KBase/apps/omreegalozMediaPermutations/data/Example_Data/new_compounds.tsv'
    cmpnd_list_d2 = tsv_to_d2_list(filename)[1:]
    logging.debug(cmpnd_list_d2)
    media_list_d3 = media_permutations(cmpnd_list_d2)
    #media_list_to_dir_files(media_list_d3,dir_path)
    list_of_media_ids
    #print(media_list_d3)


def main():
    test()

main()

