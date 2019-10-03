# Parsing file

import csv
import math
import numpy


def tsv_to_d2_list(file_name):
        d2_list = []
        with open(file_name) as list_values:
                list_reader = csv.reader(list_values, delimiter='\t')
                for list_object in list_reader:
                        d2_list.append(list_object)
        return d2_list

