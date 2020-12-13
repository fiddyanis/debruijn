#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Anis Fiddy"
__copyright__ = "CY Tech"
__credits__ = ["Anis Fiddy"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Anis Fiddy"
__email__ = "fiddyanis@eisti.eu"
__status__ = "Finished"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file,'rt') as f:
        sentences = f.readlines()
        #print(sentences)
        for sentence in sentences:
            sentence = sentence.replace("\n", "")
            for letter in sentence:
                if letter in ['+','@','J']:
                    break
            else:
                yield sentence


def cut_kmer(read, kmer_size):
    for i in range(len(read)+1-kmer_size):
        yield(read[i:i+kmer_size])

def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    for sentence in read_fastq(fastq_file):
        
        for kmer in cut_kmer(sentence, kmer_size):
            if kmer in kmer_dict:
                kmer_dict[kmer] = kmer_dict.get(kmer,0) + 1
            else:
                
                kmer_dict[kmer] = 1
    return kmer_dict

def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for kmer in kmer_dict.keys():
        graph.add_weighted_edges_from([(kmer[:-1], kmer[1:], kmer_dict[kmer])])
    return graph

def get_starting_nodes(graph):
    return [n for n,d in graph.in_degree() if d==0]

def get_sink_nodes(graph):
    return [node for node,out_degree in graph.out_degree() if out_degree == 0]

def get_contigs(graph, starting_nodes, ending_nodes):
    contigsRes = []
    for entry_node in starting_nodes:
        for end_node in ending_nodes:
            allPaths = list(nx.all_simple_paths(graph, entry_node, end_node))
            if len(allPaths) > 0:
                path = allPaths[0][0] + "".join([allPaths[0][i][-1] for i in range(1,len(allPaths[0]))])
                contigsRes.append((path, len(path)))
    return contigsRes

def fill(text, width=80):
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contigs_list, output_file):
    with open(output_file, "w") as f:
        for i in range(len(contigs_list)):
            f.write(">contig_"+str(i)+" len="+str(contigs_list[i][1])+"\n"+fill(contigs_list[i][0])+"\n")
        

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        for node in path[1 : -1]:
            graph.remove_node(node)
        if delete_sink_node:
            graph.remove_node(path[-1])
        if delete_entry_node:
            graph.remove_node(path[0])
    return graph

def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):

    if path_list == []:
        return graph
    bestId = 0
    for i in range(1,len(path_list)):
        if (weight_avg_list[bestId] == weight_avg_list[i] and path_length[bestId] < path_length[i]) or weight_avg_list[bestId] < weight_avg_list[i]:
            bestId = i
        elif path_length[bestId] == path_length[i]:
            bestId = random.choice([bestId, i])
    path_list.remove(path_list[bestId])
    return remove_paths(graph, path_list, delete_entry_node, delete_sink_node)


def path_average_weight(graph, path):
    weights = []
    length = len(path)-1
    for i in range(length):
        weights.append(graph[path[i]][path[i+1]]["weight"])
    #print(f"weights : {weights}")
    return sum(weights)/length


def solve_bubble(graph, ancestor_node, descendant_node):
    allPaths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    path_list = []
    path_length = []
    weight_avg_list = []
    
    for path in allPaths:

        avg = path_average_weight(graph, path)
        path_list.append(path)
        weight_avg_list.append(avg)
        path_length.append(len(path))
    return select_best_path(graph, path_list, path_length, weight_avg_list,False, False)
    

def simplify_bubbles(graph):
    for node1 in get_starting_nodes(graph):


        for node2 in get_sink_nodes(graph):
            graph =   solve_bubble(graph, node1, node2)

    return graph

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass
def main():
    """
    Main program function
    """
    print("1) Lecture du fichier et construction du graphe")
    args = get_arguments()
    
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)

    graph = build_graph(kmer_dict)

    print("2) Résolution des bulles")
    graph = simplify_bubbles(graph)
    print("3) Résolution des pointes d’entrée et de sortie")
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)

    print("4) Ecriture du/des contigs")
    contigs = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs, args.output_file)

    
if __name__ == '__main__':
    main()
