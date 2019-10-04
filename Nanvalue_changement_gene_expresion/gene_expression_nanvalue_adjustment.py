"""
Version: 1.0
Author: Guannan
"""
# -*- coding: utf-8 -*-
import pandas as pd
import networkx as nx
import numpy as np
import argparse
import time

def getArgs():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-inputfile', required=True)
    parser.add_argument('-outputfile', required=True)
    return parser.parse_args()


def create_network(inputfile, outputfile):
    Graph = nx.read_gexf(inputfile)
    nodesNAN = [x for x, y in Graph.nodes(data=True) if np.isnan(y['gene_exp'])]
    num = len(nodesNAN)
    for i in range(num):
        nan_neighbors = [n for n in Graph[nodesNAN[i]]]
        num_nan_neighbors = len(nan_neighbors)
        temp2 = ['0.0'] * num_nan_neighbors
        for j in range(num_nan_neighbors):
            temp1 = Graph.nodes[nan_neighbors[j]]['gene_exp']
            temp2[j] = temp1
        sum1 = np.sum(temp2)
        if sum1 > 0.0:
            sum1 = 1.0
        elif sum1 < 0.0:
            sum1 = -1.0
        else:
            sum1 = 0.0
        Graph.nodes[nodesNAN[i]]['gene_exp'] = sum1
    nx.write_gexf(Graph, outputfile)


if __name__ == "__main__":
    args = getArgs()
    create_network(args.inputfile, args.outputfile)
    start = time.time()
    end = time.time()

    print('time elapsed:' + str(end - start))
