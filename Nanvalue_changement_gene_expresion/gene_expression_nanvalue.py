
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
    nodesNAN = [x for x, y in Graph.nodes(data=True) if y['Gene-exp'] == 'nan']  # find nodes Gene-exp=nan
    num = len(nodesNAN)
    for i in range(num):
        nan_neighbors = [n for n in Graph[nodesNAN[i]]]
        j = len(nan_neighbors)
        if j == 1:
            a = (Graph.nodes[nan_neighbors[0]]['Gene-exp'])
            Graph.nodes[nodesNAN[i]]['Gene-exp'] = a
    nx.write_gexf(Graph, outputfile)


if __name__ == "__main__":
    args = getArgs()
    create_network(args.inputfile, args.outputfile)
    start = time.time()
    end = time.time()

    print('time elapsed:' + str(end - start))
