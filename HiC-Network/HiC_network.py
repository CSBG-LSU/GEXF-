import pandas as pd
import networkx as nx
import numpy as np
import argparse
import time

def getArgs():
   parser = argparse.ArgumentParser('python')
   parser.add_argument('-inter',required=True,help='interchromosomal contacts input file')
   parser.add_argument('-intra',required=True,help='intrachromosomal contacts input file')
   parser.add_argument('-bias',required=True,help='bias input file')
   parser.add_argument('-output',required=True,help='network output file')
   return parser.parse_args()


def HIC_network(intercontacts, intracontacts,bias, output):

    interchromosomal_contacts = pd.read_csv(intercontacts, delim_whitespace=True)
    intrachromosomal_contacts = pd.read_csv(intracontacts, delim_whitespace=True)
    all_contacts = intrachromosomal_contacts.append(interchromosomal_contacts, ignore_index=True)
    all_contacts['InteractorA'] = all_contacts[['chr1', 'chunk1_start', 'chunk1_end']].astype(str).apply('-'.join, 1)
    all_contacts['InteractorB'] = all_contacts[['chr2', 'chunk2_start', 'chunk2_end']].astype(str).apply('-'.join, 1)
    all_contacts['p-value'] = all_contacts['p-value'].apply(np.log10)
    all_contacts['p-value'] = - all_contacts['p-value']
    all_contacts['q-value'] = all_contacts['q-value'].apply(np.log10)
    all_contacts['q-value'] = - all_contacts['q-value']
    m = all_contacts.loc[all_contacts['p-value'] != np.inf, 'p-value'].max()
    n = all_contacts.loc[all_contacts['q-value'] != np.inf, 'q-value'].max()
    all_contacts['p-value'] = all_contacts['p-value'].replace(np.inf, m+1)
    all_contacts['q-value'] = all_contacts['q-value'].replace(np.inf, n+1)
    unique_fragments = np.unique(all_contacts[['InteractorA', 'InteractorB']].values)
    list_of_ufragments = unique_fragments.tolist()
    fragments_dataframe = pd.DataFrame(list_of_ufragments)
    fragments_dataframe['index'] = np.arange(1, len(fragments_dataframe) + 1)
    fragments_dataframe.columns = ['Interactor', 'Nodes']
    mapping_file = fragments_dataframe[['Interactor', 'Nodes']]
    bias_file = pd.read_csv(bias, delim_whitespace=True)
    bias_file['Interactor'] = bias_file[['chr', 'chunk_start', 'chunk_end']].astype(str).apply('-'.join, 1)
    H = mapping_file['Nodes'].tolist()
    lol = mapping_file.groupby('Interactor').Nodes.apply(list).to_dict()
    all_contacts['InteractorA'] = all_contacts['InteractorA'].astype(str)
    all_contacts['InteractorB'] = all_contacts['InteractorB'].astype(str)
    all_contacts['Nodes1'] = all_contacts['InteractorA'].map(lol)
    all_contacts['Nodes2'] = all_contacts['InteractorB'].map(lol)
    all_contacts['Nodes1'] = all_contacts['Nodes1'].str.get(0)
    all_contacts['Nodes2'] = all_contacts['Nodes2'].str.get(0)
    q = all_contacts[['Nodes1', 'Nodes2', 'contactCount','p-value', 'q-value']]
    G = nx.from_pandas_edgelist( q, 'Nodes1', 'Nodes2', ['contactCount', 'p-value', 'q-value'])
    chr_fragments = [i.split('-', 1)[0] for i in list_of_ufragments]
    chunk_fragments = [i.split('-', 1)[1] for i in list_of_ufragments]
    chunk_start_fragments = [i.split('-', 1)[0] for i in chunk_fragments]
    chunk_end_fragments = [i.split('-', 1)[1] for i in chunk_fragments]
    node_attr1 = dict(zip(H, chr_fragments))
    nx.set_node_attributes(G, node_attr1, 'chr')
    node_attr2 = dict(zip(H, chunk_start_fragments))
    nx.set_node_attributes(G, node_attr2, 'chunk_start')
    node_attr3 = dict(zip(H, chunk_end_fragments))
    nx.set_node_attributes(G, node_attr3, 'chunk_end')
    nx.write_gexf(G, output)

if __name__ == "__main__":   
   args = getArgs()    
   network = HIC_network(args.inter, args.intra, args.bias, args.output)  
   start = time.time()
   end = time.time()   
   print ('time elapsed:' + str(end - start))

