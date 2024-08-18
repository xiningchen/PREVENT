import numpy as np
import networkx as nx
import importlib
import Towlson_group_code.controllability as ctrb
import pandas as pd
import pickle as pkl
import statistics as stats
from scipy.stats import mannwhitneyu, linregress
import os

importlib.reload(ctrb)

TIMES = ['bl', 'Y1', 'Y3', 'Y5']
PICKLE_PATH = '../../PREVENT_Study/pickles/'
FIGURE_PATH = '../../PREVENT_Study/figures/'
ICN_COLORS = {'auditory': 'lightpink',
              'dorsal attention': 'green',
              'somatosensory': 'cornflowerblue',
              'fronto-parietal': 'orange',
              'other': 'grey',
              'default mode': 'firebrick',
              'visual': 'darkorchid',
              'cingulo-opercular': 'mediumturquoise',
              'ventral attention': "magenta",
              'somatosensory-motor': "royalblue"}


def rank_nodes(G, attr='weight'):
    """
    Given a graph G and some nodal attribute, rank all nodes from highest rank (n) to smallest rank (0).
    Store the rank value as a new node attribute. If the nodal attribute contains 0s, then there is no specific order in
    to rank these nodes. This will cause the ranking to fail. All such nodes will be ranked with 0.
    :param G: Networkx graph object
    :param attr: Nodal attribute key
    :return: whether the ranking succeeded (boolean); the updated graph G; if ranking failed, a list of nodes where the
    rank could not be assigned.
    """
    if attr == 'weight':
        seq = sorted([(node, G.degree(node, weight=attr)) for node in G], reverse=True, key=lambda x: x[1])
    else:
        seq = sorted(G.nodes(data=attr), reverse=True, key=lambda x: x[1])
    rank = np.arange(0, len(seq))
    rank = rank[::-1]
    nodeRank = {}
    success = True
    badRegions = []
    attrTag = attr + 'Rank'
    offset = 0
    for i, e in enumerate(seq):
        if e[1] == 0:
            badRegions.append(e[0])
            nodeRank[e[0]] = {attrTag: 0}
            offset += 1
            success = False
            continue
        nodeRank[e[0]] = {attrTag: rank[i-offset]}

    nx.set_node_attributes(G, nodeRank)
    return success, G, badRegions


def load_meta_data():
    """
    Get meta data of participants in the study.
    :return: meta data frame, list of network region names
    """
    meta_data = pd.read_excel('../../PREVENT_Study/data/ages_and_diagnoses.xlsx', index_col=None, header=0)
    ids = []
    for i in meta_data['Subject ID']:
        id = str(i)
        if len(id) < 3:
            id = '0'*(3-len(id))+id
        ids.append(id)
    meta_data['Subject ID'] = ids
    meta_data = meta_data.set_index(['Subject ID'])
    DIAGNOSIS = {'T': 'P', 'C': 'HC'}
    meta_data['C/T'] = meta_data['C/T'].replace(DIAGNOSIS)
    regionList = pd.read_excel('../../PREVENT_Study/data/NodeList.xlsx', index_col=None, header=None)
    return meta_data, list(regionList[0])


def get_icn_map():
    """
    Using a pre-defined node to FN spreadsheet, create a pickle of node names to ICN mapping. Only need to run once for
    each time the spreadsheet updates/changes.
    :return: ../PREVENT_study/pickles/node_icn_map.pkl
    """
    if os.path.exists('../../PREVENT_study/pickles/node_icn_map.pkl'):
        with open('../../PREVENT_study/pickles/node_icn_map.pkl', 'rb') as f:
            node_icn_map = pkl.load(f)
            icn_node_map = pkl.load(f)
        return node_icn_map, icn_node_map
    else:
        fn_map = pd.read_excel('../../PREVENT_study/data/Functional_networks_anatomy - clean.xlsx', usecols='A:C',
                               header=0, index_col=None)
        node_icn_map = {}
        icns = set()
        for row in fn_map.iterrows():
            row = row[1]
            nn = row['Short Name'].replace("&", "_and_")
            networks = row['Network(s)'].split(",")
            networks_list = [n.strip().lower() for n in networks]
            icns.update(networks_list)
            if '_' in nn:
                node_icn_map[f'ctx_lh_{nn}'] = networks_list
                node_icn_map[f'ctx_rh_{nn}'] = networks_list
            else:
                node_icn_map[f'Right-{nn}'] = networks_list
                node_icn_map[f'Left-{nn}'] = networks_list

        node_df = pd.read_excel('../PREVENT_study/data/RegionList.xlsx', usecols='A', header=None, index_col=None)
        node_list = list(node_df[0])
        dne = []
        for nn in node_list:
            if nn in node_icn_map.keys():
                continue
            else:
                dne.append(nn)
        print(f"There are {len(node_list)} Original node regions and {len(node_icn_map)} node-icn mappings.")
        print(f"Could not match {len(dne)} node regions. ")
        print(*dne, sep="\n")

        # icn to node map
        node_list = [n for n in node_list]
        icn_node_map = {f: [] for f in icns}
        skips = []
        for nn, icn_list in node_icn_map.items():
            if nn not in node_list:
                skips.append(nn)
                continue
            for f in icn_list:
                icn_node_map[f].append(nn)
        total_nodes = set()
        for f, nn in icn_node_map.items():
            total_nodes.update(nn)
        print(f"There are {len(total_nodes)} nodes. Should have {len(node_list)} nodes. ")
        print(f"Skipped {len(skips)} nodes. ")

        for s in skips:
            r = node_icn_map.pop(s, None)
            print(r)

        with open('../../PREVENT_study/pickles/node_icn_map.pkl', 'wb') as f:
            pkl.dump(node_icn_map, f)
            pkl.dump(icn_node_map, f)
        return node_icn_map, icn_node_map


def linear_model(data_df, xlabel, ylabel, remove_outliers=3):
    if remove_outliers > 0:
        # remove x-axis outliers
        y = list(data_df[ylabel])
        stddev = stats.stdev(y)
        mean = stats.mean(y)
        model_df = data_df[(data_df[ylabel] < mean + 3 * stddev) & (data_df[ylabel] > mean - 3 * stddev)]

        # remove y-axis outliers
        # x = list(data_df[ylabel])
        # stddev = stats.stdev(x)
        # print("y std: ", stddev)
        # model_df = model_df[model_df[ylabel] < 3 * stddev]
    else:
        model_df = data_df
    x_fit = np.array(model_df[xlabel])
    y_fit = np.array(model_df[ylabel])
    if len(x_fit) == 0 or len(y_fit) == 0:
        print(f"No data for {xlabel}, {ylabel}")
        return None
    slope, intercept, _, p, stderr = linregress(x_fit, y_fit)
    return slope, intercept, p, stderr, len(x_fit)


def get_avg_node_metric(G, node_list, metric):
    values = [G.nodes[node][metric] for node in node_list]
    return sum(values) / len(values)

# ---------------------------------------------- WRONG DEFINITION ---- SEE brain_network.py
# def get_icn_intra_connectivity(G, node_list):
#     stuff = nx.get_edge_attributes(G, "weight")
#     w = 0
#     num = 0
#     for edge, v in stuff.items():
#         if (edge[0] in node_list) and (edge[1] in node_list):
#             w += v
#             num += 1
#     return w / num
#
#
# def get_icn_inter_connectivity(G, node_list1, node_list2):
#     stuff = nx.get_edge_attributes(G, "weight")
#     w = 0
#     num = 0
#     for edge, v in stuff.items():
#         if (edge[0] in node_list1) and (edge[1] in node_list2):
#             w += v
#             num += 1
#         elif (edge[1] in node_list1) and (edge[0] in node_list2):
#             w += v
#             num += 1
#     return w/num


def get_global_weight(G):
    stuff = nx.get_edge_attributes(G, "weight")
    return sum(stuff.values())/len(stuff)
# def load_cognitive_data():
#     """
#     Use:  metadata, cog_columns = load_cognitive_data()
#     """
#     aged = {'bl': 0, 'y1': 1, 'y3': 3, 'y5': 5}
#     metadata, _ = loadMetaFiles()
#     DIAGNOSIS = {'T': 'P', 'C': 'HC'}
#     metadata['C/T'] = metadata['C/T'].replace(DIAGNOSIS)
#     cog_columns = ['bvmt_total_recall', 'bvmt_delayed_recall', 'ravlt_list_a_delay_recall', 'tmt_trail_a_time', 'tmt_trail_b_time', 'wais_r_total', 'naart_total']
#     for time in ['bl', 'y1', 'y3', 'y5']:
#         behavior_pd = pd.read_excel('../PREVENT_Study/data/PREVENT_Data_Xining_Oct_2022.xlsx', sheet_name=time, index_col=None, header=0)
#         ids = []
#         ages = []
#         for i in behavior_pd['ID']:
#             id = str(i)
#             if len(id) < 3:
#                 id = '0'*(3-len(id))+id
#             ids.append(id)
#             ages.append(metadata['Age'][id] + aged[time])
#         times = [time for _ in range(len(behavior_pd['ID']))]
#         new_columns = {'Subject ID': ids, 'Age': ages, 'Time': times, 'C/T': list(metadata['C/T'])}
#         for c in cog_columns:
#             new_columns[c] = list(behavior_pd[c])
#
#         # new_pd = pd.DataFrame(new_columns)
#         if time == 'bl':
#             new_metadata = pd.DataFrame(new_columns)
#         else:
#             new_pd = pd.DataFrame(new_columns)
#             new_metadata = pd.concat([new_metadata, new_pd], axis=0)
#         new_metadata = new_metadata.reset_index(drop=True)
#     return new_metadata, cog_columns
