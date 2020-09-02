# -*- coding: utf-8 -*-
"""
@Created at 2020/8/14 16:57
@Author: Kurt
@file:DAG.py
@Desc:
"""
from random import randint, random
import networkx as nx
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import matplotlib.pyplot as plt


def gen_DAG(min_width=3, max_width=5, min_length=1, max_length=5,
            percent=0.3, save="DAG.txt"):
    """
    Randomly generate a DAG type supply chain
    :param min_width: minimal number of nodes of each echelon
    :param max_width: maximal number of nodes of each echelon
    :param min_length: minimal length of the supply chain
    :param max_length: maximal length of the supply chain
    :param percent: probability of generating a link
    :param save: file path
    :return: number of nodes. num of edges
    todo: fix the bug that the disconnected graph is generated when probability is small. Hint: make sure more than one link between two lays
    """
    nodes = min_width + randint(0, max_width - min_width)
    print("initial nodes: {}".format(nodes))
    heights = min_length + randint(0, max_length - min_length)
    print("Heights: {}".format(heights))
    num_of_edges = 0
    node_set = set()
    with open(save, 'w') as f:
        for i in range(0, heights + 1):  # no action in loop 0, so increase 1.
            new_nodes = min_width + randint(0, max_width - min_width)
            # print("Weights: {}，loop {}".format(new_nodes, i))
            print("Gen new nodes: {}".format(new_nodes))
            for j in range(0, nodes):
                for k in range(0, new_nodes):
                    if random() < percent:
                        num_of_edges += 1
                        node_set.add(j)
                        node_set.add(k+nodes)
                        print("{}--->{}".format(j, k + nodes))
                        f.write(str(j))
                        f.write("\t")
                        f.write(str(k + nodes))
                        f.write('\n')
            nodes = nodes + new_nodes
    num_of_nodes = len(node_set)
    return num_of_nodes, num_of_edges


def vis_graph(file):
    G = nx.DiGraph()
    with open(file, 'r') as f:
        links = [line.strip().split('\t') for line in f]
        G.add_edges_from(links)

    write_dot(G, 'DAG.dot')
    pos = graphviz_layout(G, prog='dot')
    nx.draw(G, pos, with_labels=True, font_weight='bold')
    plt.show()


if __name__ == "__main__":
    num_of_nodes, node_of_edges = gen_DAG(min_width=200, max_width=300, min_length=5, max_length=7,
                                          percent=0.05, save="DAG.txt")
    print("Number of nodes: {}, number of edges: {}".format(num_of_nodes, node_of_edges))
    # vis_graph("DAG.txt")
