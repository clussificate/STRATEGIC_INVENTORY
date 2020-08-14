# -*- coding: utf-8 -*-
"""
@Created at 2020/8/14 16:57
@Author: Kurt
@file:DAG.py
@Desc:
"""
from random import randint
import networkx as nx
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import matplotlib.pyplot as plt


def gen_DAG(min_weight=3, max_weight=5, min_height=1, max_height=5,
            percent=30, save="DAG.txt"):
    """
    Randomly generate a DAG type supply chain
    :param min_weight: minimal number of nodes of each echelon
    :param max_weight: maximal number of nodes of each echelon
    :param min_height: minimal length of the supply chain
    :param max_height: maximal length of the supply chain
    :param percent: probability of generating a link
    :param save: file path
    :return: None
    """
    nodes = 0
    heights = min_height + randint(0, max_height - min_height)
    print("Heights: {}".format(heights))
    with open(save, 'w') as f:
        for i in range(0, heights + 1):  # no action in loop 0, so increase 1.
            new_nodes = min_weight + randint(0, max_weight - min_weight)
            print("Weights: {}，loop {}".format(new_nodes, i))

            for j in range(0, nodes):
                for k in range(0, new_nodes):
                    if randint(0, 100) < percent:
                        print("{}--->{}".format(j, k + nodes))
                        f.write(str(j))
                        f.write("\t")
                        f.write(str(k + nodes))
                        f.write('\n')
            nodes = nodes + new_nodes


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
    gen_DAG(min_weight=5, max_weight=10, min_height=5, max_height=6,
            percent=10, save="DAG.txt")
    # vis_graph("DAG.txt")
