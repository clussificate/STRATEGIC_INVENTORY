# -*- coding: utf-8 -*-
"""
@Created at 2020/8/14 17:07
@Author: Kurt
@file:BOMGraph.py
@Desc:
"""
from random import randint
from utils import DefinedException
import random
from pprint import pprint

random.seed(1994)


def gen_gurateed_service_times():
    return randint(1, 5)


def gen_leadtime():
    return randint(1, 20)


def gen_random_holding_cost():
    return randint(1, 20)


class BOMGraph:

    def __init__(self, filename=None):
        self.filename = filename
        self.nodes = {}
        self.supply_nodes = []
        self.demand_nodes = []

        self.gen_BOMGraph()

    def gen_BOMGraph(self):
        with open(self.filename, 'r') as f:
            for line in f:
                source, sink = line.strip().split('\t')
                if source not in self.nodes.keys():
                    res = {'sink': [sink], 'source': []}
                    self.nodes[source] = res
                else:
                    self.nodes[source]['sink'].append(sink)

                if sink not in self.nodes.keys():
                    res = {'source': [source], 'sink': []}
                    self.nodes[sink] = res
                else:
                    self.nodes[sink]['source'].append(source)

        self.update_info()

    def update_info(self):
        for node, info in self.nodes.items():
            lead_time = gen_leadtime()
            self.nodes[node]['lead_time'] = lead_time
            self.nodes[node]['holding_cost'] = 0

            if not info['source']:
                self.supply_nodes.append(node)

            if not info['sink']:
                self.demand_nodes.append(node)
                self.nodes[node]['demand_service_time'] = gen_gurateed_service_times()

        self.ascent_holding_cost()

    def topo_sort(self):
        """
        topological sorting for DAG
        """
        in_degrees = {key: len(val['source']) for key, val in self.nodes.items()}
        sort = []
        roots = self.supply_nodes[:]
        while roots:
            root = roots.pop()

            sort.append(root)
            for sink_node in self.nodes[root]['sink']:
                in_degrees[sink_node] -= 1

                if in_degrees[sink_node] == 0:
                    roots.insert(0, sink_node)
        num_non_zero = len([val for val in in_degrees.values() if val > 0])
        if num_non_zero:
            # if have non zero in-degrees, detect a circle
            raise DefinedException("please input a DAG.")
        else:
            return sort

    def ascent_holding_cost(self):
        topo_sort = self.topo_sort()
        # print(topo_sort)
        # print(self.nodes)
        for node in topo_sort:
            # print("current node {}".format(node))
            if self.nodes[node]["source"]:
                for source_node in self.nodes[node]["source"]:
                    self.nodes[node]["holding_cost"] += self.nodes[source_node]["holding_cost"]

                self.nodes[node]["holding_cost"] += randint(0, 5)  # positive additional costs
            else:
                self.nodes[node]["holding_cost"] = gen_random_holding_cost()


if __name__ == "__main__":
    bom = BOMGraph("DAG.txt")
    pprint(bom.nodes)
    print(len(bom.nodes))
    print(bom.supply_nodes)
    print(bom.demand_nodes)
