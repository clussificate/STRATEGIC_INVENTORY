# -*- coding: utf-8 -*-
"""
@Created at 2020/8/14 17:07
@Author: Kurt
@file:BOMGraph.py
@Desc:
"""
from random import randint
from utils import DefinedException


def gen_gurateed_service_times():
    return randint(1, 10)


def gen_leadtime():
    return randint(1, 20)


def gen_holding_cost(node, mode):
    if str.lower(mode) == "random" or "r":
        return randint(1, 20)
    elif str.lower(mode) == 'increase' or "i":
        "TODO"
        pass
    else:
        raise DefinedException("Incorrect mode")


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
            holding_cost = gen_holding_cost(node, mode='r')
            self.nodes[node]['holding_cost'] = holding_cost

            if not info['source']:
                self.supply_nodes.append(node)

            if not info['sink']:
                self.demand_nodes.append(node)
                self.nodes[node]['demand_service_time'] = gen_gurateed_service_times()


if __name__ == "__main__":
    bom = BOMGraph("DAG.txt")
    print(bom.nodes)
    print(len(bom.nodes))
