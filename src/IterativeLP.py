# -*- coding: utf-8 -*-
"""
@Created at 2020/8/14 16:57
@Author: Kurt
@file:IterativeLP.py
@Desc:
"""
from gurobipy import GRB
import gurobipy as gp
import numpy as np
from BOMGraph import BOMGraph
from utils import DefinedException
import time
from collections import defaultdict

true_function = np.sqrt


class IterativeLP:
    def __init__(self, nodes, epsilon=0.01):
        self.epsilon = epsilon
        self.iteration = 0
        self.nodes = nodes

        self.node_to_label = defaultdict(int)
        self.label_to_node = defaultdict(str)

        for node in nodes:
            num = len(self.node_to_label)
            self.node_to_label[node] = num
            self.label_to_node[num] = node

        self.N = len(self.node_to_label)

        # initial alpha for each node
        self.alpha = dict(zip(self.label_to_node.keys(), [1] * self.N))

        # initial model
        self.model = gp.Model("IterativeLP")
        self.model.setParam("OutputFlag", False)

        # add variables
        # inbound service time
        self.SI = self.model.addVars(self.N, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='SI')
        # outbound service time
        self.S = self.model.addVars(self.N, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='S')

        self.model.update()

        # add constraints
        for j, info in self.nodes.items():
            node_id = self.node_to_label[j]

            self.model.addConstr(self.SI[node_id] + info["lead_time"] - self.S[node_id], GRB.GREATER_EQUAL, 0.0,)

            if not info['sink']:
                self.model.addConstr(self.S[node_id], GRB.LESS_EQUAL, info['demand_service_time'])

            for source_node in info['source']:
                source_node_id = self.node_to_label[source_node]
                self.model.addConstr(self.SI[node_id] - self.S[source_node_id], GRB.GREATER_EQUAL, 0.0)

        # add objective function
        self.obj = None
        self.optimal_value = None
        self.update_error = None
        # print(self.model.getVars())

    def termination_criterion(self):

        flag = False
        err = 0
        for j, info in self.nodes.items():
            node_id = self.node_to_label[j]
            net_replenishment_period = self.SI[node_id].x + info['lead_time'] - self.S[node_id].x

            # print("Net replenishment: {} ".format(net_replenishment_period))
            # print("alpha:{}".format(self.alpha[i]))
            # print("holding_cost:{}".format(info['holding_cost']))
            # print("current error:{}".format(info['holding_cost'] * abs(
            #     self.alpha[i] * net_replenishment_period - truth_function(net_replenishment_period))))

            err += info['holding_cost'] * abs(
                self.alpha[node_id] * net_replenishment_period - true_function(net_replenishment_period))
            # print("Cum error: {}".format(error))
        print("Current error: {} of iteration: {}".format(err, self.iteration))

        if err <= self.epsilon / self.N:
            flag = True
        if self.update_error == err:
            flag = True
            print("No reduced error, iteration ends.")
        self.update_error = err
        return flag

    def iteration_process(self):
        while True:
            print("*********************New Iteration***********************")
            self.iteration += 1
            print("Current iter: {}".format(self.iteration))
            # print(self.iteration)
            # print(self.nodes)
            # print(self.model.getConstrs())
            self.optimization()

            if self.model.status == GRB.OPTIMAL:
                # print("Solution: \n {}".format(self.model.getVars()))
                print("Current optimal solution of approximation function: {}".format(
                    self.model.getObjective().getValue()))
                print("Current optimal solution of true function: {}".format(self.cal_optimal_value()))
                # parse_results(self)
                if self.termination_criterion():
                    self.optimal_value = self.cal_optimal_value()
                    break
                self.update_para()
                # self.update_model()   # all static constraints, we don't need to update model
            else:
                raise DefinedException("No solutions.")

    def optimization(self):
        self.obj = gp.LinExpr()
        for j, info in self.nodes.items():
            node_id = self.node_to_label[j]
            SI = self.SI[node_id]
            S = self.S[node_id]
            lead_time = info['lead_time']
            holding_cost = info['holding_cost']
            alpha = self.alpha[node_id]
            self.obj += holding_cost * alpha * (SI + lead_time - S)
        self.model.setObjective(self.obj, GRB.MINIMIZE)
        self.model.optimize()

    def update_para(self):
        for j, info in self.nodes.items():
            node_id = self.node_to_label[j]
            net_replenishment_period = self.SI[node_id].x + info['lead_time'] - self.S[node_id].x
            # print("-------------------------")
            # print("Current X: {}".format(net_replenishment_period))
            # print("Current alpha: {}".format(self.alpha[node_id]))
            if self.alpha[node_id] * net_replenishment_period == true_function(net_replenishment_period):
                continue
            else:
                self.alpha[node_id] = true_function(net_replenishment_period) / net_replenishment_period

            # print("Updated alpha: {}".format(self.alpha[j]))

    def cal_optimal_value(self):
        optimal_value = 0
        for j, info in self.nodes.items():
            node_id = self.node_to_label[j]
            net_replenishment_period = round(self.SI[node_id].x + info['lead_time'] - self.S[node_id].x, 3)
            # print(net_replenishment_period)
            # print("current j:{}, net x:{}".format(j, net_replenishment_period))
            optimal_value += info['holding_cost'] * true_function(net_replenishment_period)

        return optimal_value

    def update_model(self):
        self.model.remove(self.model.getConstrs())
        self.model.remove(self.model.getVars())

        # inbound service time
        self.SI = self.model.addVars(self.N, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='SI')
        # outbound service time
        self.S = self.model.addVars(self.N, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='S')

        # add constraints
        for j, info in self.nodes.items():
            node_id = self.node_to_label[j]

            self.model.addConstr(self.SI[node_id] + info["lead_time"] - self.S[node_id], GRB.GREATER_EQUAL, 0.0, )

            if not info['sink']:
                self.model.addConstr(self.S[node_id], GRB.LESS_EQUAL, info['demand_service_time'])

            for source_node in info['source']:
                source_node_id = self.node_to_label[source_node]
                self.model.addConstr(self.SI[node_id] - self.S[source_node_id], GRB.GREATER_EQUAL, 0.0)

        self.model.update()


def parse_results(instance: IterativeLP) -> None:
    with open("lp solution.txt", "w") as f:
        for j, info in instance.nodes.items():
            node_id = instance.node_to_label[j]
            SI = instance.SI[node_id]
            S = instance.S[node_id]
            # print("Node: {}, SI:{}, S: {}".format(j, SI.x, S.x))
            # print("Net replenishment period: {}".format(SI.x+info['lead_time']-S.x))
            f.write("{}\t{}\n".format(j, SI.x + info['lead_time'] - S.x))


if __name__ == "__main__":
    Nodes = BOMGraph("DAG.txt").nodes

    start = time.time()
    ILP = IterativeLP(nodes=Nodes)
    ILP.iteration_process()
    parse_results(ILP)
    print("Optimal value: {}".format(ILP.optimal_value))
    print("Used cpu time：{}".format(time.time() - start))
