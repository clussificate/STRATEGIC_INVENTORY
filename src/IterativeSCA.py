# -*- coding: utf-8 -*-
"""
@Created at 2020/8/18 20:09
@Author: Kurt
@file:IterativeSCA.py
@Desc:
"""
from gurobipy import GRB
import gurobipy as gp
import numpy as np
from BOMGraph import BOMGraph
from utils import DefinedException
import time
import random
from collections import defaultdict

true_function = np.sqrt


def cal_coefficient(x):
    """
    :return: derivation and intercept of negative sqrt function
    """
    if x == 0:
        return -999999, 0
    derivation = - (1 / 2) * np.power(x, -1 / 2)
    intercept = - (1 / 2) * np.power(x, 1 / 2)
    return derivation, intercept


class IterativeSCA:
    def __init__(self, nodes, epsilon=0.1, **start_points):
        self.epsilon = epsilon
        self.iteration = 0
        self.nodes = nodes
        self.start_point = start_points["start_points"]
        # print("sp: {}, info: {}".format(self.start_point[0], self.start_point[1]))

        self.node_to_label = defaultdict(int)
        self.label_to_node = defaultdict(str)

        for node in nodes:
            num = len(self.node_to_label)
            self.node_to_label[node] = num
            self.label_to_node[num] = node

        self.N = len(self.node_to_label)

        if isinstance(self.start_point[0], (int, float)):
            initial_alpha, initial_beta = cal_coefficient(self.start_point[0])
            # initial parameters
            self.alpha = dict(zip(self.label_to_node.keys(), [initial_alpha] * self.N))
            self.beta = dict(zip(self.label_to_node.keys(), [initial_beta] * self.N))
        elif str.lower(self.start_point[0]) == "random":
            initial_params = [cal_coefficient(random.uniform(1, 50)) for _ in self.label_to_node]
            initial_alphas = [x[0] for x in initial_params]
            initial_betas = [x[1] for x in initial_params]
            self.alpha = dict(zip(self.label_to_node.keys(), initial_alphas))
            self.beta = dict(zip(self.label_to_node.keys(), initial_betas))
        elif str.lower(self.start_point[0]) == "warm":
            self.alpha = defaultdict(float)
            self.beta = defaultdict(float)
            for node_id in self.label_to_node:
                node = self.label_to_node[node_id]
                self.alpha[node_id] = cal_coefficient(self.start_point[1][node])[0]
                self.beta[node_id] = cal_coefficient(self.start_point[1][node])[1]
        else:
            raise DefinedException("Incorrect parameter for start_point")

        # print("initial alphas: {}".format(self.alpha))
        # print("initial betas: {}".format(self.beta))

        # initial model
        self.model = gp.Model("IterativeSCA")
        self.model.setParam("OutputFlag", False)

        # add variables
        # inbound service time
        self.SI = self.model.addVars(self.N, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='SI')
        # outbound service time
        self.S = self.model.addVars(self.N, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='S')
        # initial CA
        self.t = self.model.addVars(self.N, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='t')

        # add constraints
        # Iteration constraints
        self.c1 = self.model.addConstrs((- self.t[j] -
                                         (self.alpha[j] *
                                          (self.SI[j] + self.nodes[self.label_to_node[j]]['lead_time'] - self.S[j])
                                          + self.beta[j]) <= 0.0
                                         for j in self.label_to_node))

        # fixed constraints
        for j, info in self.nodes.items():
            node_id = self.node_to_label[j]

            self.model.addConstr(self.SI[node_id] + info["lead_time"] - self.S[node_id], GRB.GREATER_EQUAL, 0.0, )

            if not info['sink']:
                self.model.addConstr(self.S[node_id], GRB.LESS_EQUAL, info['demand_service_time'])

            for source_node in info['source']:
                source_node_id = self.node_to_label[source_node]
                self.model.addConstr(self.SI[node_id] - self.S[source_node_id], GRB.GREATER_EQUAL, 0.0)

        self.obj = None
        self.optimal_value = None
        self.update_error = None

    def iter_process(self):
        while True:
            print("************************* New Iteration ***********************************")
            self.iteration += 1
            print("Current iter: {}".format(self.iteration))
            self.optimization()

            if self.model.status == GRB.OPTIMAL:
                print("Current optimal solution of approximation function: {}".format(
                    self.model.getObjective().getValue()))
                print("Current optimal solution of true function: {}".format(self.cal_optimal_value()))
                # parse_results(self)
                # print("Solution: \n {}".format(self.model.getVars()))

                # if self.termination_criterion("always"):
                if self.termination_criterion():
                    self.optimal_value = self.cal_optimal_value()
                    break
                self.update_para()
                self.update_model()  # need to update model, because params alpha and beta are changed.
            else:
                raise DefinedException("Model infeasible.")

    def optimization(self):
        self.obj = gp.LinExpr()
        for j, info in self.nodes.items():
            node_id = self.node_to_label[j]
            holding_cost = info['holding_cost']
            t = self.t[node_id]
            self.obj += holding_cost * t
        self.model.setObjective(self.obj, GRB.MINIMIZE)
        self.model.optimize()

    def termination_criterion(self, always=None):
        if always == "always":
            return False

        flag = False
        err = 0
        for j, info in self.nodes.items():
            node_id = self.node_to_label[j]
            net_replenishment_period = self.SI[node_id].x + info['lead_time'] - self.S[node_id].x
            appr_function = self.t[node_id].x
            err += info['holding_cost'] * abs(appr_function - true_function(round(net_replenishment_period, 3)))
        print("Current error: {} of iteration: {}".format(err, self.iteration))

        if err <= self.epsilon / self.N:
            flag = True
        if self.update_error == err:
            flag = True
            print("No reduced error, iteration ends.")
        self.update_error = err
        return flag

    def update_para(self):
        for j, info in self.nodes.items():
            node_id = self.node_to_label[j]
            net_replenishment_period = round(self.SI[node_id].x + info['lead_time'] - self.S[node_id].x, 3)

            # print("---------------------------------------")
            # print("Current node: {}".format(j))
            # print("Current Net X:{}".format(net_replenishment_period))
            # print("previous alpha: {}".format(self.alpha[node_id]))
            # print("previous beta: {}".format(self.beta[node_id]))

            if abs(self.alpha[node_id] * net_replenishment_period + self.beta[node_id]) - true_function(
                    net_replenishment_period) < 0.01:
                continue
            else:
                self.alpha[node_id], self.beta[node_id] = cal_coefficient(net_replenishment_period)

            # print("updated alpha: {}".format(self.alpha[node_id]))
            # print("update beta: {}".format(self.beta[node_id]))

    def update_model(self):
        # Iteration constraints
        self.model.remove(self.c1)

        # Iteration constraints
        self.c1 = self.model.addConstrs((- self.t[j] -
                                         (self.alpha[j] *
                                          (self.SI[j] + self.nodes[self.label_to_node[j]]['lead_time'] - self.S[j])
                                          + self.beta[j]) <= 0.0
                                         for j in self.label_to_node))
        self.model.update()

    def cal_optimal_value(self):
        optimal_value = 0
        for j, info in self.nodes.items():
            node_id = self.node_to_label[j]
            net_replenishment_period = round(self.SI[node_id].x + info['lead_time'] - self.S[node_id].x, 3)
            # print(net_replenishment_period)
            # print("current j:{}, net x:{}".format(j, net_replenishment_period))
            optimal_value += info['holding_cost'] * true_function(net_replenishment_period)

        return optimal_value


def parse_results(instance: IterativeSCA) -> None:
    for j, info in instance.nodes.items():
        node_id = instance.node_to_label[j]
        SI = instance.SI[node_id]
        S = instance.S[node_id]
        # print("Node: {}, SI:{:.3f}, S: {:.3f}".format(j, SI.x, S.x))
        # print("Net replenishment period: {:.3f}".format(SI.x + info['lead_time'] - S.x))


if __name__ == "__main__":
    Nodes = BOMGraph("DAG.txt").nodes

    start_point = defaultdict(float)
    with open("lp solution.txt", "r") as f:
        for line in f.readlines():
            node, net_replenishment = line.strip().split("\t")
            start_point[node] = float(net_replenishment)

    IterativeSCA = IterativeSCA(nodes=Nodes, start_points=("random", start_point))
    # IterativeSCA = IterativeSCA(nodes=Nodes, start_points=(1, None))

    start = time.time()
    IterativeSCA.iter_process()
    IterativeSCA.cal_optimal_value()
    parse_results(IterativeSCA)
    print("optimal value: {}".format(IterativeSCA.optimal_value))
    print("Used cpu time：{}".format(time.time() - start))
