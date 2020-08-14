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

truth_function = np.sqrt


class IterativeLP:
    def __init__(self, nodes, epsilon=0.01):
        self.epsilon = epsilon
        self.iteration = 0
        self.nodes = nodes
        self.N = len(self.nodes)
        # initial alpha for each node
        self.alpha = dict(zip(self.nodes.keys(), [1] * self.N))

        # initial model
        self.model = gp.Model("IterativeLP")
        # outbound service time, S
        [self.model.addVar(lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="S_" + node)
         for node in self.nodes.keys()]
        #  inbound service time, SI
        [self.model.addVar(lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="SI_" + node)
         for node in self.nodes.keys()]

        self.model.update()

        # add constraints
        for j, info in self.nodes.items():
            S = self.model.getVarByName("S_" + j)
            SI = self.model.getVarByName("SI_" + j)
            lead_time = info['lead_time']
            self.model.addConstr(S - SI, GRB.LESS_EQUAL, lead_time)

            if not info['sink']:
                S_bar = info['demand_service_time']
                self.model.addConstr(S, GRB.LESS_EQUAL, S_bar)

            for i in info['source']:
                S_i = self.model.getVarByName("S_" + i)
                self.model.addConstr(SI - S_i, GRB.GREATER_EQUAL, 0)

        # add objective function
        self.obj = None
        self.optimal_value = None
        # print(self.model.getVars())
        self.iteration_process()

    def termination_criterion(self):
        # if self.iteration == 0:
        #     return False

        flag = False
        error = 0
        for i, info in self.nodes.items():
            net_replenishment_period = (
                    self.model.getVarByName("SI_" + i).x + info['lead_time'] - self.model.getVarByName("S_" + i).x)

            # print("Net replenishment: {} ".format(net_replenishment_period))
            # print("alpha:{}".format(self.alpha[i]))
            # print("holding_cost:{}".format(info['holding_cost']))
            # print("current error:{}".format(info['holding_cost'] * abs(
            #     self.alpha[i] * net_replenishment_period - truth_function(net_replenishment_period))))

            error += info['holding_cost'] * abs(
                self.alpha[i] * net_replenishment_period - truth_function(net_replenishment_period))
            # print("Cum error: {}".format(error))
        print("Current error: {} of iteration: {}".format(error, self.iteration))
        if error <= self.epsilon / self.N:
            flag = True

        return flag

    def iteration_process(self):
        while True:
            self.iteration += 1
            print(self.iteration)
            self.optimization()

            if self.model.status == GRB.OPTIMAL:
                print("Iteration: {}".format(self.iteration))
                print("Solution: \n {}".format(self.model.getVars()))
                if self.termination_criterion():
                    self.optimal_value = self.cal_optimal_value()
                    break
                self.update_alpha()
            else:
                raise DefinedException("No solutions.")

    def optimization(self):
        self.obj = gp.LinExpr()
        for i, info in self.nodes.items():
            SI = self.model.getVarByName("SI_" + i)
            S = self.model.getVarByName("S_" + i)
            lead_time = info['lead_time']
            holding_cost = info['holding_cost']
            alpha = self.alpha[i]
            self.obj += holding_cost * alpha * (SI + lead_time - S)
        self.model.setObjective(self.obj, GRB.MINIMIZE)

        self.model.optimize()

    def update_alpha(self):
        for i, info in self.nodes.items():
            # print(self.model.getVarByName("S_" + i).x)
            net_replenishment_period = (
                    self.model.getVarByName("SI_" + i).x + info['lead_time'] - self.model.getVarByName("S_" + i).x)
            if self.alpha[i] * net_replenishment_period == truth_function(net_replenishment_period):
                continue
            else:
                self.alpha[i] = truth_function(net_replenishment_period) / net_replenishment_period

    def cal_optimal_value(self):
        optimal_value = 0
        for i, info in self.nodes.items():
            net_replenishment_period = (
                    self.model.getVarByName("SI_" + i).x + info['lead_time'] - self.model.getVarByName("S_" + i).x)
            optimal_value += truth_function(net_replenishment_period)

        return optimal_value


if __name__ == "__main__":
    print("time:")
    Nodes = BOMGraph("DAG.txt").nodes
    # print(Nodes)
    ILP = IterativeLP(nodes=Nodes)
    print("Optimal value: {}".format(ILP.optimal_value))
