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
        self.model.setParam("OutputFlag", False)

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
            self.model.addConstr(SI + lead_time - S, GRB.GREATER_EQUAL, 0.0, name="c1_" + j)

            if not info['sink']:
                S_bar = info['demand_service_time']
                self.model.addConstr(S, GRB.LESS_EQUAL, S_bar, name="c2_" + j)

            for i in info['source']:
                S_i = self.model.getVarByName("S_" + i)
                self.model.addConstr(SI - S_i, GRB.GREATER_EQUAL, 0.0, name="c3_" + i + "_" + j)

        # add objective function
        self.obj = None
        self.optimal_value = None
        # print(self.model.getVars())

    def termination_criterion(self):

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
            print("*********************New Iteration***********************")
            self.iteration += 1
            print("Current iter: {}".format(self.iteration))
            # print(self.iteration)
            # print(self.nodes)
            # print(self.model.getConstrs())
            self.optimization()

            if self.model.status == GRB.OPTIMAL:
                # print("Solution: \n {}".format(self.model.getVars()))
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
            SI = self.model.getVarByName("SI_" + j)
            S = self.model.getVarByName("S_" + j)
            lead_time = info['lead_time']
            holding_cost = info['holding_cost']
            alpha = self.alpha[j]
            self.obj += holding_cost * alpha * (SI + lead_time - S)
        self.model.setObjective(self.obj, GRB.MINIMIZE)
        self.model.optimize()

    def update_para(self):
        for j, info in self.nodes.items():
            # print(self.model.getVarByName("S_" + i).x)
            net_replenishment_period = (
                    self.model.getVarByName("SI_" + j).x + info['lead_time'] - self.model.getVarByName("S_" + j).x)
            print("-------------------------")
            print("Current X: {}".format(net_replenishment_period))
            print("Current alpha: {}".format(self.alpha[j]))
            if self.alpha[j] * net_replenishment_period == truth_function(net_replenishment_period):
                continue
            else:
                self.alpha[j] = truth_function(net_replenishment_period) / net_replenishment_period

            print("Updated alpha: {}".format(self.alpha[j]))

    def cal_optimal_value(self):
        optimal_value = 0
        for j, info in self.nodes.items():
            net_replenishment_period = (
                    self.model.getVarByName("SI_" + j).x + info['lead_time'] - self.model.getVarByName("S_" + j).x)
            # print("current j:{}, net x:{}".format(j, net_replenishment_period))
            optimal_value += info['holding_cost'] * truth_function(round(net_replenishment_period, 3))

        return optimal_value

    def update_model(self):
        self.model.remove(self.model.getConstrs())
        self.model.remove(self.model.getVars())

        # outbound service time, S
        [self.model.addVar(lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="S_" + node)
         for node in self.nodes.keys()]

        #  inbound service time, SI
        [self.model.addVar(lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="SI_" + node)
         for node in self.nodes.keys()]

        self.model.update()
        for j, info in self.nodes.items():
            S = self.model.getVarByName("S_" + j)
            SI = self.model.getVarByName("SI_" + j)
            lead_time = info['lead_time']
            self.model.addConstr(S - SI, GRB.LESS_EQUAL, lead_time, name="c1_" + j)

            if not info['sink']:
                S_bar = info['demand_service_time']
                self.model.addConstr(S, GRB.LESS_EQUAL, S_bar, name="c2_" + j)

            for i in info['source']:
                S_i = self.model.getVarByName("S_" + i)
                self.model.addConstr(SI - S_i, GRB.GREATER_EQUAL, 0, name="c3_" + i + "_" + j)

        self.model.update()


def parse_results(instance: IterativeLP) -> None:
    with open("Lp solution.txt", "w") as f:
        for j, info in instance.nodes.items():
            SI = instance.model.getVarByName("SI_"+j)
            S = instance.model.getVarByName("S_" + j)
            print("Node: {}, SI:{}, S: {}".format(j, SI.x, S.x))
            print("Net replenishment period: {}".format(SI.x+info['lead_time']-S.x))
            f.write("{}\t{}\n".format(j, SI.x+info['lead_time']-S.x))


if __name__ == "__main__":

    Nodes = BOMGraph("DAG.txt").nodes
    start = time.time()
    ILP = IterativeLP(nodes=Nodes)
    ILP.iteration_process()
    print("Optimal value: {}".format(ILP.optimal_value))
    parse_results(ILP)
    print("Used cpu time：{}".format(time.time()-start))
