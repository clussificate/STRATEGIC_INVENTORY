# -*- coding: utf-8 -*-
"""
@Created at 2020/8/15 10:11
@Author: Kurt
@file:IterativeMIP.py
@Desc:
"""
from gurobipy import GRB
import gurobipy as gp
import numpy as np
from BOMGraph import BOMGraph
from utils import DefinedException
from IterativeLP import IterativeLP
import time

truth_function = np.sqrt
M = 9999


def cal_coefficient(x):
    """
    :return: derivation and intercept of sqrt function
    """
    if x == 0:
        return 0, 0
    derivation = (1/2) * np.power(x, -1/2)
    intercept = 1/2 * np.power(x, 1/2)
    return derivation, intercept


class IterativeMIP(IterativeLP):
    def __init__(self, nodes, epsilon=0.1):
        super().__init__(nodes, epsilon)
        self.alpha = dict(zip(self.nodes.keys(), [1 / 2] * self.N))
        self.beta = dict(zip(self.nodes.keys(), [1 / 2] * self.N))
        self.y = self.model.addVars(self.N, vtype=GRB.BINARY, name='y')

        self.model.setParam("OutputFlag", False)
        self.model.update()
        for j, info in self.nodes.items():
            S = self.model.getVarByName("S_" + j)
            SI = self.model.getVarByName("SI_" + j)
            lead_time = info['lead_time']
            self.model.addConstr((SI+lead_time-S), GRB.LESS_EQUAL, M*self.y[int(j)])
        self.update_error = None

    def termination_criterion(self):

        flag = False
        error = 0
        # print("current alpha: {}".format(self.alpha))
        # print("current beta: {}".format(self.beta))
        for j, info in self.nodes.items():
            net_replenishment_period = round((
                    self.model.getVarByName("SI_" + j).x + info['lead_time'] - self.model.getVarByName("S_" + j).x),3)

            y = self.y[int(j)].x

            error += info['holding_cost'] * abs(
                (self.alpha[j] * net_replenishment_period + y * self.beta[j]) - truth_function(net_replenishment_period))
            # print("Cum error: {}".format(error))
        print("Current error: {} of iteration: {}".format(error, self.iteration))
        if error <= self.epsilon / self.N:
            flag = True
        if self.update_error == error:
            flag = True
            print("No reduced error, iteration ends.")
        self.update_error = error
        return flag


    def iteration_process(self):
        while True:
            print("************************* New Iteration ***********************************")
            self.iteration += 1
            self.optimization()
            print("Current iter: {}".format(self.iteration))

            if self.model.status == GRB.OPTIMAL:
                print("Current optimal solution of true function: {}".format(self.cal_optimal_value()))
                # print("Current appr obj value: {}".format(self.model.objVal))
                # print("Current true obj value: {}".format(self.cal_optimal_value()))
                # print("Solution: \n {}".format(self.model.getVars()))
                if self.termination_criterion():
                    self.optimal_value = self.cal_optimal_value()
                    break
                self.update_para()
            else:
                raise DefinedException("No solutions.")

    def optimization(self):
        self.obj = gp.LinExpr()
        for j, info in self.nodes.items():
            SI = self.model.getVarByName("SI_" + j)
            S = self.model.getVarByName("S_" + j)
            y = self.y[int(j)]
            lead_time = info['lead_time']
            holding_cost = info['holding_cost']
            alpha = self.alpha[j]
            beta = self.beta[j]

            self.obj += holding_cost * (alpha * (SI + lead_time - S) + y * beta)
        self.model.setObjective(self.obj, GRB.MINIMIZE)

        self.model.optimize()

    def update_para(self):
        for j, info in self.nodes.items():
            net_replenishment_period = round((
                    self.model.getVarByName("SI_" + j).x + info['lead_time'] - self.model.getVarByName("S_" + j).x), 3)
            y = self.y[int(j)].x
            print("---------------------------------------")
            print("Current Net X:{}".format(net_replenishment_period))
            print("Current y: {}".format(y))
            print("previous alpha: {}".format(self.alpha[j]))
            print("previous beta: {}".format(self.beta[j]))

            if abs(self.alpha[j] * net_replenishment_period + y * self.beta[j]) - truth_function(net_replenishment_period) < 0.01:
                continue
            else:
                self.alpha[j], self.beta[j] = cal_coefficient(net_replenishment_period)

            print("updated alpha: {}".format(self.alpha[j]))
            print("update beta: {}".format(self.beta[j]))

    def cal_optimal_value(self):
        optimal_value = 0
        for i, info in self.nodes.items():
            net_replenishment_period = round((
                    self.model.getVarByName("SI_" + i).x + info['lead_time'] - self.model.getVarByName("S_" + i).x),3)

            optimal_value += info['lead_time'] * truth_function(net_replenishment_period)

        return optimal_value


def parse_result(instance:IterativeMIP) -> None:
    for j, info in instance.nodes.items():
        SI = instance.model.getVarByName("SI_"+j)
        S = instance.model.getVarByName("S_" + j)
        print("Node: {}, SI:{}, S: {}".format(j, SI.x, S.x))


if __name__ == "__main__":

    Nodes = BOMGraph("DAG.txt").nodes
    start = time.time()
    # print(Nodes)
    IMIP = IterativeMIP(nodes=Nodes)
    IMIP.iteration_process()
    print("Optimal value: {}".format(IMIP.optimal_value))
    parse_result(IMIP)
    print("Used cpu time：{}".format(time.time() - start))
