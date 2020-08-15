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

truth_function = np.sqrt
M = 9999


def cal_coefficient(x):
    """
    :return: derivation and intercept of sqrt function
    """
    if x ==0:
        return 0, 0
    derivation = (1/2) * np.power(x, -1/2)
    intercept = 1/2 * np.power(x, 1/2)
    return derivation, intercept


class IterativeMIP(IterativeLP):
    def __init__(self, nodes, epsilon=0.1):
        super().__init__(nodes, epsilon)
        self.alpha = dict(zip(self.nodes.keys(), [1 / 2] * self.N))
        self.beta = dict(zip(self.nodes.keys(), [1 / 2] * self.N))
        self.model.addVar(vtype=GRB.BINARY, name='y')
        self.model.setParam("OutputFlag", False)
        self.model.update()
        for j, info in self.nodes.items():
            S = self.model.getVarByName("S_" + j)
            SI = self.model.getVarByName("SI_" + j)
            lead_time = info['lead_time']
            self.model.addConstr((SI+lead_time-S), GRB.LESS_EQUAL, M*self.model.getVarByName("y"))

    def termination_criterion(self):

        flag = False
        error = 0
        # print("current alpha: {}".format(self.alpha))
        # print("current beta: {}".format(self.beta))
        for i, info in self.nodes.items():
            net_replenishment_period = (
                    self.model.getVarByName("SI_" + i).x + info['lead_time'] - self.model.getVarByName("S_" + i).x)

            y = self.model.getVarByName("y").x

            error += info['holding_cost'] * abs(
                (self.alpha[i] * net_replenishment_period + y * self.beta[i]) - truth_function(net_replenishment_period))
            # print("Cum error: {}".format(error))
        print("Current error: {} of iteration: {}".format(error, self.iteration))
        if error <= self.epsilon / self.N:
            flag = True

        return flag

    def iteration_process(self):
        k=0
        while True and k<3:
            self.iteration += 1
            k+=1
            self.optimization()

            if self.model.status == GRB.OPTIMAL:
                print("Current appr obj value: {}".format(self.model.objVal))
                print("Current true obj value: {}".format(self.cal_optimal_value()))
                # print("Solution: \n {}".format(self.model.getVars()))
                if self.termination_criterion():
                    self.optimal_value = self.cal_optimal_value()
                    break
                self.update_para()
            else:
                raise DefinedException("No solutions.")

    def optimization(self):
        self.obj = gp.LinExpr()
        for i, info in self.nodes.items():
            SI = self.model.getVarByName("SI_" + i)
            S = self.model.getVarByName("S_" + i)
            y = self.model.getVarByName("y")
            lead_time = info['lead_time']
            holding_cost = info['holding_cost']
            alpha = self.alpha[i]
            beta = self.beta[i]

            if self.iteration > 1:
                print("y: {}".format(y))
                print("X: {}".format(SI.x+lead_time-S.x))
            self.obj += holding_cost * (alpha * (SI + lead_time - S) + y * beta)
        self.model.setObjective(self.obj, GRB.MINIMIZE)

        self.model.optimize()

    def update_para(self):
        for i, info in self.nodes.items():
            # print(self.model.getVarByName("S_" + i).x)
            net_replenishment_period = (
                    self.model.getVarByName("SI_" + i).x + info['lead_time'] - self.model.getVarByName("S_" + i).x)
            y = self.model.getVarByName("y").x

            if (self.alpha[i] * net_replenishment_period + y * self.beta[i]) == truth_function(net_replenishment_period):
                continue
            else:
                self.alpha[i], self.beta[i] = cal_coefficient(net_replenishment_period)

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
    IMIP = IterativeMIP(nodes=Nodes)
    IMIP.iteration_process()
    print("Optimal value: {}".format(IMIP.optimal_value))