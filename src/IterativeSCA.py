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
from IterativeLP import IterativeLP
import time

true_function = np.sqrt
M = 9999


def cal_coefficient(x):
    """
    :return: derivation and intercept of sqrt function
    """
    if x == 0:
        return 0, 0
    derivation = - (1 / 2) * np.power(x, -1 / 2)
    intercept = - (1 / 2) * np.power(x, 1 / 2)
    return derivation, intercept


class IterativeSCA:
    def __init__(self, nodes, epsilon=0.1):
        self.epsilon = epsilon
        self.iteration = 0
        self.nodes = {int(k): v for k, v in nodes.items()}
        self.N = len(self.nodes)

        # initial parameters
        self.alpha = dict(zip(self.nodes.keys(), [-1 / 2] * self.N))
        self.beta = dict(zip(self.nodes.keys(), [-1 / 2] * self.N))

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

        self.model.update()
        # add constraints

        # Iteration constraints
        self.c1 = self.model.addConstrs((- self.t[j] -
                                         (self.alpha[j] *
                                          (self.SI[j] + self.nodes[j]['lead_time'] - self.S[j])
                                          + self.beta[j]) <= 0
                                         for j in self.nodes.keys()))

        # unchanged constraints
        for j, info in self.nodes.items():

            self.model.addConstr(self.SI[j] + info["lead_time"] - self.S[j], GRB.GREATER_EQUAL, 0.0,
                                 name="r[1]_" + str(j))

            if not info['sink']:
                self.model.addConstr(self.S[j], GRB.LESS_EQUAL, info['demand_service_time'], name="r[2]_" + str(j))

            for i in info['source']:
                i = int(i)
                self.model.addConstr(self.SI[j] - self.S[i], GRB.GREATER_EQUAL, 0.0, "r[3]_" + str(j) + "_" + str(i))

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
            holding_cost = info['holding_cost']
            t = self.t[j]
            self.obj += holding_cost * t
        self.model.setObjective(self.obj, GRB.MINIMIZE)

        self.model.optimize()

    def termination_criterion(self):

        flag = False
        err = 0
        for j, info in self.nodes.items():
            net_replenishment_period = self.SI[j].x + info['lead_time'] - self.S[j].x
            appr_function = self.alpha[j] * self.t[j].x + self.beta[j]
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
            net_replenishment_period = round(self.SI[j].x + info['lead_time'] - self.S[j].x, 3)

            print("---------------------------------------")
            print("Current Net X:{}".format(net_replenishment_period))
            print("previous alpha: {}".format(self.alpha[j]))
            print("previous beta: {}".format(self.beta[j]))

            if abs(self.alpha[j] * net_replenishment_period + self.beta[j]) - true_function(
                    net_replenishment_period) < 0.01:
                continue
            else:
                self.alpha[j], self.beta[j] = cal_coefficient(net_replenishment_period)

            print("updated alpha: {}".format(self.alpha[j]))
            print("update beta: {}".format(self.beta[j]))

    def cal_optimal_value(self):
        optimal_value = 0
        for j, info in self.nodes.items():
            net_replenishment_period = round(self.SI[j].x + info['lead_time'] - self.S[j].x,  3)
            # print(net_replenishment_period)
            # print("current j:{}, net x:{}".format(j, net_replenishment_period))
            optimal_value += info['lead_time']*true_function(net_replenishment_period)

        return optimal_value


def parse_results(instance: IterativeSCA) -> None:
    for j, info in instance.nodes.items():
        SI = instance.SI[j]
        S = instance.S[j]
        print("Node: {}, SI:{}, S: {}".format(j, SI.x, S.x))


if __name__ == "__main__":
    Nodes = BOMGraph("DAG.txt").nodes
    # print(Nodes)
    IterativeSCA = IterativeSCA(nodes=Nodes)
    start = time.time()
    IterativeSCA.iter_process()
    IterativeSCA.cal_optimal_value()
    print("optimal value: {}".format(IterativeSCA.optimal_value))
    parse_results(IterativeSCA)
    print("Used cpu time：{}".format(time.time() - start))
    # print(IterativePW.model.getConstrs())


