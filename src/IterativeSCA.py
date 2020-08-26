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

true_function = np.sqrt
M = 9999


def cal_coefficient(x):
    """
    :return: derivation and intercept of negative sqrt function
    """
    if x == 0:
        return 0, 0
    derivation = - (1 / 2) * np.power(x, -1 / 2)
    intercept = - (1 / 2) * np.power(x, 1 / 2)
    return derivation, intercept


class IterativeSCA:
    def __init__(self, nodes, epsilon=0.1, start_point=None):
        self.epsilon = epsilon
        self.iteration = 0
        self.nodes = {int(k): v for k, v in nodes.items()}
        self.N = len(self.nodes)

        self.start_point = start_point
        if isinstance(self.start_point, (int, float)):
            initial_alpha, initial_beta = cal_coefficient(self.start_point)

            # initial parameters
            self.alpha = dict(zip(self.nodes.keys(), [initial_alpha] * self.N))
            self.beta = dict(zip(self.nodes.keys(), [initial_beta] * self.N))

        elif str.lower(self.start_point) == "random":
            initial_params = [cal_coefficient(random.uniform(1, 50)) for _ in self.nodes.keys()]
            initial_alphas = [x[0] for x in initial_params]
            initial_betas = [x[1] for x in initial_params]
            self.alpha = dict(zip(self.nodes.keys(), initial_alphas))
            self.beta = dict(zip(self.nodes.keys(), initial_betas))
        else:
            raise DefinedException("Incorrect parameter for start_point")

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
                                          + self.beta[j]) <= 0.0
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
                print("Current optimal solution of approximation function: {}".format(self.model.getObjective().getValue()))
                print("Current optimal solution of true function: {}".format(self.cal_optimal_value()))
                # parse_results(self)
                # print("Solution: \n {}".format(self.model.getVars()))
                if self.termination_criterion():
                    self.optimal_value = self.cal_optimal_value()
                    break
                self.update_para()
                self.update_model()  # need to update model, because params alpha and beta are changed.
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
            print("Current node: {}".format(j))
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
            net_replenishment_period = round(self.SI[j].x + info['lead_time'] - self.S[j].x, 3)
            # print(net_replenishment_period)
            # print("current j:{}, net x:{}".format(j, net_replenishment_period))
            optimal_value += info['holding_cost'] * true_function(net_replenishment_period)

        return optimal_value

    def update_model(self):
        # Iteration constraints
        self.model.remove(self.c1)

        # Iteration constraints
        self.c1 = self.model.addConstrs((- self.t[j] -
                                         (self.alpha[j] *
                                          (self.SI[j] + self.nodes[j]['lead_time'] - self.S[j])
                                          + self.beta[j]) <= 0.0
                                         for j in self.nodes.keys()))
        self.model.update()


def parse_results(instance: IterativeSCA) -> None:
    for j, info in instance.nodes.items():
        SI = instance.SI[j]
        S = instance.S[j]
        print("Node: {}, SI:{}, S: {}".format(j, SI.x, S.x))
        print("Net replenishment period: {}".format(SI.x + info['lead_time'] - S.x))


if __name__ == "__main__":
    Nodes = BOMGraph("DAG.txt").nodes

    IterativeSCA = IterativeSCA(nodes=Nodes, start_point=1)
    start = time.time()
    IterativeSCA.iter_process()
    IterativeSCA.cal_optimal_value()
    print("optimal value: {}".format(IterativeSCA.optimal_value))
    parse_results(IterativeSCA)
    print("Used cpu time：{}".format(time.time() - start))

