# -*- coding: utf-8 -*-
"""
@Created at 2020/8/15 18:04
@Author: Kurt
@file:IterativePW.py
@Desc:
"""
import numpy as np
import gurobipy as gp
from gurobipy import GRB
from BOMGraph import BOMGraph
from utils import DefinedException

true_function = np.sqrt


class IterativePW:
    def __init__(self, nodes, epsilon=0.01):
        self.epsilon = epsilon
        self.iteration = 0
        self.nodes = {int(k): v for k, v in nodes.items()}
        self.N = len(self.nodes)

        # initial parameters
        self.R = [2] * self.N
        self.alpha = dict(zip(self.nodes.keys(), [[true_function(9999) / 9999, 0]] * self.N))
        self.beta = dict(zip(self.nodes.keys(), [[0, true_function(9999)]] * self.N))
        self.M = dict(zip(self.nodes.keys(), [[0, 9999, 99999999]] * self.N))

        # initial model
        self.model = gp.Model("IterativeLP")
        self.model.setParam("OutputFlag", True)
        # inbound service time
        self.SI = self.model.addVars(self.N, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='SI')
        # outbound service time
        self.S = self.model.addVars(self.N, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='SI')

        # piece-wide function variables
        ind_variables = [(n, t) for n, r in zip(range(self.N), self.R) for t in range(r)]

        self.z = self.model.addVars(ind_variables, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='z')
        self.u = self.model.addVars(ind_variables, vtype=GRB.BINARY, name='u')

        self.model.update()

        # add constraints
        # notice: for j,info in self.nodes.items() is not allowed.
        # z.sum(j, "*") can be replaced by sum([z[j,r] for r in self.R])
        self.pw_c0 = self.model.addConstrs(
            (self.SI[j] + self.nodes[j]["lead_time"] - self.S[j] == self.z.sum(j, "*")
             for j in self.nodes.keys()))

        self.pw_c1 = self.model.addConstrs(
            (sum([self.u[j, r] for r in range(self.R[j])]) <= 1
             for j in self.nodes.keys()))

        self.pw_c2 = self.model.addConstrs(
            (self.z[j, r] <= self.M[j][r + 1] * self.u[j, r]
             for j in self.nodes.keys() for r in range(self.R[j])))

        self.pw_c3 = self.model.addConstrs(
            (self.z[j, r] >= self.M[j][r] * self.u[j, r]
             for j in self.nodes.keys() for r in range(self.R[j])))

        for j, info in self.nodes.items():

            self.model.addConstr(self.SI[j] + info["lead_time"] - self.S[j], GRB.GREATER_EQUAL, 0.0)

            if not info['sink']:
                self.model.addConstr(self.S[j], GRB.LESS_EQUAL, info['demand_service_time'])

            for i in info['source']:
                i = int(i)
                self.model.addConstr(self.SI[j] - self.S[i], GRB.GREATER_EQUAL, 0.0)

        self.model.update()
        self.obj = None
        self.optimal_value = None

    def iter_process(self):
        while True:
            self.iteration += 1
            self.optimization()
            print("Current iter: {}".format(self.iteration))

            if self.model.status == GRB.OPTIMAL:
                # print("Solution: \n {}".format(self.model.getVars()))
                if self.termination_criterion():
                    self.optimal_value = self.cal_optimal_value()
                    break
                self.update_para()
                self.update_model()

            else:
                raise DefinedException("No solutions.")

    def optimization(self):
        """ Objective functions are updated here"""
        self.obj = gp.LinExpr()
        for j, info in self.nodes.items():
            holding_cost = info['lead_time']
            val = 0
            for r in range(self.R[j]):
                val += self.u[j, r] * self.beta[j][r] + self.z[j, r] * self.alpha[j][r]
            self.obj += holding_cost * val
        self.model.setObjective(self.obj, GRB.MINIMIZE)
        self.model.optimize()
        # print(self.model.getVars())

    def termination_criterion(self):
        flag = False
        error = 0
        for j, info in self.nodes.items():
            net_replenishment_period = self.SI[j].x + info['lead_time'] - self.S[j].x
            appr_function = sum(
                [self.beta[j][r] * self.u[j, r].x + self.alpha[j][r] * self.z[j, r].x for r in range(self.R[j])])
            error += info['holding_cost'] * abs(appr_function - true_function(net_replenishment_period))
        print("Current error: {} of iteration: {}".format(error, self.iteration))
        if error <= self.epsilon / self.N:
            flag = True
        return flag

    def cal_optimal_value(self):
        optimal_value = 0
        for j, info in self.nodes.items():
            net_replenishment_period = self.SI[j].x + info['lead_time'] - self.S[j].x
            optimal_value = true_function(net_replenishment_period)

        return optimal_value

    def update_model(self):
        self.model.remove([self.z, self.u])
        self.update_variables()

        self.model.remove([self.pw_c0, self.pw_c1, self.pw_c2, self.pw_c3])
        self.update_constraints()

        self.model.update()

    def update_para(self):
        for j, info in self.nodes.items():
            # print("current j:{} of R{}".format(j, self.R[j]))
            # print([self.u[j, r].x for r in range(self.R[j])])
            interval_index = [r for r in range(self.R[j]) if self.u[j, r].x == 1][0]
            net_replenishment_period = self.SI[j].x + info['lead_time'] - self.S[j].x

            if abs(net_replenishment_period - self.M[j][interval_index]) <= 0.01:
                continue

            self.R[j] += 1

            new_alpha = (true_function(net_replenishment_period) - true_function(self.M[j][interval_index])) / \
                        (net_replenishment_period - self.M[j][interval_index])
            new_beta = true_function(net_replenishment_period) - new_alpha * net_replenishment_period

            new_alpha_ = (true_function(self.M[j][interval_index + 1]) - true_function(net_replenishment_period)) / \
                         (self.M[j][interval_index + 1] - net_replenishment_period)
            new_beta_ = true_function(net_replenishment_period) - new_alpha_ * net_replenishment_period

            self.alpha[j][interval_index] = new_alpha
            self.beta[j][interval_index] = new_beta

            self.alpha[j].insert(interval_index + 1, new_alpha_)
            self.beta[j].insert(interval_index + 1, new_beta_)

            self.M[j].insert(interval_index+1, net_replenishment_period)

    def update_constraints(self):
        self.pw_c0 = self.model.addConstrs(
            (self.SI[j] + self.nodes[j]["lead_time"] - self.S[j] == sum([self.z[j, r] for r in range(self.R[j])])
             for j in self.nodes.keys()))

        self.pw_c1 = self.model.addConstrs(
            (sum([self.u[j, r] for r in range(self.R[j])]) <= 1
             for j in self.nodes.keys()))

        self.pw_c2 = self.model.addConstrs(
            (self.z[j, r] <= self.M[j][r + 1]
             for j in self.nodes.keys() for r in range(self.R[j])))

        self.pw_c3 = self.model.addConstrs(
            (self.z[j, r] >= self.M[j][r]
             for j in self.nodes.keys() for r in range(self.R[j])))

    def update_variables(self):
        # piece-wide function variables
        ind_variables = [(n, t) for n, r in zip(range(self.N), self.R) for t in range(r)]

        self.z = self.model.addVars(ind_variables, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='z')
        self.u = self.model.addVars(ind_variables, vtype=GRB.BINARY, name='u')


if __name__ == "__main__":
    print("time:")
    Nodes = BOMGraph("DAG.txt").nodes
    IterativePW = IterativePW(nodes=Nodes)
    IterativePW.iter_process()
