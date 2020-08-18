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
import time

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
        self.M = dict(zip(self.nodes.keys(), [[0, 999, 9999]] * self.N))

        # initial model
        self.model = gp.Model("IterativeLP")
        self.model.setParam("OutputFlag", False)
        # inbound service time
        self.SI = self.model.addVars(self.N, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='SI')
        # outbound service time
        self.S = self.model.addVars(self.N, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='S')

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

            self.model.addConstr(self.SI[j] + info["lead_time"] - self.S[j], GRB.GREATER_EQUAL, 0.0, name="r[1]_"+str(j))

            if not info['sink']:
                self.model.addConstr(self.S[j], GRB.LESS_EQUAL, info['demand_service_time'], name="r[2]_"+str(j))

            for i in info['source']:
                i = int(i)
                self.model.addConstr(self.SI[j] - self.S[i], GRB.GREATER_EQUAL, 0.0, "r[3]_"+str(j)+"_"+str(i))

        self.model.update()
        # print(self.model.getConstrs())
        # print(self.model.getVars())
        self.obj = None
        self.optimal_value = None
        self.update_error = None

    def iter_process(self):
        while True:
            print("******************************* New Iteration ********************************************")
            self.iteration += 1
            self.optimization()
            print("Current iter: {}".format(self.iteration))

            if self.model.status == GRB.OPTIMAL:
                # print("Solution: \n {}".format(self.model.getVars()))
                print("Current optimal solution of true function: {}".format(self.cal_optimal_value()))
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
            holding_cost = info['holding_cost']
            val = 0
            for r in range(self.R[j]):
                val += self.u[j, r] * self.beta[j][r] + self.z[j, r] * self.alpha[j][r]
            self.obj += holding_cost * val
        self.model.setObjective(self.obj, GRB.MINIMIZE)
        self.model.optimize()
        # print(self.model.getVars())

    def termination_criterion(self):
        flag = False
        err = 0
        for j, info in self.nodes.items():
            net_replenishment_period = self.SI[j].x + info['lead_time'] - self.S[j].x
            appr_function = sum(
                [self.beta[j][r] * self.u[j, r].x + self.alpha[j][r] * self.z[j, r].x for r in range(self.R[j])])
            # if net_replenishment_period == np.nan: print(net_replenishment_period)
            err += info['holding_cost'] * abs(appr_function - true_function(max(net_replenishment_period, 0)))
        print("Current error: {} of iteration: {}".format(err, self.iteration))

        if err <= self.epsilon / self.N:
            flag = True
        if self.update_error == err:
            flag = True
            print("No reduced error, iteration ends.")
        self.update_error = err
        return flag

    def cal_optimal_value(self):
        optimal_value = 0
        for j, info in self.nodes.items():
            net_replenishment_period = round(max(self.SI[j].x + info['lead_time'] - self.S[j].x, 0), 3)
            # print(net_replenishment_period)
            # print("current j:{}, net x:{}".format(j, net_replenishment_period))
            optimal_value += info['lead_time']*true_function(net_replenishment_period)

        return optimal_value

    def update_model(self):
        self.model.remove([self.z, self.u])
        self.update_variables()

        self.model.remove([self.pw_c0, self.pw_c1, self.pw_c2, self.pw_c3])
        self.update_constraints()

        self.model.update()

    def update_para(self):
        for j, info in self.nodes.items():
            print("-----------------------")
            print("current j:{} of R{}".format(j, self.R[j]))
            print("current u:{}".format([self.u[j, r].x for r in range(self.R[j])]))

            interval_index = [r for r in range(self.R[j]) if abs(self.u[j, r].x - 1) <= 0.001][0]
            print("current index {}".format(interval_index))
            print("Current z {}".format([self.z[j, r].x for r in range(self.R[j])]))

            net_replenishment_period = self.SI[j].x + info['lead_time'] - self.S[j].x

            print("Net X: {}".format(net_replenishment_period))
            print("Previous Node: {}".format(self.M[j][interval_index]))
            print("Previous M: {}".format(self.M[j]))

            if abs(net_replenishment_period - self.M[j][interval_index]) <= 0.1:
                continue

            if abs(net_replenishment_period - self.M[j][interval_index + 1]) <= 0.1:
                continue

            self.R[j] += 1

            # not take round operator here.
            new_alpha = (true_function(net_replenishment_period) - true_function(self.M[j][interval_index])) / \
                        (net_replenishment_period - self.M[j][interval_index])
            new_beta = true_function(net_replenishment_period) - new_alpha * net_replenishment_period

            new_alpha_ = (true_function(self.M[j][interval_index + 1]) - true_function(net_replenishment_period)) / \
                         (self.M[j][interval_index + 1] - net_replenishment_period)
            new_beta_ = true_function(net_replenishment_period) - new_alpha_ * net_replenishment_period

            # dict insert;
            # notice: alpha[j][interval_index] = new_alpha will insert new value to all keys.
            # use copy method
            update_alpha = self.alpha[j].copy()
            update_alpha[interval_index] = round(new_alpha, 3)
            update_alpha.insert(interval_index + 1, round(new_alpha_, 3))
            self.alpha[j] = update_alpha

            update_beta = self.beta[j].copy()
            update_beta[interval_index] = round(new_beta, 3)
            update_beta.insert(interval_index + 1, round(new_beta_, 3))
            self.beta[j] = update_beta

            update_M = self.M[j].copy()
            update_M.insert(interval_index + 1, net_replenishment_period)
            self.M[j] = update_M

            if sorted(update_M) != update_M:
                raise DefinedException("Incorrect M")
            if len(set(update_M)) != len(update_M):
                raise DefinedException("Incorrect M")

            print("update para alpha:{}".format(self.alpha[j]))
            print("update para beta:{}".format(self.beta[j]))
            print("update para M:{}".format(self.M[j]))

    def update_constraints(self):
        self.pw_c0 = self.model.addConstrs(
            (self.SI[j] + self.nodes[j]["lead_time"] - self.S[j] == sum([self.z[j, r] for r in range(self.R[j])])
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

    def update_variables(self):
        # piece-wide function variables
        ind_variables = [(n, t) for n, r in zip(range(self.N), self.R) for t in range(r)]

        self.z = self.model.addVars(ind_variables, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name='z')
        self.u = self.model.addVars(ind_variables, vtype=GRB.BINARY, name='u')


def parse_results(instance: IterativePW) -> None:
    for j, _ in instance.nodes.items():
        SI = instance.SI[j]
        S = instance.S[j]
        print("Node: {}, SI:{}, S: {}".format(j, round(SI.x, 3), round(S.x,3)))


if __name__ == "__main__":
    Nodes = BOMGraph("DAG.txt").nodes
    # print(Nodes)
    IterativePW = IterativePW(nodes=Nodes)
    start = time.time()
    IterativePW.iter_process()
    IterativePW.cal_optimal_value()
    print("optimal value: {}".format(IterativePW.optimal_value))
    parse_results(IterativePW)
    print("Used cpu time：{}".format(time.time() - start))
    # print(IterativePW.model.getConstrs())

