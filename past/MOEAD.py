import itertools
import random
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt


class MOEAD:
    def __init__(self, toolbox, cxpb, mutpb, m, lambda_, neighbor_size=2, criteria=200):
        self.toolbox = toolbox                                          # 工具
        self.population = []                                            # 种群
        self.population_size = len(lambda_)                             # 种群大小
        self.neighbor_size = neighbor_size                              # 邻居个数
        self.cxpb = cxpb                                                # 交叉概率
        self.mutpb = mutpb                                              # 变异概率
        self.m = m                                                      # 目标函数个数
        self.criteria = criteria                                        # 遗传代数
        self.weights = lambda_                                          # 权重向量
        self.z_ = []                                                    # 理想点
        self.ep = []                                                    # 外部种群
        self.min_fit = []
        self.f = []

    # 执行程序
    def execute(self):
        # 初始化
        self.init_neighbors()                                           # 初始化邻居
        self.init_population()                                          # 初始化种群
        self.init_ideal_point()                                         # 初始化参考点（理想点）
        self.init_ep()                                                  # 初始化外部种群
        x1_ = [self.population[i].fitness.values[0] for i in range(self.population_size)]
        y1_ = [self.population[i].fitness.values[1] for i in range(self.population_size)]
        plt.scatter(x1_, y1_)
        plt.show()
        t = 0
        # 开始迭代
        while self.criteria > 0:
            t += 1
            if t == 10:
                t = 0
                x1_ = [self.population[i].fitness.values[0] for i in range(self.population_size)]
                y1_ = [self.population[i].fitness.values[1] for i in range(self.population_size)]
                plt.scatter(x1_, y1_)
                plt.show()
            # 更新
            for i in range(self.population_size):
                flag = False                                            # 标志是否进行了交叉或变异

                # 交叉操作
                j, k = self.get_random()                                # 获取随机数
                offspring = self.crossover(self.population[k],
                                           self.population[j])          # 获取交叉产生的后代
                if offspring:                                           # 如果进行了交叉操作，则将标志位置为True
                    flag = True
                    fit = self.toolbox.evaluate(offspring)
                    del offspring.fitness.values
                    offspring.fitness.values = fit
                    # if self.get_fitness(offspring, self.weights[i]) < self.get_fitness(self.population[i], self.weights[i]):
                    #     self.population[i] = deepcopy(offspring)

                # 变异操作
                offspring1 = self.mutate(self.population[i])
                # 如果进行了变异操作，则使用新个体替换旧个体
                if offspring1:
                    fit = self.toolbox.evaluate(offspring1)
                    del offspring1.fitness.values
                    offspring1.fitness.values = fit
                    # if self.get_fitness(offspring, self.weights[i]) < self.get_fitness(self.population[i], self.weights[i]):
                    #     self.population[i] = deepcopy(offspring)
                    self.update_ideal_point(offspring1)             # 更新参考点
                    self.update_neighbor(self.B[i], offspring1)     # 更新邻居
                    self.update_ep(offspring1)                      # 更新外部种群
                elif flag:
                    self.update_ideal_point(offspring)             # 更新参考点
                    self.update_neighbor(self.B[i], offspring)     # 更新邻居
                    self.update_ep(offspring)                      # 更新外部种群
            fit = [self.get_fitness(self.population[i], self.weights[i]) for i in range(self.population_size)]
            # fit = [self.population[i].fitness.values[0] for i in range(self.population_size)]
            # fits = [self.population[i].fitness.values[1] for i in range(self.population_size)]
            self.min_fit.append(min(fit))
            # self.f.append(min(fits))
            self.criteria = self.criteria - 1                           # 停止条件减1
        return self.ep, self.min_fit, self.f, self.population

    # 初始化权重向量 单纯形质心设计
    def init_weights(self):
        m = self.m
        vectors = []
        c1 = 1
        c2 = 2
        vector = [0 for i in range(m)]
        numbers = [i for i in range(m)]                                 # 保存向量的下标以便生成排列 如m=2时，number=[0,1]
        for i in range(1, m+1):
            a = 1/i
            combs = list(itertools.combinations(numbers, i))            # 得到排列组合的结果
            for j in combs:
                for n in j:
                    vector[n] = a                                       # 根据结果将向量适当位置修改成对应值，获得质心向量
                vectors.append(vector)                                  # 将得到的向量添加到向量集合中
                vector = [0 for i in range(m)]                          # 将向量重新置零
        self.weights = vectors
        self.population_size = len(vectors)                             # 种群大小与权重向量个数相同

    # 初始化邻居
    def init_neighbors(self):
        vectors = self.weights                                          # 将权重向量赋值给vectors
        K = self.neighbor_size                                          # 将邻居个数赋值给K
        B = []                                                          # 初始化邻居矩阵
        v = [np.array(vectors[i]) for i in range(len(vectors))]
        # 计算所有向量之间的距离
        distances = [
            [np.sqrt(np.sum((v[i] - v[j])**2)) for i in range(len(vectors))]
            for j in range(len(vectors))]
        for i in range(len(distances)):
            distances[i][i] = 1                                         # 向量与自身的距离置为1，以免影响找其他最近的向量
        # 根据距离矩阵选择最近的邻居
        for distance in distances:
            b = []
            for i in range(K):                                          # 找K个邻居
                index = distance.index(min(distance))                   # 找到当前离个体最近的点的下标
                b.append(index)                                         # 将下标添加到b
                distance[index] = 1                                     # 将最小距离改为1，避免重复选择
            B.append(b)                                                 # 将代表邻居的数组添加到B中
        self.B = B                                                      # 得到所有向量的邻居

    # 初始化种群
    def init_population(self):
        self.population = self.toolbox.population(n=self.population_size)
        fitness = map(self.toolbox.evaluate, self.population)
        for ind, fit in zip(self.population, fitness):
            ind.fitness.values = fit                                    # 函数值

    # 初始化参考点
    def init_ideal_point(self):
        z = [0 for i in range(self.m)]                                  # 初始化参照点
        # 得到参考点
        for i in range(self.m):
            zi = [fit.fitness.values[i] for fit in self.population]                    # 将每个个体的第i个目标函数值，即得到目标函数值矩阵的第i列
            z[i] = min(zi)                                              # 得到所有个体中第i个目标函数最小的值
        self.z_ = z                                                     # 赋值

    # 初始化外部种群
    def init_ep(self):
        for ind in self.population:
            self.update_ep(ind)

    # 计算适应度
    # 输入：个体，对应的权重向量
    # 输出：个体适应度
    def get_fitness(self, individual, lambda_):
        values = individual.fitness.values                              # 获取个体的目标函数值
        z_ = self.z_                                                    # 获取参考点
        fit_list = []                                                   # 初始化列表，用于保存每个目标函数值与参考点的差值与权重向量的乘积
        # 使用切比雪夫法计算适应度
        for i in range(self.m):
            fit = abs(values[i] - z_[i]) * lambda_[i]                   # 计算第i个目标函数值与参考点的差值与权重向量的乘积
            fit_list.append(fit)                                        # 将得到的值添加到列表中
        fitness = max(fit_list)                                         # 得到个体的适应度
        fitness = sum([values[i] * lambda_[i] for i in range(self.m)])  # 加权法求适应度
        return fitness

    # 产生两个不相同的随机数
    # 输出：两个随机数
    def get_random(self):
        i = j = random.randint(0, self.neighbor_size)                    # 先产生一个随机数
        # 在产生另一个不同的随机数前，不能跳出循环
        while i == j:
            j = random.randint(0, self.neighbor_size)                   # 产生另一个随机数
        return i, j

    # 交叉
    # 输入： 两个父代个体
    # 输出： 新个体
    def crossover(self, f1, f2):
        r = random.random()                                             # 产生一个0-1的随机数
        # 如果随机数小于等于交叉概率，则进行交叉操作
        if r <= self.cxpb:
            offspring = self.toolbox.mate(f1, f2)                       # 根据交叉得到个体
            return offspring[0]
        # 否则不进行交叉操作
        else:
            return None

    # 变异
    # 输入：父个体
    # 输出：新个体
    def mutate(self, individual):
        r = random.random()                                             # 产生一个0-1的随机数
        # 如果随机数小于等于变异概率，则进行变异操作
        if r < self.mutpb:
            offspring = self.toolbox.mutate(individual, 21, 0, 1, 0.5)
            return offspring[0]
        # 否则不进行变异操作
        return None

    # 更新参考点
    # 输入：新个体
    def update_ideal_point(self, individual):
        value = individual.fitness.values                               # 得到新个体的目标函数值
        # 更新参考点
        for i in range(self.m):
            # 如果新个体的第i个目函数比参考点的第i个值更小，则更新参考点
            if value[i] < self.z_[i]:
                self.z_[i] = value[i]

    # 更新相邻解
    # 输入：相邻解， 新个体
    def update_neighbor(self, neighbors, individual):
        n = 0
        for i in range(self.neighbor_size):
            weight = self.weights[i]                                    # 得到当前这个邻居的权重向量
            # 如果新个体的适应度值小于邻居i，则更新邻居i
            if self.get_fitness(individual, weight) - self.get_fitness(self.population[neighbors[i]], weight) > 0.00001 and n <= 1:
                self.population[neighbors[i]] = deepcopy(individual)
                n += 1

    # 更新外部种群
    # 输入：新个体
    def update_ep(self, individual):
        number = 0                                                      # 记录不支配individual的个体数目
        ids = []                                                        # 记录被支配的id
        for ind in self.ep:
            m = 0                                                       # 记录ind比individual差的函数个数
            for i in range(self.m):
                if individual.fitness.values[i] < ind.fitness.values[i]:
                    m += 1
            # 如果individual的每个函数值都优于ind，那么就删除ind, 添加number
            if m == self.m:
                ids.append(self.ep.index(ind))
            elif m > 0:
                number += 1
        for i in range(0, len(ids)):
            self.ep.pop(ids[-i-1])
        if number == len(self.ep):
            self.ep.append(individual)

# parse_doc_string(MOEAD.__init__)
