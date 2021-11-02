from deap import creator, base, tools
import math
import numpy as np
import matplotlib.pyplot as plt
from . import MOEAD
# 邻居
NER_ITEM = 20
# 种群大小
MAX_ITEM = 50
# 权重向量个数
MAX_WEIGHT = 50


# 产生随机浮点数
# 输入：上界、下界
# 输出：介于上下界间的随机浮点数
def rnd(i=0, j=1):
    # random.seed(64)
    f = np.random.random() * abs(i - j) + i
    # print(f)
    return f


# SCH问题的目标函数值计算
def SCH(individual):
    x_ = individual[0]
    f2_ = (x_ - 5) ** 2
    if x_ <= 1:
        f1_ = -x_
    elif 1 < x_ < 3:
        f1_ = -2 + x_
    elif 3 < x_ <= 4:
        f1_ = 4-x_
    else:
        f1_ = -4 + x_
    return f1_, f2_


# DEB问题的目标函数值计算
def DEB(individual):
    x_value = individual[0]
    y_value = individual[1]
    f1_ = x_value
    f2_ = (1+10*y_value)*(1-(x_value/(1+10*y_value))**2-(x_value/(1+10*y_value)*math.sin(8*math.pi*x_value)))
    return f1_, f2_


# KUR问题的目标函数值计算
def KUR(individual):
    f1_ = sum([-10 * math.e ** ((-0.2) * math.sqrt(individual[i] ** 2 + individual[i+1] ** 2))
               for i in range(2)])
    f2_ = sum([abs(individual[i]) ** 0.8 + 5 * math.sin(individual[i] ** 3)
               for i in range(3)])
    return f1_, f2_


# ZDT1问题的目标函数值计算
def ZDT1(individual):
    x_ = [individual[i] for i in range(len(individual))]
    f1_ = x_[0]
    g_ = 1 + 9 * sum(x_[1:]) / (len(x_)-1)
    f2_ = g_ * (1 - math.sqrt(f1_ / g_))
    return f1_, f2_


# ZDT2问题
def ZDT2(individual):
    x_ = [individual[i] for i in range(len(individual))]
    f1_ = x_[0]
    g_ = 1 + 9 * sum(x_[1:]) / (len(x_)-1)
    f2_ = g_ * (1 - (f1_ / g_) ** 2)
    return f1_, f2_


# ZDT3问题
def ZDT3(individual):
    x_ = [individual[i] for i in range(len(individual))]
    f1_ = x_[0]
    g_ = 1 + 9 * sum(x_[1:]) / (len(x_)-1)
    f2_ = g_ * (1 - math.sqrt(f1_ / g_) - x_[0] / g_ * math.sin(10 * math.pi * x_[0]))
    return f1_, f2_


# ZDT4问题
def ZDT4(individual):
    x_ = [individual[i] for i in range(len(individual))]
    f1_ = x_[0]
    g_ = 1 + 10 * len(x_)-1 + sum([x_[i] ** 2 - 10 * math.cos(4 * math.pi * x_[i]) for i in range(1, len(x_))])
    f2_ = g_ * (1 - math.sqrt(f1_ / g_))
    return f1_, f2_


# DTLZ1问题
def DTLZ1(individual):
    x_ = [individual[i] for i in range(len(individual))]
    g_ = 100 * (len(x_) - 2) + 100 * sum([(x_[i] - 0.5) **2 - math.cos(20 * math.pi * (x_[i] * 0.5))
                                          for i in range(2, len(x_))])
    f1_ = (1 + g_) * x_[0] * x_[1]
    f2_ = (1 + g_) * x_[0] * (1 - x_[1])
    f3_ = (1 + g_) * (1 - x_[0])
    return f1_, f2_, f3_


# DTLZ2问题
def DTLZ2(individual):
    x_ = [individual[i] for i in range(len(individual))]
    g_ = sum([x_[i] ** 2 for i in range(2, len(x_))])
    f1_ = (1 + g_) * math.cos(x_[0] * math.pi / 2) * math.cos(x_[1] * math.pi / 2)
    f2_ = (1 + g_) * math.cos(x_[0] * math.pi / 2) * math.sin(x_[1] * math.pi / 2)
    f3_ = (1 + g_) * math.sin(x_[0] * math.pi / 2)
    return f1_, f2_, f3_


# 计算权重向量
def get_vectors():
    vectors = []
    H = 120
    for i in range(1, H):
        vectors.append([i/H, (H-i)/H])
    return vectors


# 读取数据
def get_data(name):
    data = []
    f_name = './data/' + name + '.txt'                              # 根据需要的数据集名称合成相应文件路径
    f = open(f_name, 'r')                                           # 打开文件
    for line in f:                                                  # 按行读取数据
        l = line.split()                                            # 将一行的数据按空格划分
        data.append([float(l[i]) for i in range(len(l))])           # 将数据转化为浮点数添加到数组中
    f.close()                                                       # 关闭文件
    return data


# 计算两点间的欧几里得距离
def e_distance(point1, point2):
    dist = math.sqrt(sum([(point1[i] - point2[i])**2
                          for i in range(len(point1))]))
    return dist


# 计算 D-metric
def d_metric(pareto, pop):
    distance = [min([e_distance(pareto[i], [pop[j].fitness.values[0], pop[j].fitness.values[1]])
                    for j in range(len(pop))])
                for i in range(len(pareto))]                        # 计算理想帕累托前沿中的点与结果中的所有点的距离，取最小值
    result = sum(distance) / len(pareto)                            # 根据公式得到最后的值
    return result


def main(evaluate, objective=2, individual_size=2, x_min=0, x_max=1):
    if objective == 2:
        creator.create("Fitness", base.Fitness, weights=(-1.0, -1.0))
    elif objective == 3:
        creator.create("Fitness", base.Fitness, weights=(-1.0, -1.0, -1.0))
    else:
        print("objective number error")
        exit(-1)

    creator.create("Individual", list, fitness=creator.Fitness)
    toolbox = base.Toolbox()
    toolbox.register("attr_item", rnd, x_min, x_max)                # 属性生成器(邻居)
    toolbox.register("individual", tools.initRepeat, creator.Individual,
                     toolbox.attr_item, n=individual_size)          # 注册个体
    toolbox.register("population", tools.initRepeat,
                     list, toolbox.individual)                      # 注册种群
    toolbox.register("evaluate", evaluate)                          # 注册评价（目标）函数
    toolbox.register("mate", tools.cxTwoPoint)                      # 注册交叉操作
    toolbox.register("mutate", tools.mutPolynomialBounded)          # 注册变异操作

    v = get_vectors()                                               # 得到权重向量
    x = [v[i][0] for i in range(len(v))]                            # 绘制权重向量
    y = [v[i][1] for i in range(len(v))]
    plt.scatter(x, y)
    plt.show()
    ea = MOEAD(toolbox, 0.8, 0.2, objective, v, criteria=200)
    ep, fit, f, pop = ea.execute()
    return ep, fit, f, pop


if __name__ == "__main__":
    # ep, fit, f, pop = main(DEB, 2, 2)
    # ep, fit, f, pop = main(KUR, 2, 3, x_max=5, x_min=-5)
    ep, fit, f, pop = main(ZDT1, 2, 30)
    # ep, fit, f, pop = main(ZDT2, 2, 30)
    # ep, fit, f, pop = main(ZDT3, 2, 30)
    plt.plot(fit)                                                   # 绘制算法过程中最小适应度的变化
    plt.show()

    x1_ = [pop[i][0] for i in range(len(pop))]
    y1_ = [pop[i][1] for i in range(len(pop))]
    plt.scatter(x1_, y1_)
    plt.show()

    data = get_data('ZDT1')                                         # 获取理想帕累托前沿数据
    p1 = [d[0] for d in data]
    p2 = [d[1] for d in data]
    plt.scatter(p1, p2)
    f1 = [ep[i].fitness.values[0] for i in range(len(ep))]          # 获取所有个体第一个目标函数值
    f2 = [ep[i].fitness.values[1] for i in range(len(ep))]          # 获取所有个体第二个目标函数值
    pop1 = [pop[i].fitness.values[0] for i in range(len(pop))]          # 获取所有个体第一个目标函数值
    pop2 = [pop[i].fitness.values[1] for i in range(len(pop))]          # 获取所有个体第二个目标函数值
    # plt.ylim(-0.5, 1)                                                 # 设置y轴的范围
    plt.scatter(f1, f2)                                             # 绘制得到的帕累托前沿
    plt.scatter(pop1, pop2)
    plt.show()

    # print(d_metric(data, ep))2
