import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def get_data(name):
    data = []
    # f_name = 'D://Study/Code/python/MOEAD-master/src/vector_csv_file/' + name + '.txt'                              # 根据需要的数据集名称合成相应文件路径
    f_name = './data/ZDT1.txt'
    f = open(f_name, 'r')                                           # 打开文件
    for line in f:                                                  # 按行读取数据
        l = line.split()                                            # 将一行的数据按空格划分
        data.append([float(l[i]) for i in range(len(l))])           # 将数据转化为浮点数添加到数组中
    f.close()                                                       # 关闭文件
    return data


# 计算 D-metric
def d_metric(pareto, pop):
    distance = [min([e_distance(pareto[i], pop[j])
                    for j in range(len(pop))])
                for i in range(len(pareto))]                        # 计算理想帕累托前沿中的点与结果中的所有点的距离，取最小值
    result = sum(distance) / len(pareto)                            # 根据公式得到最后的值
    return result


# 计算两点间的欧几里得距离
def e_distance(point1, point2):
    dist = sum([(point1[i] - point2[i])**2
                          for i in range(len(point1))])
    return math.sqrt(dist)
    # return dist


def compute(pop):
    pareto = get_data("ZDT1")
    return d_metric(pareto, pop)


def draw(data):
    if data.dsim == 1:
        plt.scatter()
    if data.dsim == 2:
        plt.scatter()


def vis(pop):
    pop_x = []
    pop_y = []
    for i in range(len(pop)):
        pop_x.append(pop[i][0])
        pop_y.append(pop[i][1])
        # plt.plot([0, weight[i][0]], [0, weight[i][1]], color='r', linestyle='-')
    plt.scatter(pop_x, pop_y, s=20, c='b', alpha=0.5)
    plt.show()


def vis1(pop, p):
    pop_x = []
    pop_y = []
    for i in range(len(pop)):
        pop_x.append(pop[i][0])
        pop_y.append(pop[i][1])
    plt.scatter(pop_x, pop_y, s=50, c='b', alpha=0.5)
    pop_x = []
    pop_y = []
    for i in range(len(p)):
        pop_x.append(p[i][0])
        pop_y.append(p[i][1])
        # plt.plot([0, weight[i][0]], [0, weight[i][1]], color='r', linestyle='-')
    plt.scatter(pop_x, pop_y, s=20, c='r', alpha=1)
    plt.show()

def vis_3d(pop, p):
    #定义坐标轴
    fig = plt.figure()
    ax1 = plt.axes(projection='3d')
    xd = []
    yd = []
    zd = []
    for i in range(len(pop)):
        xd.append(pop[i][0])
        yd.append(pop[i][1])
        zd.append(pop[i][2])
    ax1.scatter3D(xd, yd, zd, s=50, cmap='Blues', alpha=0.5)  #绘制散点图
    xd = []
    yd = []
    zd = []
    for i in range(len(p)):
        xd.append(p[i][0])
        yd.append(p[i][1])
        zd.append(p[i][2])
    ax1.scatter3D(xd, yd, zd, s=20, cmap='r', alpha=1)  #绘制散点图
    plt.show()


def visualize(pop, p1, p2, p3):
    # pop_x = []
    # pop_y = []
    # for i in range(len(pop)):
    #     pop_x.append(pop[i][0])
    #     pop_y.append(pop[i][1])
    #     # plt.plot([0, weight[i][0]], [0, weight[i][1]], color='r', linestyle='-')
    # plt.scatter(pop_x, pop_y, s=35, c='b', alpha=0.5)
    # pop_x = []
    # pop_y = []
    # for i in range(len(p1)):
    #     pop_x.append(p1[i][0])
    #     pop_y.append(p1[i][1])
    # plt.scatter(pop_x, pop_y, s=25, c='r', alpha=0.5)
    pop_x = []
    pop_y = []
    for i in range(len(p2)):
        pop_x.append(p2[i][0])
        pop_y.append(p2[i][1])
    plt.scatter(pop_x, pop_y, s=70, c='g', alpha=0.5)
    pop_x = []
    pop_y = []
    for i in range(len(p3)):
        pop_x.append(p3[i][0])
        pop_y.append(p3[i][1])
    plt.scatter(pop_x, pop_y, s=30, c='r', alpha=1)
    plt.show()
    # plt.scatter(p_x, p_y, s=20, c='y',alpha=0.5)
