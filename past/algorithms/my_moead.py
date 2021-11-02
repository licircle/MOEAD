import numpy as np
from scipy.spatial.distance import cdist
import math

from pymoo.algorithms.genetic_algorithm import GeneticAlgorithm
from pymoo.factory import get_decomposition
from pymoo.operators.crossover.simulated_binary_crossover import SimulatedBinaryCrossover
from pymoo.operators.mutation.polynomial_mutation import PolynomialMutation
from pymoo.operators.sampling.random_sampling import FloatRandomSampling
from pymoo.util.display import MultiObjectiveDisplay
from pymoo.util.misc import set_if_none
from pymoo.model.population import Population

from pymoo.factory import get_visualization
from pymoo.visualization.scatter import Scatter

import matplotlib.pyplot as plt
# =========================================================================================================
# Implementation
# =========================================================================================================

class M_MOEAD(GeneticAlgorithm):
    # 初始化
    def __init__(self,
                 ref_dirs,
                 n_neighbors=20,
                 decomposition='auto',
                 prob_neighbor_mating=0.9,
                 max_angle = 5,
                 min_angle = 1,
                 display=MultiObjectiveDisplay(),
                 **kwargs):
        """

        Parameters
        ----------
        ref_dirs                    # 权重向量
        n_neighbors                 # 邻居
        decomposition               # 分解方法
        prob_neighbor_mating        # 从邻居中选择父代的概率
        display
        kwargs
        """

        self.n_neighbors = n_neighbors
        self.prob_neighbor_mating = prob_neighbor_mating
        self.decomposition = decomposition
        # self.external_pop = []

        set_if_none(kwargs, 'pop_size', len(ref_dirs))                                  # 权重向量的个数等同于种群大小
        set_if_none(kwargs, 'sampling', FloatRandomSampling())                          #
        set_if_none(kwargs, 'crossover', SimulatedBinaryCrossover(prob=1.0, eta=20))    # 交叉算子：模拟二进制交叉，交叉概率：100%
        set_if_none(kwargs, 'mutation', PolynomialMutation(prob=None, eta=20))          # 变异算子：多项式变异
        set_if_none(kwargs, 'survival', None)
        set_if_none(kwargs, 'selection', None)

        super().__init__(display=display, **kwargs)

        # initialized when problem is known
        self.ref_dirs = ref_dirs                                                        # 权重向量
        self.last_angle = [0 for i in range(self.pop_size)]                             # 上一代种群与权重向量的角度
        self.history_flag = [0 for i in range(self.pop_size)]                           # 角度差符合最小值的连续代数
        self.flag_dis = [[0 for j in range(len(ref_dirs[0]))] for i in range(self.pop_size)]
        # self.external_pop = []
        # self.min_angle = self._distance(self.ref_dirs[0], ref_dirs[1])

        self.max_angle = 10 * self._angle1(self.ref_dirs[0], ref_dirs[1])               # 权重向量与个体间的最大角度
        self.min_angle = 5 * self._angle1(self.ref_dirs[0], ref_dirs[1])                # 最小角

        if self.ref_dirs.shape[0] < self.n_neighbors:                                   # 邻居大小最大为种群大小
            print("Setting number of neighbours to population size: %s" % self.ref_dirs.shape[0])
            self.n_neighbors = self.ref_dirs.shape[0]

        # neighbours includes the entry by itself intentionally for the survival method
        # 根据权重向量间的距离排序并取最近的K个作为邻居
        self.neighbors = np.argsort(cdist(self.ref_dirs, self.ref_dirs), axis=1, kind='quicksort')[:, :self.n_neighbors]

    def _initialize(self):

        if isinstance(self.decomposition, str):

            # set a string
            decomp = self.decomposition

            # 根据目标函数个数选择分解方法
            # for one or two objectives use tchebi otherwise pbi
            if decomp == 'auto':
                if self.problem.n_obj <= 2:
                    decomp = 'tchebi'
                else:
                    decomp = 'pbi'

            # set the decomposition object
            self._decomposition = get_decomposition(decomp)

        else:
            self._decomposition = self.decomposition

        super()._initialize()
        # 设置参考点
        self.ideal_point = np.min(self.pop.get("F"), axis=0)

        self.last_angle = self._angle_A()
        #---------------------------------------------------------------------------------------------------------------
        # 调整个体与权重向量的关联关系
        # for i in range(0, len(self.ref_dirs)):
        #     i_FV = self._decomposition.do(self.pop[i:].get("F"), weights=self.ref_dirs[i, :], ideal_point=self.ideal_point)
        #     min_w = np.argmin(i_FV)
        #     temp = self.pop[i]
        #     self.pop[i] = self.pop[min_w]
        #     self.pop[min_w] = temp
        #---------------------------------------------------------------------------------------------------------------

    def _next(self):
        repair, crossover, mutation = self.mating.repair, self.mating.crossover, self.mating.mutation

        # retrieve the current population
        pop = self.pop

        # iterate for each member of the population in random order
        for i in np.random.permutation(len(pop)):

            # all neighbors of this individual and corresponding weights
            N = self.neighbors[i, :]                            # 获取个体i的邻居

            # 获取一个随机数，根据其大小决定交配池是由邻居还是整个种群组成，并从中随机选择父代
            if np.random.random() < self.prob_neighbor_mating:
                parents = N[np.random.permutation(self.n_neighbors)][:crossover.n_parents]
            else:
                parents = np.random.permutation(self.pop_size)[:crossover.n_parents]


            # do recombination and create an offspring
            off = crossover.do(self.problem, pop, parents[None, :])         # 交叉
            off = mutation.do(self.problem, off)                            # 变异
            off = off[np.random.randint(0, len(off))]                       # 随机选择一个产生的个体作为后代

            # repair first in case it is necessary                          # 修复个体（当超出可行域时）
            if repair:
                off = self.repair.do(self.problem, off, algorithm=self)

            # evaluate the offspring
            self.evaluator.eval(self.problem, off)                          # 计算函数值

            # update the ideal point
            self.ideal_point = np.min(np.vstack([self.ideal_point, off.F]), axis=0)     # 更新参考点

            # calculate the decomposed values for each neighbor
            FV = self._decomposition.do(pop[N].get("F"), weights=self.ref_dirs[N, :], ideal_point=self.ideal_point)     # 计算邻居聚合函数值
            off_FV = self._decomposition.do(off.F[None, :], weights=self.ref_dirs[N, :], ideal_point=self.ideal_point)  # 计算新个体与邻居的权重的聚合函数值

            # get the absolute index in F where offspring is better than the current F (decomposed space)
            # 更新邻居
            I = np.where(off_FV < FV)[0]
            pop[N[I]] = off

            #--------------新代码---------------------------------------------
            # n = np.random.random()

            # 根据个体是否更新作为个体活跃程度
            # for i in I:
            #     self.flag_dis[i] = 1


            # for i in I:
            #     self._update(pop[N[i]])
        # for i in range(self.pop_size):
        #     if self.flag_dis[i] == 0:
        #         self.history_flag[i] += 1
        #     else:
        #         self.history_flag[i] = 0
        #     self.flag_dis[i] = 0

        print(self.n_gen)

        # n = np.random.random()
        # if n < (self.n_gen/201):
        #     new_w, new_p = self._found()
        #     if len(new_p) > 0:
        #         self._delete_weight(new_w, new_p)

        #-----------------------------------------------------新代码----------------------------------------------------
        # 版本1
        # # 求小生境半径
        # r = self._radius_R()
        # # 求解种群中个体的密度
        # D = self._crowd(r)
        # # 求解角度
        # A = self._angle_A()
        # if self.n_gen >= 10 and self.n_gen <= 100:
        #
        #     # 替换个体
        #     max_angle_id, ex_id = self._exchange_id(D, A)
        #     e = self.external_pop[ex_id]
        #     self.external_pop.pop(ex_id)
        #     self._update(self.pop[max_angle_id])
        #     self.pop[max_angle_id] = e
        #
        #
        # # 调整权重向量
        # if self.n_gen >= 100 and self.n_gen % 20 == 0:
        #     self._weigt_adjust(A)
        #     get_visualization("scatter").add(self.ref_dirs).show()
        #---------------------------------------------------------------------------------------------------------------

        # 版本2
        # 计算个体与权重向量的角度
        # 更新个体与邻居的距离的标志位
        # self._get_add_area()
        A = self._angle_A()
        for i in range(self.pop_size):
            a = math.fabs(self.last_angle[i]-A[i])
            if a < self.min_angle:
                self.history_flag[i] += 1
            else:
                self.history_flag[i] = 0

        # 获得角度超过阈值的个体id
        p_id = self._weigt_adjust(A)
        self.last_angle = A
        # 调整种群
        if len(p_id)> 0:
            print("len: ", len(p_id))
            # 判断是否有断开的个体
            # new_w, new_p = self._dis_neighbor(len(p_id))
            new_w, new_p = self._sort_crowd(len(p_id))
            self._adjust(p_id, new_w, new_p)
            self.neighbors = np.argsort(cdist(self.ref_dirs, self.ref_dirs), axis=1, kind='quicksort')[:, :self.n_neighbors]
        if self.n_gen == 199:
            get_visualization("scatter").add(self.ref_dirs).show()
    #-------------------------------------------------------------------------------------------------------------------
    # 支配关系
    # 如果F1不被F2支配，则返回True
    #反之则返回False
    def _dominate(self, F1, F2):
        for i in range(self.problem.n_obj):
            if F1[i] < F2[i]:
                return True
        return False


    # 新权重以及个体替换旧个体
    def _adjust(self, p_id, new_w, new_p):
        for i in range(len(p_id)):
            self.pop[p_id[i]] = new_p[i]
            self.ref_dirs[p_id[i]] = new_w[i]
            self.history_flag[i] = 0
            self.last_angle[p_id[i]] = 0



    # 需要被调整的权重向量
    def _weigt_adjust(self, A):
        p_id = []
        for i in range(self.pop_size):
            if A[i] > self.max_angle * (self.n_gen+2)/self.n_gen and self.history_flag[i] >= 20:
                self.history_flag[i] = 0
                p_id.append(i)
        if len(p_id)>0:
            print("p_id", len(p_id))
            p = [self.pop[i].get("F") for i in p_id]
            re = [self.ref_dirs[i] for i in p_id]
            # plot = Scatter()
            # plot.add(self.pop.get("F"), color="red")
            # plot.add(self.ref_dirs, color="blue")
            # plot.show()
            # plot.add(np.array(p), color="green")
            # plot.add(np.array(re), color="yellow")
            # plot.show()
        return p_id

    # 根据个体的拥挤度调整权重
    def _sort_crowd(self, length):
        r = self._radius_R()
        D = self._crowd(r)
        min_list = np.argsort(D, kind='quicksort')[:length]
        new_p = []
        new_w = []
        for i in min_list:
            new_w.append(self._p_weight(self.pop[i]))
            new_p.append(self.pop[i])
        return new_w, new_p

    # 根据个体生成权重
    def _p_weight(self, p):
        sum_f = sum([p.get("F")[j]-self.ideal_point[j] for j in range(self.problem.n_obj)])
        new_ref = [(p.get("F")[j] - self.ideal_point[j]) / sum_f for j in range(self.problem.n_obj)]
        return new_ref

    # 根据个体的权重与相邻权重生成新权重
    def _weight_p(self, p_id):
        p1 = self.pop[p_id]
        p2 = self.pop[self.neighbors[p_id][1]]
        fv = [(p1.get("F")[i] + p2.get("F")[i])/2 for i in range(self.problem.n_obj)]
        sum_f = sum([fv[j]-self.ideal_point[j] for j in range(self.problem.n_obj)])
        new_ref = [(fv[j] - self.ideal_point[j]) / sum_f for j in range(self.problem.n_obj)]
        return new_ref

    # 求每个个体与其相对应的权重向量的角度
    def _angle_A(self):
        # 保存角度
        A = []
        for i in range(self.pop_size):
            A.append(self._angle(self.pop[i].get("F"), self.ref_dirs[i]))
        return A

    # 计算两点与参考点形成的两个向量间的角度
    # 输入：两点
    # 输出：角度
    def _angle(self, v1, v2):
        dis1 = math.sqrt(sum([(v1[i]-self.ideal_point[i])**2 for i in range(self.problem.n_obj)]))
        dis2 = math.sqrt(sum([(v2[i])**2 for i in range(self.problem.n_obj)]))
        s = sum([(v1[i]-self.ideal_point[i]) * v2[i] for i in range(self.problem.n_obj)])
        a = s / (dis1 * dis2)
        if a > 1:
            a = 1
        angle = math.degrees(math.acos(a))
        return angle

    def _angle1(self, v1, v2):
        dis1 = math.sqrt(sum([(v1[i])**2 for i in range(len(v1))]))
        dis2 = math.sqrt(sum([(v2[i])**2 for i in range(len(v2))]))
        s = sum([v1[i] * v2[i] for i in range(len(v2))])
        angle = math.degrees(math.acos(s / (dis1 * dis2)))
        return angle

    # 求个体的拥挤度
    # 输入：小生境半径
    def _crowd(self, r):
        D = [1 for i in range(self.pop_size)]
        for i in range(self.pop_size):
            for j in range(i+1, self.pop_size):
                d = self._distance(self.pop[i].get("F"), self.pop[j].get("F"))
                if 0 < d < r:
                    D[i] *= d / r
                    D[j] *= d / r
                elif d == 0:
                    D[i] *= 0.01 * r
                    D[j] *= 0.01 * r
        for i in range(self.pop_size):
            D[i] = 1 - D[i]
        return D


    # 计算半径
    # 输入：种群
    # 输出：集合中所有解到其第k个最近解的距离的中值
    def _radius_R(self):
        # 保存个体间的欧式距离
        dis = []
        # 计算距离
        for i in range(self.pop_size):
            for j in range(self.n_neighbors):
                d = self._distance(self.pop[i].get("F"), self.pop[self.neighbors[i][j]].get("F"))
                dis.append(d)
        # 以距离中的平均值作为小生境半径
        r = np.median(dis)
        return r


    # 求两点间的欧几里得距离
    def _distance(self, point1, point2):
        s = sum([(point1[i] - point2[i])**2 for i in range(len(point1))])
        return math.sqrt(s)


# parse_doc_string(MOEAD.__init__)
