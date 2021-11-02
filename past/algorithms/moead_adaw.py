import numpy as np
import math
from scipy.spatial.distance import cdist

from pymoo.algorithms.genetic_algorithm import GeneticAlgorithm
from pymoo.factory import get_decomposition
from pymoo.operators.crossover.simulated_binary_crossover import SimulatedBinaryCrossover
from pymoo.operators.mutation.polynomial_mutation import PolynomialMutation
from pymoo.operators.sampling.random_sampling import FloatRandomSampling
from pymoo.util.display import MultiObjectiveDisplay
from pymoo.util.misc import set_if_none
from pymoo.visualization.scatter import Scatter
from pymoo.factory import get_visualization
# =========================================================================================================
# Implementation
# =========================================================================================================


class MOEAD_AdaW(GeneticAlgorithm):
    # 初始化
    def __init__(self,
                 ref_dirs,
                 n_neighbors=20,
                 decomposition='auto',
                 prob_neighbor_mating=0.9,
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

        set_if_none(kwargs, 'pop_size', len(ref_dirs))                                  # 权重向量的个数等同于种群大小
        set_if_none(kwargs, 'sampling', FloatRandomSampling())                          #
        set_if_none(kwargs, 'crossover', SimulatedBinaryCrossover(prob=1.0, eta=20))    # 交叉算子：模拟二进制交叉，交叉概率：100%
        set_if_none(kwargs, 'mutation', PolynomialMutation(prob=None, eta=20))          # 变异算子：多项式变异
        set_if_none(kwargs, 'survival', None)
        set_if_none(kwargs, 'selection', None)

        super().__init__(display=display, **kwargs)

        # initialized when problem is known
        self.ref_dirs = ref_dirs

        if self.ref_dirs.shape[0] < self.n_neighbors:                                   # 邻居大小最大为种群大小
            print("Setting number of neighbours to population size: %s" % self.ref_dirs.shape[0])
            self.n_neighbors = self.ref_dirs.shape[0]

        # neighbours includes the entry by itself intentionally for the survival method
        # 根据权重向量间的距离排序并取最近的K个作为邻居
        self.neighbors = np.argsort(cdist(self.ref_dirs, self.ref_dirs), axis=1, kind='quicksort')[:, :self.n_neighbors]

        # 外部种群
        self.EP = []
        self.EP_size = int(1.5 * self.pop_size)
        # self.EP_size = self.pop_size
        # 小生境半径
        self.r = -1

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

        # 初始化外部种群
        for i in range(self.pop_size):
            self._update(self.pop[i])


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

            # 更新外部种群
            self._update(off)
        # 维持外部种群的大小
        if len(self.EP) > self.EP_size and self.n_gen < 180:
            self._maintain_EP()

        # 调整权重
        if self.n_gen < 180 and self.n_gen % 10 == 0:
            EP_F = [self.EP[i].get("F") for i in range(len(self.EP))]
            # get_visualization("scatter").add(np.array(EP_F)).show()
            new_p, new_w = self._add_weight()
            print(len(new_w))
            self._delete_weight(new_p, new_w)
            # 更新邻居
            self.neighbors = np.argsort(cdist(self.ref_dirs, self.ref_dirs), axis=1, kind='quicksort')[:, :self.n_neighbors]
        print(self.n_gen)
        if self.n_gen == 199:
            get_visualization("scatter").add(self.pop.get("F")).show()
            # for i in range(self.pop_size):
            #     self.pop[i] = self.EP[i]
            # get_visualization("scatter").add(self.pop.get("F")).show()

    # 更新外部种群
    def _update(self, p):
        if p in self.EP:
            return
        i = 0
        while i < len(self.EP):
            flag = 0
            for n in range(self.problem.n_obj):
                if p.get("F")[n] > self.EP[i].get("F")[n]:
                    flag += 1
            if flag == 0:
                self.EP.pop(i)
            elif flag == self.problem.n_obj:
                return
            else:
                i += 1
        self.EP.append(p)

    # 维持外部种群
    def _maintain_EP(self):
        while len(self.EP) > self.EP_size:
            crowd = self._crowd(self.EP)
            id = crowd.index(max(crowd))
            self.EP.pop(id)
        # 删除EP中的最拥挤的个体，直至EP的规模=阈值
        # if len(self.EP) > self.EP_size:
        #     l = len(self.EP) - self.EP_size
        #     # 求解个体拥挤度
        #     crowd = self._crowd(self.EP)
        #     # 找拥挤度最大的个体
        #     i = np.argsort(crowd, kind='quicksort')[:l]
        #     # 删除该个体
        #     p = []
        #     for j in range(len(self.EP)):
        #         if j not in i:
        #              p.append(self.EP[j])
        #     self.EP = p
        return

    # 计算拥挤度
    def _crowd(self, pop):
        # 外部种群的大小
        size = len(pop)
        # 外部种群个体之间的欧几里得距离（二维数组）
        dis = []
        # 计算距离
        for i in range(size):
            d = []
            for j in range(size):
                if i != j:
                    # 计算欧几里得距离
                    d.append(self._distance(pop[i].get("F"), pop[j].get("F")))
                else:
                    d.append(0)
            dis.append(d)
        dis = np.array(dis)
        # 小生境半径
        r = self._radius(dis, pop)
        # 计算拥挤度
        crowd = []
        for i in range(size):
            c = 1
            for j in range(size):
                if i != j:
                    if dis[i][j] <= r:
                        c *= dis[i][j] / r
            c = 1 - c
            crowd.append(c)
        return crowd

    # 求小生境半径
    def _radius(self, dis, pop):
        # 外部种群的大小
        ep_size = len(pop)
        # 取前m+1个最近距离的id（包括自己）
        min_d = np.argsort(dis, axis=1, kind='quicksort')[:, :self.problem.n_obj+1]
        # 取值
        min_v = []
        for i in range(ep_size):
            for j in range(1, self.problem.n_obj+1):
                min_v.append(dis[i][min_d[i][j]])
        # 小生境半径r
        r = np.median(min_v)
        self.r = r
        return r

    # 增加权重
    def _add_weight(self):
        # 新添加的个体与权重
        new_p = []
        new_w = []
        flag = True
        for i in range(len(self.EP)):
            for j in range(self.pop_size):
                d = self._distance(self.EP[i].get("F"), self.pop[j].get("F"))
                if d < self.r:
                    flag = False
                    break
            if flag:
                w = self._create_weight(self.EP[i])
                for j in range(self.pop_size):
                    fp = self._decomposition.do(self.EP[i].get("F"), weights=w, ideal_point=self.ideal_point)
                    fq = self._decomposition.do(self.pop[j].get("F"), weights=w, ideal_point=self.ideal_point)
                    if fp > fq:
                        flag = False
                        break
                    else:
                        if sum(self.EP[i].get("F")) > sum(self.pop[j].get("F")):
                            flag = False
                            break
                if flag:
                    new_p.append(self.EP[i])
                    new_w.append(w)
            flag = True
        for i in new_p:
            self.EP.pop(self.EP.index(i))
        return new_p, new_w

    # 根据个体生成权重向量
    def _create_weight(self, p):
        w = []
        v = []
        for i in range(self.problem.n_obj):
            vi = p.get("F")[i] - self.ideal_point[i]
            v.append(vi)
        for i in range(self.problem.n_obj):
            w.append(v[i]/sum(v))
        return w

    # 删除权重
    def _delete_weight(self, new_p, new_w):
        pop = []
        for i in self.pop:
            pop.append(i)
        number = len(new_w)
        for i in range(number):
            max_c = 1
            max_id = -1
            for j in range(self.pop_size):
                # 找到最多相同的个体
                if pop.count(pop[j]) > max_c:
                    max_c = pop.count(pop[j])
                    max_id = j
            # 如果没有相同的个体，则删除拥挤度最大的个体
            if max_c == 1:
                # 计算拥挤度
                crowd = self._crowd(self.pop)
                id = crowd.index(max(crowd))
                # 将拥挤度最大的个体删除，添加外部个体。
                self.pop[id] = new_p[i]
                self.ref_dirs[id] = new_w[i]
            # 如果有多个重复的个体，删除其中聚合函数值最多的个体和权重向量
            else:
                max_g = 0
                max_gid = -1
                id = -1
                # 遍历相同的个体
                for j in range(max_c):
                    # 找下一个相同个体的id
                    id = pop.index(pop[max_id], id+1)
                    # 计算聚合函数值
                    FV = self._decomposition.do(pop[id].get("F"), weights=self.ref_dirs[id], ideal_point=self.ideal_point)
                    # 找到聚合函数最大的id
                    if FV > max_g:
                        max_g = FV
                        max_gid = id
                # 删除个体、权重，添加新个体与新权重
                self.pop[max_gid] = new_p[i]
                self.ref_dirs[max_gid] = new_w[i]
        return

    # 求两点间的欧几里得距离
    def _distance(self, point1, point2):
        s = sum([(point1[i] - point2[i])**2 for i in range(len(point1))])
        return math.sqrt(s)

# parse_doc_string(MOEAD.__init__)
