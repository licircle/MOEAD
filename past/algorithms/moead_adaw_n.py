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
from pymoo.factory import get_visualization


# =========================================================================================================
# Implementation
# =========================================================================================================

class MOEAD_AdaW_N(GeneticAlgorithm):
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
        self.EP_size = self.pop_size
        # 小生境半径
        self.r = -1
        self.evolution = [0 for i in range(self.pop_size)]                              # 个体与上一代的优化程度
        self.flag = 0                                                                   # 不活跃代数
        self.max_angle = 10 * self._angle1(self.ref_dirs[0], ref_dirs[1])               # 权重向量与个体间最大的角度
        self.pool = [i for i in range(self.pop_size)]

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
        # for i in np.random.permutation(len(pop)):
        for i in self.pool:
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
             # 计算新个体与旧个体间的聚合函数值的差值
            for i in I:
                self.evolution[i] = FV[i] - off_FV[i]
            pop[N[I]] = off

            # 更新外部种群
            self._update(off)
        # if self.n_gen == 199:
        #     plot = Scatter()
        #     plot.add(self.pop.get("F"), color="red")
        #     plot.add(self.ref_dirs, color="blue")
        #     plot.show()

        # 维持外部种群的大小
        if len(self.EP) > self.EP_size:
            self._maintain_EP()
        s = sum(self.evolution) / self.pop_size
        if s <= 0.0001:
            self.flag += 1
        else:
            self.flag = 0
            # 调整权重
        # 每隔10代调整交配池中个体的比例
        # 将计算资源更多分配在未找到非支配解的方向
        # if self.n_gen > 100 and self.n_gen % 10 == 0:
        #     self.pool = [i for i in range(self.pop_size)]
        #     I = NonDominatedSorting().do(F=self.pop.get("F"), only_non_dominated_front=True)   # 获得种群中的非支配解
        #     for j in range(self.pop_size):
        #         if j not in I:
        #             self.pool.append(j)
        # if self.flag >= 5:
        if self.n_gen < 180 and self.n_gen % 10 == 0:
            # if self._exchange():
            #     # 更新邻居
            #     print(-1)
            #     self.neighbors = np.argsort(cdist(self.ref_dirs, self.ref_dirs), axis=1, kind='quicksort')[:, :self.n_neighbors]
            new_p, new_w = self._add_weight()
            # print(len(new_w))
            self._delete_weight(new_p, new_w)
            self.flag = 0
        print(self.n_gen)
        if self.n_gen == 199:
            get_visualization("scatter").add(self.ref_dirs).show()
            get_visualization("scatter").add(self.pop.get("F")).show()
            for i in range(self.pop_size):
                self.pop[i] = self.EP[i]
            get_visualization("scatter").add(self.pop.get("F")).show()

    # 计算种群中所有个体与相应权重向量间的角度
    # 输出角度大于最大角度的个体
    def _indiv_max(self):
        ids = []
        for i in range(self.pop_size):
            a = self._angle_FW(self.pop[i].get("F"), self.ref_dirs[i])
            if a > self.max_angle:
                ids.append(i)
        return ids

    # 进行权重向量以及个体的替换
    # 输出：如果进行了替换，返回True，否则返回False
    def _exchange(self):
        pop_ids = self._indiv_max()
        if len(pop_ids) > 0:
            new_w = []
            new_p = []
            ext = self._ext_id(len(pop_ids))
            if len(ext) > 0:
                for i in range(len(ext)):
                    new_w.append(self._create_weight(self.EP[ext[i]]))
                    new_p.append(self.EP[ext[i]])
                if len(ext) == len(pop_ids):
                    self._update_pop(new_w, new_p, pop_ids)
                elif len(ext) > 0:
                    FV = self._decomposition.do(self.pop[pop_ids].get("F"), weights=self.ref_dirs[pop_ids], ideal_point=self.ideal_point)
                    ids = np.argsort(FV, kind='quicksort')[:len(ext)]
                    p_id = []
                    for j in ids:
                        p_id.append(pop_ids[j])
                    self._update_pop(new_w, new_p, p_ids)
            return True
        else:
            return False
    # 使用外部种群中的个体替换种群中与权重向量角度大的个个体
    # 输入：将被替换的个体的个数
    # 输出：外部种群的个体id
    def _ext_id(self, length):
        ids = []
        for i in range(len(self.EP)):
            if self.EP[i] not in self.pop:
                ids.append(i)
        if len(ids) > length:
            F = [sum(self.EP[j].get("F")) for j in range(len(ids))]
            ind = np.argsort(F, kind='quicksort')[:length]
            r = []
            for j in ind:
                r.append(ids[j])
            return r
        else:
            return ids

    # 替换个体
    # 输入：新权重，新个体，旧个体
    def _update_pop(self, new_w, new_p, ids):
        for i in range(len(new_w)):
            self.ref_dirs[ids[i]] = new_w[i]
            self.pop[ids[i]] = new_p[i]
    # 更新外部种群
    def _update(self, p):
        if p in self.EP:
            return
        i = 0
        while i < len(self.EP):
            flag = 0
            for n in range(self.problem.n_obj):
                if p.get("F")[n] >= self.EP[i].get("F")[n]:
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
        crowd = self._crowd(self.EP)
        # l = len(self.EP) - self.EP_size
        # I = np.argsort(crowd, kind='quicksort')[:l]
        # F = [self.EP[i].get("F") for i in I]
        # EP_F = [self.EP[i].get("F") for i in range(len(self.EP))]
        # comp.vis_3d(EP_F, F)
        while len(self.EP) > self.EP_size:
            id = crowd.index(min(crowd))
            self.EP.pop(id)
            crowd = self._crowd(self.EP)
        # 删除EP中的最拥挤的个体，直至EP的规模=阈值
        # EP_F = [self.EP[i].get("F") for i in range(len(self.EP))]
        # get_visualization("scatter").add(np.array(EP_F)).show()
        return

    # 计算拥挤度
    def _crowd(self, pop):
        # 外部种群的大小
        size = len(pop)
        # 外部种群个体之间的欧几里得距离（二维数组）
        crowd = [0 for i in range(size)]
         # 根据每个函数值将解进行排序（从小到大）
        dim_sort = []
        for j in range(self.problem.n_obj):
            fi = [pop[i].get("F")[j] for i in range(size)]
            dim = np.argsort(fi, kind='quicksort')[:]
            dim = dim.tolist()
            dim_sort.append(dim)
         # 保存每个个体其相邻个体id
        inds = [[0 for i in range(self.problem.n_obj)] for j in range(size)]
        # 计算每个个体与其相邻非支配解的距离
        for i in range(size):
            for j in range(self.problem.n_obj):
                ind = dim_sort[j].index(i)
                if ind != 0:
                    inds[i][j] = dim_sort[j][ind - 1]
                    distance = self._distance(pop[i].get("F"), pop[dim_sort[j][ind-1]].get("F"))
                    crowd[i] += distance / self.problem.n_obj
        self.r = np.mean(crowd)
        return crowd

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

    # 计算个体与其对应的权重向量间的角度
    # 输入：个体的函数值，权重向量
    # 返回值：个体与其权重向量间的角度
    def _angle_FW(self, F, w):
        dis1 = math.sqrt(sum([(F[i]-self.ideal_point[i])**2 for i in range(self.problem.n_obj)]))
        dis2 = math.sqrt(sum([(w[i])**2 for i in range(self.problem.n_obj)]))
        s = sum([(F[i]-self.ideal_point[i]) * w[i] for i in range(self.problem.n_obj)])
        a = s / (dis1 * dis2)
        if a > 1:
            a = 1
        angle = math.degrees(math.acos(a))
        return angle

    # 计算两点与参考点形成的两个向量间的角度
    # 输入：两点
    # 输出：角度
    def _angle(self, v1, v2):
        dis1 = math.sqrt(sum([(v1[i]-self.ideal_point[i])**2 for i in range(self.problem.n_obj)]))
        dis2 = math.sqrt(sum([(v2[i]-self.ideal_point[i])**2 for i in range(self.problem.n_obj)]))
        s = sum([(v1[i]-self.ideal_point[i]) * (v2[i]-self.ideal_point[i]) for i in range(self.problem.n_obj)])
        a = s / (dis1 * dis2)
        if a > 1:
            a = 1
        angle = math.degrees(math.acos(a))
        return angle

    # 计算两权重向量间的角度
    # 输入：两个权重向量
    # 输出：权重向量间的角度
    def _angle1(self, v1, v2):
        dis1 = math.sqrt(sum([(v1[i])**2 for i in range(len(v1))]))
        dis2 = math.sqrt(sum([(v2[i])**2 for i in range(len(v2))]))
        s = sum([v1[i] * v2[i] for i in range(len(v2))])
        a = s / (dis1 * dis2)
        if a>1:
            a = 1
        angle = math.degrees(math.acos(a))
        return angle
# parse_doc_string(MOEAD.__init__)
