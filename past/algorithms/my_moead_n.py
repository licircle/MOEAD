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
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting

from pymoo.factory import get_visualization

from my_moead.Utils import Compare as comp


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
        self.flag_dis = [[0 for j in range(len(ref_dirs[0]))] for i in range(self.pop_size)]    #
        self.evolution = [0 for i in range(self.pop_size)]                              # 个体与上一代的优化程度
        self.flag = 0                                                                   # 不活跃代数

        # self.external_pop = []
        # self.min_angle = self._distance(self.ref_dirs[0], ref_dirs[1])

        self.max_angle = 10 * self._angle1(self.ref_dirs[0], ref_dirs[1])
        # self.min_angle = 5 * self._angle1(self.ref_dirs[0], ref_dirs[1])

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

        # self.last_angle = self._angle_A()

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
            # 计算新个体与旧个体间的聚合函数值的差值
            for i in I:
                self.evolution[i] = FV[i] - off_FV[i]
            pop[N[I]] = off

            #--------------新代码---------------------------------------------

        s = sum(self.evolution) / self.pop_size
        # print(s)
        if s <= 0.0001:
            self.flag += 1
        else:
            self.flag = 0
        if self.flag == 5:
            print(self.n_gen)
            self._angle_dis()
            # A = self._angle_A()
            # p_id = self._weigt_adjust(A)
            # # 调整种群
            # if len(p_id)> 0:
            #     print("len: ", len(p_id))
            #     P = self.ref_dirs[p_id]
            #     comp.vis(P)
            #     # 判断是否有断开的个体
            #     # new_w, new_p = self._dis_neighbor(len(p_id))
            #     new_w, new_p = self._sort_crowd(len(p_id))
            #     self._adjust(p_id, new_w, new_p)
            #     self.neighbors = np.argsort(cdist(self.ref_dirs, self.ref_dirs), axis=1, kind='quicksort')[:, :self.n_neighbors]

        if self.n_gen == 199:
            get_visualization("scatter").add(self.ref_dirs).show()

    #-------------------------------------------------------------------------------------------------------------------

    #
    def _angle_dis(self):
        I = NonDominatedSorting().do(F=self.pop.get("F"), only_non_dominated_front=True)   # 获得种群中的非支配解
        P = self.pop[I]
        get_visualization("scatter").add(P.get("F")).show()
        dim_sort = []
        # 根据每个函数值将解进行排序（从小到大）
        for j in range(self.problem.n_obj):
            fi = [P[i].get("F")[j] for i in range(len(P))]
            dim = np.argsort(fi, kind='quicksort')[:]
            dim = dim.tolist()
            dim_sort.append(dim)
        # 保存每个个体与相邻非支配解的角度
        angels = [[-1 for i in range(self.problem.n_obj)] for j in range(len(I))]
        # 保存距离
        distances = [[-1 for i in range(self.problem.n_obj)] for j in range(len(I))]
        # 保存每个个体其相邻个体id
        inds = [[-1 for i in range(self.problem.n_obj)] for j in range(len(I))]
        # 计算每个个体与其相邻非支配解的角度
        for i in range(len(I)):
            for j in range(self.problem.n_obj):
                ind = dim_sort[j].index(i)
                if ind != 0:
                    inds[i][j] = dim_sort[j][ind - 1]
                    angels[i][j] = self._angle(P[i].get("F"), P[dim_sort[j][ind-1]].get("F"))
                    distances[i][j] = self._distance(P[i].get("F"), P[dim_sort[j][ind-1]].get("F"))
        n = 0       # 非-1的个体数
        ang = 0     #
        dis = 0
        # 计算平均角度及平均距离
        for i in range(len(I)):
            for j in range(self.problem.n_obj):
                if angels[i][j] != -1:
                    n += 1
                    ang += angels[i][j]
                    dis += distances[i][j]
        mean_ang = ang / n
        mean_dis = dis / n
        crowds = [[1 for i in range(self.problem.n_obj)] for j in range(len(I))]
        add_area = []
        del_area = []
        Parent = []
        i1 = []
        i2 = []
        i3 = []
        i4 = []
        # 找到平均角度的个体对，判断个体间的领域是未探索到还是非连续的
        for i in range(len(I)):
            flag1 = True
            flag2 = True
            # 遍历所有个体距离和角度
            for j in range(self.problem.n_obj):
                # 保留角度或距离较大的值，相乘作为拥挤度（越大越稀疏）
                if angels[i][j] > mean_ang:
                    crowds[i][j] *= angels[i][j] / mean_ang
                if distances[i][j] > mean_dis:
                    crowds[i][j] *= distances[i][j] / mean_dis
                # 如果拥挤度大于，表明该个体与临近的非支配个体距离或角度超出了平均水平
                if crowds[i][j] > 1:
                # if 1:
                    # 判断这样的个体对间是否存在个体或是否存在权重向量
                    for k in range(self.pop_size):
                        if k != I[i] and k not in I[inds[i]]:
                            a = self._angle(P[i].get("F"), self.pop[k].get("F"))
                            if 0 < a < angels[i][j]:
                                a1 = self._angle(P[j].get("F"), self.pop[k].get("F"))
                                if 0 < a1 < angels[i][j]:
                                    if k not in Parent:
                                        Parent.append(k)
                                    flag1 = False
                    if I[i] not in self.neighbors[I[inds[i][j]]][:self.problem.n_obj+1]:
                        flag2 = False
                    if not flag1 and flag2:
                        i3.append(I[i])
                    # 如果既不存在个体也不存在权重，则表明该区域未被探索
                    elif flag1 and flag2:
                        # i1 += 1
                        if [I[i], I[inds[i][j]]] not in add_area and [I[inds[i][j]], I[i]] not in add_area:
                            add_area.append([I[i], I[inds[i][j]]])
                    # 如果没有解，但是有权重向量，则表明帕累托前沿为非连续的
                    # 找到两个个体间对应的权重向量
                    elif flag1 and not flag2:
                        # i2 += 1
                        a = self._angle1(self.ref_dirs[I[i]], self.ref_dirs[I[inds[i][j]]])
                        for k in range(1, self.n_neighbors):
                            if k != I[i] and k != I[inds[i][j]]:
                                a1 = self._angle1(self.ref_dirs[I[i]], self.ref_dirs[self.neighbors[I[i]][k]])
                                if a1 < a:
                                    a2 = self._angle1(self.ref_dirs[I[inds[i][j]]], self.ref_dirs[self.neighbors[I[i]][k]])
                                    if a2 < a and self.neighbors[I[i]][k] not in del_area:
                                        del_area.append(self.neighbors[I[i]][k])
                    else:
                        i4.append(I[i])
        # print(i1, i2, i3, i4)
        add = []
        PP = self.pop[Parent]
        # comp.vis(PP.get("F"))
        for i in range(len(add_area)):
            if add_area[i][0] not in add:
                add.append(add_area[i][0])
            if add_area[i][1] not in add:
                add.append(add_area[i][1])
        W = self.ref_dirs[del_area]
        Pop = self.pop[add]
        comp.visualize(self.pop.get("F"), PP.get("F"), P.get("F"), Pop.get("F"))
        # get_visualization("scatter").add(Pop.get("F")).show()
        # get_visualization("scatter").add(W).show()
        del_p = []
        for i in range(len(del_area)):
            if self._angle2(self.pop[del_area[i]].get("F"), self.ref_dirs[del_area[i]]) > self.max_angle:
                del_p.append(del_area[i])
        print(len(del_area))
        print(len(add))
        print(len(Parent))
        if len(add_area) > len(del_area):
            print(0)
        elif len(add_area) < len(del_area):
            fv = [self._decomposition.do(self.pop[i].get("F"), weights=self.ref_dirs[i, :], ideal_point=self.ideal_point) for i in range(len(del_area))]
            s_id = np.argsort(fv, kind='quicksort')[:len(add_area)]
            new_w, new_p = self._weight_p(add_area)
            print("change:", len(add_area))
            self._adjust(s_id, new_w, new_p)
            self.neighbors = np.argsort(cdist(self.ref_dirs, self.ref_dirs), axis=1, kind='quicksort')[:, :self.n_neighbors]

    # 新权重以及个体替换旧个体
    def _adjust(self, p_id, new_w, new_p):
        for i in range(len(p_id)):
            self.pop[p_id[i]] = new_p[i]
            self.ref_dirs[p_id[i]] = new_w[i]
            self.history_flag[i] = 0


    # 需要被调整的权重向量
    def _weigt_adjust(self, A):
        p_id = []
        for i in range(self.pop_size):
            if A[i] > self.max_angle * (self.n_gen+2)/self.n_gen:
                # self.history_flag[i] = 0
                p_id.append(i)
        # if len(p_id)>0:
        #     p = [self.pop[i].get("F") for i in p_id]
        #     re = [self.ref_dirs[i] for i in p_id]
        return p_id

    # 根据个体的权重与相邻权重生成新权重
    def _weight_p(self, p_id):
        new_w = []
        new_p = []
        for i in p_id:
            w1 = self.ref_dirs[i[0]]
            w2 = self.ref_dirs[i[1]]
            # fv = [(p1.get("F")[i] + p2.get("F")[i])/2 for i in range(self.problem.n_obj)]
            # sum_f = sum([fv[j]-self.ideal_point[j] for j in range(self.problem.n_obj)])
            # new_ref = [(fv[j] - self.ideal_point[j]) / sum_f for j in range(self.problem.n_obj)]
            new_ref = [(w1[k]+w2[k]) / 2 for k in range(self.problem.n_obj)]
            new_fv1 = self._decomposition.do(self.pop[i[0]].get("F"), weights=new_ref, ideal_point=self.ideal_point)
            new_fv2 = self._decomposition.do(self.pop[i[1]].get("F"), weights=new_ref, ideal_point=self.ideal_point)
            if new_fv1 < new_fv2:
                new_p.append(self.pop[i[0]])
            else:
                new_p.append(self.pop[i[1]])
            new_w.append(new_ref)
        return new_w, new_p

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


    # 求两点间的欧几里得距离
    def _distance(self, point1, point2):
        s = sum([(point1[i] - point2[i])**2 for i in range(len(point1))])
        return math.sqrt(s)

    def _angle1(self, v1, v2):
        dis1 = math.sqrt(sum([(v1[i])**2 for i in range(len(v1))]))
        dis2 = math.sqrt(sum([(v2[i])**2 for i in range(len(v2))]))
        s = sum([v1[i] * v2[i] for i in range(len(v2))])
        a = s / (dis1 * dis2)
        if a > 1:
            a = 1
        angle = math.degrees(math.acos(a))
        return angle

    def _angle2(self, v1, v2):
        dis1 = math.sqrt(sum([(v1[i]-self.ideal_point[i])**2 for i in range(self.problem.n_obj)]))
        dis2 = math.sqrt(sum([(v2[i])**2 for i in range(self.problem.n_obj)]))
        s = sum([(v1[i]-self.ideal_point[i]) * v2[i] for i in range(self.problem.n_obj)])
        a = s / (dis1 * dis2)
        if a > 1:
            a = 1
        angle = math.degrees(math.acos(a))
        return angle

    def _angle_A(self):
        # 保存角度
        A = []
        for i in range(self.pop_size):
            A.append(self._angle2(self.pop[i].get("F"), self.ref_dirs[i]))
        return A
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
    # 根据个体生成权重
    def _p_weight(self, p):
        sum_f = sum([p.get("F")[j]-self.ideal_point[j] for j in range(self.problem.n_obj)])
        new_ref = [(p.get("F")[j] - self.ideal_point[j]) / sum_f for j in range(self.problem.n_obj)]
        return new_ref
