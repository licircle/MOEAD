import numpy as np
from scipy.spatial.distance import cdist

from pymoo.algorithms.genetic_algorithm import GeneticAlgorithm
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.factory import get_decomposition
from pymoo.model.individual import Individual
from pymoo.model.population import Population
from pymoo.operators.crossover.simulated_binary_crossover import SimulatedBinaryCrossover
from pymoo.operators.mutation.polynomial_mutation import PolynomialMutation
from pymoo.operators.sampling.random_sampling import FloatRandomSampling
from pymoo.util.display import MultiObjectiveDisplay
from pymoo.util.misc import set_if_none

from pymoo.visualization.scatter import Scatter
from pymoo.operators.selection.tournament_selection import TournamentSelection

# =========================================================================================================
# Implementation
# =========================================================================================================

class MOEAD_AWA(GeneticAlgorithm):

    def __init__(self,
                 ref_dirs,
                 n_neighbors=20,
                 decomposition='auto',
                 prob_neighbor_mating=0.9,
                 display=MultiObjectiveDisplay(),
                 **kwargs):
        """

        MOEAD Algorithm.

        Parameters
        ----------
        ref_dirs
        n_neighbors
        decomposition
        prob_neighbor_mating
        display
        kwargs
        """

        self.n_neighbors = n_neighbors
        self.prob_neighbor_mating = prob_neighbor_mating
        self.decomposition = decomposition

        set_if_none(kwargs, 'pop_size', len(ref_dirs))
        set_if_none(kwargs, 'sampling', FloatRandomSampling())
        set_if_none(kwargs, 'crossover', SimulatedBinaryCrossover(prob=1.0, eta=20))
        set_if_none(kwargs, 'mutation', PolynomialMutation(prob=None, eta=20))
        set_if_none(kwargs, 'survival', None)
        set_if_none(kwargs, 'selection', None)

        super().__init__(display=display, **kwargs)

        # initialized when problem is known
        self.ref_dirs = ref_dirs
        self.PI = [1 for i in range(self.pop_size)]                 # 资源分配比例
        self.EP = []                                                # 外部种群
        self.last_value = [0 for i in range(self.pop_size)]         # 记录个体上一代的函数值
        self.I = [i for i in range(self.pop_size)]                  # 保存进行进化的个体的id

        if self.ref_dirs.shape[0] < self.n_neighbors:
            print("Setting number of neighbours to population size: %s" % self.ref_dirs.shape[0])
            self.n_neighbors = self.ref_dirs.shape[0]

        # neighbours includes the entry by itself intentionally for the survival method
        self.neighbors = np.argsort(cdist(self.ref_dirs, self.ref_dirs), axis=1, kind='quicksort')[:, :self.n_neighbors]

    def _initialize(self):

        if isinstance(self.decomposition, str):

            # set a string
            decomp = self.decomposition

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
        self.ideal_point = np.min(self.pop.get("F"), axis=0)

    def _next(self):
        repair, crossover, mutation = self.repair, self.mating.crossover, self.mating.mutation

        # retrieve the current population
        pop = self.pop
        number = 0

        # 计算并保留当前种群中个体的函数值
        for i in self.I:
            self.last_value[i] = self._decomposition.do(pop[i].get("F"), weights=self.ref_dirs[i], ideal_point=self.ideal_point)[0][0]
        # iterate for each member of the population in random order
        for i in np.random.permutation(len(pop)):

            # all neighbors of this individual and corresponding weights
            N = self.neighbors[i, :]

            if np.random.random() < self.prob_neighbor_mating:
                parents = N[np.random.permutation(self.n_neighbors)][:crossover.n_parents]
            else:
                parents = np.random.permutation(self.pop_size)[:crossover.n_parents]

            # do recombination and create an offspring
            off = crossover.do(self.problem, pop, parents[None, :])
            off = mutation.do(self.problem, off)
            off = off[np.random.randint(0, len(off))]

            # repair first in case it is necessary - disabled if instance of NoRepair
            off = repair.do(self.problem, off, algorithm=self)

            # evaluate the offspring
            self.evaluator.eval(self.problem, off)

            # update the ideal point
            self.ideal_point = np.min(np.vstack([self.ideal_point, off.F]), axis=0)

            # calculate the decomposed values for each neighbor
            FV = self._decomposition.do(pop[N].get("F"), weights=self.ref_dirs[N, :], ideal_point=self.ideal_point)
            off_FV = self._decomposition.do(off.F[None, :], weights=self.ref_dirs[N, :], ideal_point=self.ideal_point)

            # get the absolute index in F where offspring is better than the current F (decomposed space)
            I = np.where(off_FV < FV)[0]
            pop[N[I]] = off
            number += len(I)
        print(self.n_gen)

        # 根据每个个体的进化程度更新资源分配比例
        if self.n_gen % 50 == 0:
            self._alloc_res()
        if self.n_gen == 199:
            plot = Scatter()
            plot.add(self.pop.get("F"), color="red")
            plot.add(self.ref_dirs, color="blue")
            plot.show()

    # 根据WS方法生成新的权重向量
    def _init_ref(self):
        delata = 0
        for i in range(self.pop_size):
            wi = []
            for j in range(self.problem.n_obj):
                if self.ref_dirs[i][j] == 0:
                    delta = 0.00001
                    break
            s = sum([1/(self.ref_dirs[i][j] + delata) for j in range(self.problem.n_obj)])
            for j in range(self.problem.n_obj):
                self.ref_dirs[i][j] = 1 / self.ref_dirs[i][j] / s
        self.neighbors = np.argsort(cdist(self.ref_dirs, self.ref_dirs), axis=1, kind='quicksort')[:, :self.n_neighbors]

    # 分配计算资源
    def _alloc_res(self):
        for i in range(self.pop_size):
            new_value = self._decomposition.do(self.pop[i].get("F"), weights=self.ref_dirs[i], ideal_point=self.ideal_point)[0][0]
            d = (self.last_value[i] - new_value) / self.last_value[i]
            if d <= 0.001:
                self.PI[i] = (0.95 + 0.05 * d / 0.001) * self.PI[i]
        self.I = []
        remain = []
        for i in range(self.pop_size):
            if 1 in self.ref_dirs[i]:
                self.I.append(i)
            else:
                remain.append(i)
        PI = [self.PI[j] for j in remain]
        select = TournamentSelection(pressure=10).do(remain, int(self.pop_size/5-self.problem.n_obj), 1, PI=PI)
        for i in select:
            self.I.append(remain[i])
    def compare(self):
        return
# parse_doc_string(MOEAD.__init__)
