# repair, crossover, mutation = self.mating.repair, self.mating.crossover, self.mating.mutation
#
#         # get_visualization("scatter").add(self.pop.get("F")).show()
#
#         # retrieve the current population
#         pop = self.pop
#
#         #------------------------------------------------------------------------------------------------------------
#         Q = self.pop.copy()
#         # if 10 < self.n_gen < 100:
#         #     Q = self._mating_pool(Q)
#
#         # 新种群，包括父代和子代
#         # offs = Population(n_individuals=2*self.pop_size)
#         offs = []
#         #------------------------------------------------------------------------------------------------------------
#
#         # iterate for each member of the population in random order
#         for i in np.random.permutation(len(pop)):
#             # for i in range(len(pop)):
#
#             # all neighbors of this individual and corresponding weights
#             N = self.neighbors[i, :]                            # 获取个体i的邻居
#             # offs[i*2] = pop[i]
#             # offs.append(pop[i])
#             # 获取一个随机数，根据其大小决定交配池是由邻居还是整个种群组成，并从中随机选择父代
#             if np.random.random() < self.prob_neighbor_mating:
#                 parents = N[np.random.permutation(self.n_neighbors)][:crossover.n_parents]
#             else:
#                 parents = np.random.permutation(self.pop_size)[:crossover.n_parents]
#
#             # if np.random.random() < self.prob_neighbor_mating:
#             #     parents = N[np.random.permutation(self.n_neighbors)][:crossover.n_parents]
#             # else:
#             #     parents = np.random.permutation(len(Q))[:crossover.n_parents]
#
#             # do recombination and create an offspring
#             off = crossover.do(self.problem, pop, parents[None, :])         # 交叉
#             # off = crossover.do(self.problem, Q, parents[None, :])           # 交叉
#             off = mutation.do(self.problem, off)                            # 变异
#             off = off[np.random.randint(0, len(off))]                       # 随机选择一个产生的个体作为后代
#
#             # repair first in case it is necessary                          # 修复个体（当超出可行域时）
#             if repair:
#                 off = self.repair.do(self.problem, off, algorithm=self)
#
#             # evaluate the offspring
#             self.evaluator.eval(self.problem, off)                          # 计算函数值
#
#             # update the ideal point
#             self.ideal_point = np.min(np.vstack([self.ideal_point, off.F]), axis=0)     # 更新参考点
#
#             # calculate the decomposed values for each neighbor
#             # FV = self._decomposition.do(pop[N].get("F"), weights=self.ref_dirs[N, :], ideal_point=self.ideal_point)     # 计算邻居聚合函数值
#             # off_FV = self._decomposition.do(off.F[None, :], weights=self.ref_dirs[N, :], ideal_point=self.ideal_point)  # 计算新个体与邻居的权重的聚合函数值
#             off_FV = self._decomposition.do(off.F[None, :], weights=self.ref_dirs[i], ideal_point=self.ideal_point)  # 计算新个体与邻居的权重的聚合函数值
#             FV = self._decomposition.do(self.pop[i].get("F"), weights=self.ref_dirs[i], ideal_point=self.ideal_point)  # 计算新个体与邻居的权重的聚合函数值
#
#             # get the absolute index in F where offspring is better than the current F (decomposed space)
#             # 更新邻居
#             I = np.where(off_FV < FV)[0]
#             pop[N[I]] = off
#             for n in I:
#                 self._update(pop[N[n]])
#
#             #**********************************************************************************************************
#             # if off_FV < FV:
#             #     pop[i] = off
#             # offs[i*2+1] = off
#             # offs.append(off)
#
#         # 将父代与子代进行帕累托分级，根据支配关系选择新一代种群
#         # 如果不足种群大小，则加入下一级，如果超过，通过密度选择最后一级
#         # P, Q = self._level(offs)
#         # get_visualization("scatter").add(P.get("F")).show()
#         # while len(P) < self.pop_size:
#         #     new_P, Q = self._level(Q)
#         #     if len(P) + len(new_P) > self.pop_size:
#         #         P = self._density_sel(P, new_P)
#         #         break
#         # for i in range(self.pop_size):
#         #     self.pop[i] = P[i]
#
#         # get_visualization("scatter").add(Q.get("F")).show()
#         # den = self._density(offs)
#         # self._s_value(offs)
#         # get_visualization("scatter").add(self.pop.get("F")).show()
#         #***************************************************************************************************************
#         # 进化后期，判断权重向量与相关联的个体角度是否差距仍较大，
#         # 较大则表明可能PF是不连续的当前权重向量应调整为个体所在角度。
#         # if self.n_gen >= 50 and self.n_gen % 20 == 0:
#         #     self._weigt_adjust()
#         #     get_visualization("scatter").add(self.ref_dirs).show()
#         #--------------------------------------------------------------------------------------------------------------

