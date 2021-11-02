# # 计算每个个体在每个目标函数上最近的个体
#     def _found(self):
#         get_visualization("scatter").add(self.pop.get("F")).show()
#         # 保存新个体和新权重
#         new_w = []
#         new_p = []
#         # 保存每个个体每个目标函数维度与其差距最小的id与值
#         id_min = [[ 0 for j in range(self.problem.n_obj)]for i in range(self.pop_size)]
#         all_d = [[0 for j in range(self.problem.n_obj)] for i in range(self.pop_size)]
#         # 计算差值
#         for i in range(self.pop_size):
#             for k in range(self.problem.n_obj):
#                 min_m = 100000
#                 min_id = -1
#                 for j in range(self.n_neighbors):
#                     s = self.pop[i].get("F")[k] - self.pop[self.neighbors[i][j]].get("F")[k]
#                     if 0 < s < min_m:
#                         min_m = s
#                         min_id = self.neighbors[i][j]
#                 id_min[i][k] = min_id
#                 all_d[i][k] = min_m
#         # 计算每个目标函数差值的平均值
#         means = []
#         for k in range(self.problem.n_obj):
#             for i in range(self.pop_size):
#                 n = 1
#                 s = 0
#                 if id_min[i][k] != -1:
#                     s += all_d[i][k]
#                     n += 1
#             m = s / n
#             means.append(m)
#         # 保存将要添加权重的个体对
#         add_id = []
#         # 判断个体是否处于稀疏的地区
#         for i in range(self.pop_size):
#             for k in range(self.problem.n_obj):
#                 if id_min[i][k] != -1 and all_d[i][k] > means[k]:
#                     if [i, id_min[i][k]] not in add_id and [id_min[i][k], i]:
#                         add_id.append([i, id_min[i][k]])
#                     break
#         for i in range(len(add_id)):
#             if self._is_have_weight(add_id[i][0], add_id[i][1]):
#                 w, p = self._get_weight(self.pop[add_id[i][0]], self.pop[add_id[i][1]])
#                 new_w.append(w)
#                 new_p.append(p)
#         return new_w, new_p
#
#     # 根据最近权重向量计算最近距离
#     def _dis_neighbor(self):
#         min_neighgbor = np.argsort(cdist(self.ref_dirs, self.ref_dirs), axis=1, kind='quicksort')[:, 1:self.problem.n_obj+1]
#         min_d = []
#         for i in range(self.pop_size):
#             d = []
#             for j in range(self.problem.n_obj):
#                 dj = self._distance(self.pop[i].get("F"), self.pop[min_neighgbor[i][j]].get("F"))
#                 d.append(dj)
#             min_d.append(d)
#         means = []
#         for k in range(self.problem.n_obj):
#             s = 0
#             for i in range(self.pop_size):
#                 s += min_d[i][k]
#             m = s / self.pop_size
#             means.append(m)
#        # 保存将要添加权重的个体对
#         add_id = []
#         # 判断个体是否处于稀疏的地区
#         for i in range(self.pop_size):
#             for k in range(self.problem.n_obj):
#                 if min_d[i][k] > means[k]:
#                     if [i, min_neighgbor[i][k]] not in add_id and [min_neighgbor[i][k], i] not in add_id:
#                         add_id.append([i, min_neighgbor[i][k]])
#                     break
#         new_w = []
#         new_p = []
#         for i in range(len(add_id)):
#             if self._is_have_weight(add_id[i][0], add_id[i][1]):
#                 w, p = self._get_weight(self.pop[add_id[i][0]], self.pop[add_id[i][1]])
#                 new_w.append(w)
#                 new_p.append(p)
#         return new_w, new_p
#
#     # 判断两个个体间是否存在权重向量
#     def _is_have_weight(self, i1, i2):
#         # 计算两个个体关联的权重向量间的距离r
#         r = self._distance(self.ref_dirs[i1], self.ref_dirs[i2])
#         d1 = []
#         d2 = []
#         # 计算两个权重向量的邻居中小于r的权重向量
#         for i in range(self.n_neighbors):
#             if self.neighbors[i1][i] != i1:
#                 d = self._distance(self.ref_dirs[i1], self.ref_dirs[self.neighbors[i1][i]])
#                 if d < r:
#                     d1.append(self.neighbors[i1][i])
#             if self.neighbors[i1][i] != i2:
#                 d = self._distance(self.ref_dirs[i2], self.ref_dirs[self.neighbors[i2][i]])
#                 if d < r:
#                     d2.append(self.neighbors[i2][i])
#         # 判断r内的邻居是否有交集，有则代表两个个体间存在权重向量
#         for i in d1:
#             if i in d2:
#                 return False
#         return True
#
#     # 增加权重向量及个体
#     def _get_weight(self, p1, p2):
#         # 获取两个个体的目标函数值
#         v1 = p1.get("F")
#         v2 = p2.get("F")
#         # 计算新的权重向量
#         new_ref = self._create_weight(v1, v2)
#         # 计算聚合函数值
#         FV1 = self._decomposition.do(v1, weights=new_ref, ideal_point=self.ideal_point)
#         FV2 = self._decomposition.do(v2, weights=new_ref, ideal_point=self.ideal_point)
#         # 以聚合函数值小的个体作为新个体
#         if FV1 < FV2:
#             return new_ref, p1
#         else:
#             return new_ref, p2
#
#     # 根据个体生成权重向量
#     def _create_weight(self, v1, v2):
#         length = len(v1)
#         FV = [(v1[i]+v2[i]) / 2  for i in range(length)]
#         sum_f = sum([FV[j]-self.ideal_point[j] for j in range(length)])
#         new_ref = [(FV[j] - self.ideal_point[j]) / sum_f for j in range(length)]
#         return new_ref
#
#     # 删除权重
#     def _delete_weight(self, new_w, new_p):
#         pop = []
#         for i in self.pop:
#             pop.append(i)
#         number = len(new_w)
#         for i in range(number):
#             max_c = 1
#             max_id = -1
#             for j in range(self.pop_size):
#                 # 找到最多相同的个体
#                 if pop.count(pop[j]) > max_c:
#                     max_c = pop.count(pop[j])
#                     max_id = j
#             # 如果没有相同的个体，则删除拥挤度最大的个体
#             if max_c == 1:
#                 # 计算拥挤度
#                 r = self._radius_R(r)
#                 crowd = self._crowd(self.pop)
#                 id = crowd.index(max(crowd))
#                 # 将拥挤度最大的个体删除，添加外部个体。
#                 self.pop[id] = new_p[i]
#                 self.ref_dirs[id] = new_w[i]
#             # 如果有多个重复的个体，删除其中聚合函数值最多的个体和权重向量
#             else:
#                 max_g = 0
#                 max_gid = -1
#                 id = -1
#                 # 遍历相同的个体
#                 for j in range(max_c):
#                     # 找下一个相同个体的id
#                     id = pop.index(pop[max_id], id+1)
#                     # 计算聚合函数值
#                     FV = self._decomposition.do(pop[id].get("F"), weights=self.ref_dirs[id], ideal_point=self.ideal_point)
#                     # 找到聚合函数最大的id
#                     if FV > max_g:
#                         max_g = FV
#                         max_gid = id
#                 # 删除个体、权重，添加新个体与新权重
#                 self.pop[max_gid] = new_p[i]
#                 self.ref_dirs[max_gid] = new_w[i]
#         # get_visualization("scatter").add(self.ref_dirs).show()
#         return
#
#
#     # 根据密度替换个体
#     def _exchange_id(self, D, A):
#          # 记录替换的id-密度较大的个体中，与权重向量角度最大的个体的id
#             max_angle_id = -1
#             max_angle = 0
#             # 记录已经比较过的id
#             ids = []
#             # 找前K个密度最大的个体
#             for i in range(self.n_neighbors):
#                 # 记录密度最大的id和值
#                 max_d = 0
#                 max_id = -1
#                 for j in range(self.pop_size):
#                     if D[j] > max_d and j not in ids:
#                         max_d = D[j]
#                         max_id = j
#                 ids.append(max_id)
#                 # 找出K个个体中与自己权重向量夹角最大的个体
#                 if A[max_id] > max_angle:
#                     max_angle = A[max_id]
#                     max_angle_id = max_id
#             ex_id = -1
#             min_angle = 180
#             # 找出外部种群中与要被替换的个体的权重向量夹角最小的个体
#             for j in range(len(self.external_pop)):
#                 a = self._angle(self.external_pop[j].get("F"), self.ref_dirs[max_angle_id])
#                 if a < min_angle:
#                     ex_id = j
#                     min_angle = a
#                 # if a == min_angle:
#                 #     # 如果角度相等，比较聚合函数值
#                 #     v1 = 1
#             return  max_angle_id, ex_id
#
#     # 求每个个体与其相对应的权重向量的角度
#     def _angle_A(self):
#         # 保存角度
#         A = []
#         for i in range(self.pop_size):
#             A.append(self._angle(self.pop[i].get("F"), self.ref_dirs[i]))
#         return A
#
#     # 求个体的拥挤度
#     # 输入：小生境半径
#     def _crowd(self, r):
#         D = [1 for i in range(self.pop_size)]
#         for i in range(self.pop_size):
#             for j in range(1, self.pop_size):
#                 d = self._distance(self.pop[i].get("F"), self.pop[j].get("F"))
#                 if d < r:
#                     D[i] *= d / r
#                     D[j] *= d / r
#         for i in range(self.pop_size):
#             D[i] = 1 - D[i]
#         return D
#
#
#     # 计算半径
#     # 输入：种群
#     # 输出：集合中所有解到其第k个最近解的距离的中值
#     def _radius_R(self):
#         # 保存个体间的欧式距离
#         dis = []
#         # 计算距离
#         for i in range(self.pop_size):
#             for j in range(self.n_neighbors):
#                 d = self._distance(self.pop[i].get("F"), self.pop[self.neighbors[i][j]].get("F"))
#                 dis.append(d)
#         # 以距离中的平均值作为小生境半径
#         r = np.median(dis)
#         return r
#
#
#     # 求两点间的欧几里得距离
#     def _distance(self, point1, point2):
#         s = sum([(point1[i] - point2[i])**2 for i in range(len(point1))])
#         return math.sqrt(s)
#
#
#      # 调整权重向量
#     def _weigt_adjust(self, A):
#         # get_visualization("scatter").add(self.ref_dirs).show()
#         # get_visualization("scatter").add(self.pop.get("F")).show()
#         for i in range(self.pop_size):
#             if A[i] > self.max_angle:
#                 sum_f = sum([self.pop[i].get("F")[j]-self.ideal_point[j] for j in range(self.problem.n_obj)])
#                 new_ref = [(self.pop[i].get("F")[j] - self.ideal_point[j]) / sum_f for j in range(self.problem.n_obj)]
#                 for j in range(self.problem.n_obj):
#                     self.ref_dirs[i][j] = new_ref[j]
#
#
#     # 计算两点与参考点形成的两个向量间的角度
#     # 输入：两点
#     # 输出：角度
#     def _angle(self, v1, v2):
#         dis1 = math.sqrt(sum([(v1[i]-self.ideal_point[i])**2 for i in range(self.problem.n_obj)]))
#         dis2 = math.sqrt(sum([(v2[i]-self.ideal_point[i])**2 for i in range(self.problem.n_obj)]))
#         s = sum([(v1[i]-self.ideal_point[i]) * (v2[i]-self.ideal_point[i]) for i in range(self.problem.n_obj)])
#         angle = math.degrees(math.acos(s / (dis1 * dis2)))
#         return angle
#
#
#     # 维护外部种群
#     # 输入：淘汰的个体
#     def _update(self, p):
#         i = 0
#         while i < len(self.external_pop):
#             flag = 0
#             for n in range(self.problem.n_obj):
#                 if p.get("F")[n] > self.external_pop[i].get("F")[n]:
#                     flag += 1
#             if flag == 0:
#                 self.external_pop.pop(i)
#             elif flag == self.problem.n_obj:
#                 return
#             else:
#                 i += 1
#         self.external_pop.append(p)
#
#     # 计算个体与权重向量的聚合函数值，根据聚合函数值选择个体
#     # 输入：父代+子代的种群
#     def _s_value(self, pop):
#         for i in range(self.pop_size):
#             min_value = 100000
#             min_id = -1
#             for j in range(len(pop)):
#                 v = self._decomposition.do(pop[j].get("F"), weights=self.ref_dirs[i], ideal_point=self.ideal_point)
#                 if v < min_value:
#                     min_id = j
#                     min_value = v
#             self.pop[i] = pop[j]
#             pop.pop(j)
#
#
#     # 支配关系
#     # 输入：种群
#     # 输出：种群的帕累托最优集，被支配的个体集合
#     def _level(self, pop):
#         # 将种群转化为list，方便后续操作
#         p = []
#         for i in pop:
#             p.append(i)
#         # 保存帕累托最优集合
#         P = []
#         # 保存被支配的个体
#         D = []
#         i = 0
#         j = 1
#         while i < len(p):
#             flag = 0
#             if j <len(p):
#                 for n in range(self.problem.n_obj):
#                     if p[i].get("F")[n] > p[j].get("F")[n]:
#                         flag += 1
#             # i和其他个体对比结束，没有支配i的个体。将i放到P中，并从p中删除
#             else:
#                 P.append(p[i])
#                 p.pop(i)
#                 j = 1
#                 continue
#             # 个体i被个体j支配，将i放到D中，并从p中删除
#             if flag == self.problem.n_obj:
#                 D.append(pop[i])
#                 p.pop(i)
#                 j = 1
#                 continue
#             # 个体j被个体i支配，将j放到D中，并从p中删除
#             elif flag == 0:
#                 D.append(pop[j])
#                 p.pop(j)
#             else:
#                 j += 1
#         D1 = Population(n_individuals=len(D))
#         for j in range(len(D)):
#             D1[i] = D[i]
#         return P, D1
#
#
#     # 根据密度选择个体
#     # 输入：种群，最后一级的个体
#     def _density_sel(self, pop, new_p):
#         # 得到需要选择的个体个数
#         l = self.pop_size - len(pop)
#         r = self._radius(pop)
#         # 计算new_p中每个个体的密度
#         den = []
#         for i in new_p:
#             den.append(self._density_p(i, pop, r))
#         for j in range(l):
#             id = den.index(min(den))
#             pop.append(new_p[id])
#             new_p.pop(id)
#             den.pop(id)
#         return pop
#
#
#     # 计算半径
#     # 输入：种群
#     # 输出：小生境半径(所有个体间的距离的平均值)
#     def _radius(self, pop):
#         # 保存个体间的欧式距离
#         dis = []
#         # 计算距离
#         for i in range(len(pop)):
#             for j in range(len(pop)):
#                 d = self._distance(pop[i].get("F"), pop[j].get("F"))
#                 dis.append(d)
#         # 以距离中的平均值作为小生境半径
#         r = sum(dis) / (len(dis) - len(pop))
#         return r
#
#
#     # 计算个体p在种群中的密度
#     # 输入：个体p,种群pop
#     def _density_p(self, p, pop, r):
#         den = 0
#         for i in range(len(pop)):
#             d = self._distance(pop[i].get("F"), p.get("F"))
#             if d <= r:
#                 den += 1
#         return den
#
#     # 计算输入种群中的个体的密度
#     # 输入：种群
#     # 输出：种群中个体间的密度
#     def _density(self, pop):
#         # 保存个体的密度
#         den = []
#         # 保存个体间的欧式距离
#         dis = []
#         # 计算距离
#         for i in range(len(pop)):
#             for j in range(len(pop)):
#                 d = self._distance(pop[i].get("F"), pop[j].get("F"))
#                 dis.append(d)
#         # 以距离中的平均值作为小生境半径
#         r = sum(dis) / (len(dis) - len(pop))
#         for i in range(len(dis)):
#             if dis[i] <= r:
#                 # dis[i] = dis[i] / r
#                 dis[i] = 1
#             else:
#                 dis[i] = 0
#         D = 0
#         for i in range(len(pop)):
#             for j in range(len(pop)):
#                 if i != j:
#                     D += dis[i*len(pop) + j]
#             den.append(D)
#         return den
#
#     # 获得交配池
#     # 输入：种群
#     # 输出：交配池
#     def _mating_pool(self, Q):
#         den = self._density(Q)
#         for i in range(len(den)):
#             for n in range(int(len(Q) / (den[i]+1))-1):
#                 Q = Q.merge(self.pop[i])
#         return Q
#
#     #-------------------------------------------------------------------------------------------------------------------
#
#     def _dbs(self, Z, L):
#         for i in range(L):
#             Q[i] = []
#         for j in range(self.pop_size):
#             P = self._associate()
#         return
#
#     def _associate(self):
#         return
