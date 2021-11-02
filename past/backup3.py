# # 新权重以及个体替换旧个体
#     def _adjust(self, p_id, new_w, new_p):
#         # get_visualization("scatter").add(self.ref_dirs).show()
#         # get_visualization("scatter").add(self.pop.get("F")).show()
#         for i in range(len(p_id)):
#             self.pop[p_id[i]] = new_p[i]
#             self.ref_dirs[p_id[i]] = new_w[i]
#             self.history_flag[i] = 0
#             self.last_angle = self._angle(new_w[i], new_p[i].get("F"))
#         plot = Scatter()
#         plot.add(self.pop.get("F"), color="red")
#         plot.add(self.ref_dirs, color="blue")
#         plot.show()
#         # get_visualization("scatter").add(self.ref_dirs).show()
#         # get_visualization("scatter").add(self.pop.get("F")).show()
#
#
#     # 需要被调整的权重向量
#     def _weigt_adjust(self, A):
#         p_id = []
#         for i in range(self.pop_size):
#             if A[i] > self.max_angle and self.history_flag[i] >= 20:
#                 self.history_flag[i] = 0
#                 p_id.append(i)
#         if len(p_id) > 0:
#             p = [self.pop[i].get("F") for i in p_id]
#             re = [self.ref_dirs[i] for i in p_id]
#             plot = Scatter()
#             plot.add(self.pop.get("F"), color="red")
#             plot.add(self.ref_dirs, color="blue")
#             plot.show()
#             plot = Scatter()
#             plot.add(np.array(p), color="green")
#             plot.add(np.array(re), color="yellow")
#             plot.show()
#         return p_id
#
#     # 求每个个体与其相对应的权重向量的角度
#     def _angle_A(self):
#         # 保存角度
#         A = []
#         for i in range(self.pop_size):
#             A.append(self._angle(self.pop[i].get("F"), self.ref_dirs[i]))
#         return A
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
#     # 获取个体的最近邻居的差值
#     def _get_add_area(self):
#         min_neighgbor = [[self.neighbors[i][j+1] for j in range(self.problem.n_obj)] for i in range(self.pop_size)]
#         min_d = [[self._distance(self.pop[i].get("F"), self.pop[min_neighgbor[i][j]].get("F")) for j in range(self.problem.n_obj)]
#                  for i in range(self.pop_size)]
#         d = []
#         for i in range(self.pop_size):
#             for j in range(self.problem.n_obj):
#                 d.append(min_d[i][j])
#         means = np.mean(d)
#         # 判断个体是否处于稀疏的地区
#         for i in range(self.pop_size):
#             for k in range(self.problem.n_obj):
#                 if min_d[i][k] > means:
#                     self.flag_dis[i][k] += 1
#                 else:
#                     self.flag_dis[i][k] = 0
#
#     # 根据最近权重向量计算最近距离
#     def _dis_neighbor(self, length):
#         min_neighgbor = [[self.neighbors[i][j+1] for j in range(self.problem.n_obj)] for i in range(self.pop_size)]
#         add_id = []
#         for i in range(self.pop_size):
#             for k in range(self.problem.n_obj):
#                 if self.flag_dis[i][k] >= 10:
#                     if [i, min_neighgbor[i][k]] not in add_id and [min_neighgbor[i][k], i] not in add_id:
#                         add_id.append([i, min_neighgbor[i][k]])
#         new_w = []
#         new_p = []
#         a = []
#         for i in range(len(add_id)):
#             if self._is_have_weight(add_id[i][0], add_id[i][1]):
#                 a.append(add_id[i])
#         add_id = a
#         for i in range(len(add_id)):
#             w, p = self._get_weight(add_id[i][0], add_id[i][1])
#             new_w.append(w)
#             new_p.append(p)
#         if len(add_id) > length:
#             r = self._radius_R()
#             D = [self._p_crowd(r, new_p[i]) for i in range(len(add_id))]
#             max_list = np.argsort(D, kind='quicksort')[:(len(add_id)-length)]
#             w = []
#             p = []
#             for i in range(len(new_w)):
#                 if i not in max_list:
#                     w.append(new_w[i])
#                     p.append(new_p[i])
#             new_w = w
#             new_p = p
#             print("s1")
#         elif len(add_id) < length:
#             r = self._radius_R()
#             D = self._crowd(r)
#             min_list = np.argsort(D, kind='quicksort')[:(length - len(add_id))]
#             for i in min_list:
#                 sum_f = sum([self.pop[i].get("F")[j]-self.ideal_point[j] for j in range(self.problem.n_obj)])
#                 new_ref = [(self.pop[i].get("F")[j] - self.ideal_point[j]) / sum_f for j in range(self.problem.n_obj)]
#                 new_w.append(new_ref)
#                 new_p.append(self.pop[i])
#             print("s2")
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
#     def _get_weight(self, id1, id2):
#         # 获取两个个体的目标函数值
#         v1 = self.ref_dirs[id1]
#         v2 = self.ref_dirs[id2]
#         # 计算新的权重向量
#         new_ref = self._create_weight(v1, v2)
#         # 计算聚合函数值
#         FV1 = self._decomposition.do(self.pop[id1].get("F"), weights=new_ref, ideal_point=self.ideal_point)
#         FV2 = self._decomposition.do(self.pop[id2].get("F"), weights=new_ref, ideal_point=self.ideal_point)
#         # 以聚合函数值小的个体作为新个体
#         if FV1 < FV2:
#             return new_ref, self.pop[id1]
#         else:
#             return new_ref, self.pop[id2]
#
#     # 根据个体生成权重向量
#     def _create_weight(self, v1, v2):
#         new_ref = [(v1[i] + v2[i]) / 2 for i in range(len(v1))]
#         return new_ref
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
#     # 计算单个个体的拥挤度
#     def _p_crowd(self, r, p):
#         D = 1
#         for i in range(self.pop_size):
#             d = self._distance(self.pop[i].get("F"), p.get("F"))
#             if d < r:
#                 D *= d / r
#         D = 1 - D
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
