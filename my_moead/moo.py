from pymoo.algorithms.moo.moead import MOEAD
from pymoo.algorithms.moo.nsga2 import NSGA2
# from pymoo.algorithms.my_moead import M_MOEAD
# from pymoo.algorithms.my_moead_n import M_MOEAD
# from pymoo.algorithms.moead_adaw import MOEAD_AdaW

# from pymoo.algorithms.moead_adaw_n import MOEAD_AdaW_N
# from pymoo.algorithms.moead_awa import MOEAD_AWA
from pymoo.factory import get_problem, get_visualization, get_reference_directions, get_performance_indicator
from pymoo.optimize import minimize
import math


def read_data(name):
    path = "./data/"+name+".txt"
    f = open(path, 'r')                                           # 打开文件
    data = []
    for line in f:                                                  # 按行读取数据
        l = line.split()                                            # 将一行的数据按空格划分
        l = [float(l[i]) for i in range(len(l))]
        data.append(l)
    f.close()                                                       # 关闭文件
    return data


# MOEA/D
def MOEAD_test(problem, n_gen, n_problem, pop_size=300):
    if n_problem == 2:
        algorithm = MOEAD(
            get_reference_directions("uniform", n_problem, n_partitions=pop_size),
            n_neighbors=15,
            decomposition="pbi",
            prob_neighbor_mating=0.7,
            seed=1,
            verbose=True
        )
        pf = problem.pareto_front()
    else:
        ref_dirs = get_reference_directions("uniform", n_problem, n_partitions=40)
        pf = problem.pareto_front(ref_dirs)
        algorithm = MOEAD(
            get_reference_directions("uniform", n_problem, n_partitions=30),
            n_neighbors=30,
            decomposition="pbi",
            prob_neighbor_mating=0.7,
            seed=1,
            verbose=True
        )

        # get_visualization("scatter", angle=(45,45)).add(pf).show()

    res = minimize(problem, algorithm, termination=('n_gen', n_gen))

    get_visualization("scatter").add(res.F).show()
    igd = get_performance_indicator("igd", pf)
    IGD = igd.do(res.F)
    print("IGD:", IGD)
    return IGD


# def MOEADT_test(problem, n_gen, n_problem, pop_size=199):
#     if n_problem == 2:
#         algorithm = MOEAD_T(
#             get_reference_directions("uniform", n_problem, n_partitions=pop_size),
#             n_neighbors=15,
#             decomposition="pbi",
#             prob_neighbor_mating=0.7,
#             seed=1,
#             verbose=True
#         )
#         pf = problem.pareto_front()
#     else:
#         ref_dirs = get_reference_directions("uniform", n_problem, n_partitions=40)
#         pf = problem.pareto_front(ref_dirs)
#         algorithm = MOEAD_T(
#             get_reference_directions("uniform", n_problem, n_partitions=30),
#             n_neighbors=30,
#             decomposition="pbi",
#             prob_neighbor_mating=0.7,
#             seed=1,
#             verbose=True
#         )
#
#         # get_visualization("scatter", angle=(45,45)).add(pf).show()
#
#     res = minimize(problem, algorithm, termination=('n_gen', 300))
#
#     # get_visualization("scatter").add(res.F).show()
#     igd = get_performance_indicator("igd", pf)
#     IGD = igd.calc(res.F)
#     print("IGD:", IGD)
#     return IGD
#
#
# def MOEADadaw_test(problem, n_gen, n_problem, pop_size=300):
#     if n_problem == 2:
#         algorithm = MOEAD_AdaW(
#             get_reference_directions("uniform", n_problem, n_partitions=pop_size),
#             n_neighbors=15,
#             decomposition="tchebi",
#             prob_neighbor_mating=0.7,
#             seed=1,
#             verbose=True
#         )
#         pf = problem.pareto_front()
#     else:
#         algorithm = MOEAD_AdaW(
#             get_reference_directions("uniform", n_problem, n_partitions=30),
#             n_neighbors=30,
#             decomposition="pbi",
#             prob_neighbor_mating=0.7,
#             seed=1,
#             verbose=True
#         )
#
#         ref_dirs = get_reference_directions("uniform", n_problem, n_partitions=40)
#         pf = problem.pareto_front(ref_dirs)
#
#     res = minimize(problem, algorithm, termination=('n_gen', n_gen))
#
#     get_visualization("scatter").add(res.F).show()
#     igd = get_performance_indicator("igd", pf)
#     IGD = igd.calc(res.F)
#     print("IGD:", IGD)
#     return IGD

#
#
# # MOEA/D-
# def MMOEAD_test(problem, n_gen, n_problem, pop_size=300):
#     if n_problem == 2:
#         algorithm = M_MOEAD(
#             get_reference_directions("uniform", n_problem, n_partitions=pop_size),
#             n_neighbors=15,
#             decomposition="pbi",
#             prob_neighbor_mating=0.7,
#             seed=1,
#             verbose=True
#         )
#         pf = problem.pareto_front()
#     else:
#         ref_dirs = get_reference_directions("uniform", n_problem, n_partitions=40)
#         pf = problem.pareto_front(ref_dirs)
#         algorithm = M_MOEAD(
#             get_reference_directions("uniform", n_problem, n_partitions=30),
#             n_neighbors=30,
#             decomposition="pbi",
#             prob_neighbor_mating=0.7,
#             seed=1,
#             verbose=True
#         )
#
#     res = minimize(problem, algorithm, termination=('n_gen', n_gen))
#
#     get_visualization("scatter").add(res.F).show()
#     igd = get_performance_indicator("igd", pf)
#     IGD = igd.calc(res.F)
#     print("IGD:", IGD)
#     return IGD
#
#

#
# def MOEADadawn_test(problem, n_gen, n_problem, pop_size=300):
#     if n_problem == 2:
#         algorithm = MOEAD_AdaW_N(
#             get_reference_directions("uniform", n_problem, n_partitions=pop_size),
#             n_neighbors=15,
#             decomposition="pbi",
#             prob_neighbor_mating=0.7,
#             seed=1,
#             verbose=True
#         )
#         pf = problem.pareto_front()
#     else:
#         algorithm = MOEAD_AdaW_N(
#             get_reference_directions("uniform", n_problem, n_partitions=30),
#             n_neighbors=30,
#             decomposition="pbi",
#             prob_neighbor_mating=0.7,
#             seed=1,
#             verbose=True
#         )
#
#         ref_dirs = get_reference_directions("uniform", n_problem, n_partitions=40)
#         pf = problem.pareto_front(ref_dirs)
#
#     res = minimize(problem, algorithm, termination=('n_gen', n_gen))
#
#     # get_visualization("scatter").add(res.F).show()
#     igd = get_performance_indicator("igd", pf)
#     IGD = igd.calc(res.F)
#     print("IGD:", IGD)
#     return IGD
#
#
# def MOEADawa_test(problem, n_gen, n_problem, pop_size=300):
#     if n_problem == 2:
#         algorithm = MOEAD_AWA(
#             get_reference_directions("uniform", n_problem, n_partitions=pop_size),
#             n_neighbors=15,
#             decomposition="pbi",
#             prob_neighbor_mating=0.7,
#             seed=1,
#             verbose=True
#         )
#         pf = problem.pareto_front()
#     else:
#         algorithm = MOEAD_AWA(
#             get_reference_directions("uniform", n_problem, n_partitions=30),
#             n_neighbors=30,
#             decomposition="pbi",
#             prob_neighbor_mating=0.7,
#             seed=1,
#             verbose=True
#         )
#
#         ref_dirs = get_reference_directions("uniform", n_problem, n_partitions=40)
#         pf = problem.pareto_front(ref_dirs)
#
#     res = minimize(problem, algorithm, termination=('n_gen', n_gen))
#
#     # get_visualization("scatter").add(res.F).show()
#     igd = get_performance_indicator("igd", pf)
#     IGD = igd.calc(res.F)
#     print("IGD:", IGD)
#     return IGD
#

# def NSGA3_test(problem, n_gen, n_problem=2, pop_size=300):
#     if n_problem == 2:
#         algorithm = NSGA3(
#             ref_dirs=get_reference_directions("uniform", n_problem, n_partitions=pop_size),
#             pop_size=301
#         )
#         pf = problem.pareto_front()
#     else:
#         algorithm = M_MOEAD(
#             ref_dirs=get_reference_directions("uniform", n_problem, n_partitions=30),
#             pop_size=496
#         )
#
#         ref_dirs = get_reference_directions("uniform", n_problem, n_partitions=40)
#         pf = problem.pareto_front(ref_dirs)
#
#     res = minimize(problem, algorithm, termination=('n_gen', n_gen))
#
#     # get_visualization("scatter").add(res.F).show()
#     igd = get_performance_indicator("igd", pf)
#     IGD = igd.calc(res.F)
#     print("IGD:", IGD)
#     return IGD


# NSGA-II
def NSGA2_test(problem, n_gen, pop_size):

    algorithm = NSGA2(pop_size=pop_size)

    res = minimize(problem,
                   algorithm,
                   ('n_gen', n_gen),
                   seed=1,
                   verbose=False)

    # plot = Scatter()
    # plot.add(problem.pareto_front(), plot_type="line", color="black", alpha=0.7)
    # plot.add(res.F, color="red")
    # plot.show()
    # d = read_data("ZDT3")
    igd = get_performance_indicator("igd", problem.pareto_front())
    IGD = igd.calc(res.F)
    print("IGD:", IGD)
    return IGD


if __name__ == "__main__":
    problem = get_problem("zdt3")
    stand_dev = 0
    igds = []
    s = 0
    t = 10
    for i in range(t):
        igd = MOEAD_test(problem, 200, 2)
        s += igd
        igds.append(igd)
    avg = s / t
    for i in range(t):
        stand_dev += (igds[i] - avg) ** 2
    stand_dev = math.sqrt(stand_dev / t)
    with open("D:/test.txt", "w") as f:
        f.write(str(avg))
        f.write("\n")
        f.write(str(stand_dev))
        f.close()
    print("avg:", avg)
    print("standard deviation:", stand_dev)

