import numpy as np

def sort_test():
    arr = [1, 7, 4, 12, 52, 2, 6]
    arr = np.array(arr)
    index = np.argsort(arr)
    print(index)


def random_test():
    arr = [1, 3, 5, 6, 7]
    n = np.random.choice(arr)
    print(n)
# random_test()


def mat_test():
    a = np.array([1, 4, 6])
    b = a.T * a
    print(b)
    return b

def arr_add():
    a = np.array([[1, 4, 6], [2, 1, 5]])
    b = np.array([0, 3, 7])
    print(a)
    np.append(a, [b], axis=0)
    print(a)


a = [1, 2]
b = a
b[0] = 4
print(a)
sort_test()
