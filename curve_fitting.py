import numpy as np
import matplotlib.pyplot as plt


def curve_fit(x_points, y_points, order):
    """
    :param x_points: x_coordinates of data points (in our frog example, the number of decaliters of insecticide to be sprayed in a given month)
    :param y_points: y_coordinates of data points (the freqeuncy-to-quality score which measures the average frog's sexual satisfaction in that month)
    :param order: order of the polynomial you wish to fit to the data points (e.g 1-linear, 2-quadratic, 3-cubic, etc)
    :return: an array of the coefficients of the polynomial function that fits the data points with the least error
    """
    N = len(x_points)
    y_sum = sum(y_points)
    x_sums = [0] * (order*2)
    xy_prods = [0] * order

    # create an array containing all the x_sums
    for n in range(1, order*2+1):
        x_sums[n-1] = sum([i**n for i in x_points])

    # creates an array containing all the sums of xy products
    for i in range(1, order+1):
        xy_prods[i-1] = 0
        for j in range(0, N):
            xy_prods[i-1] += (x_points[j]**i) * y_points[j]
    
    # Recalll that the aim is to solve the equation Ax = b
    
    # we create an empty matrix, A, for the x_points
    matrix = np.zeros(shape=(order+1, order+1))
    matrix[0][0] = N
    
    # now we assemble the rows of the matrix

    # row 0
    for i in range(1, order+1):
        matrix[0][i] = x_sums[i-1]

    # rows 1-j
    for r in range(1, order+1):
        row_r = []
        for v in range(r-1, order+r):
            row_r.append(x_sums[v])
        matrix[r] = row_r

    # now we make b
    b = np.zeros((order + 1))
    b[0] = y_sum

    for i in range(1, order+1):
        b[i] = xy_prods[i-1]

    # solve for the coefficients
    ans = np.matmul(np.linalg.inv(matrix), b.T)
    return ans


def f(coeffs, t):
    # rev_coeffs = np.flip(coefficients, axis=0)
    eq = 0
    for i in range(len(coeffs)):
        if i == 0:
            eq += coeffs[i]
        else:
            eq += coeffs[i] * (t**i)
    return eq


if __name__ == "__main__":
    # x = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5])
    # y = np.array([0, 1.5, 3.0, 4.5, 6.0, 7.5])

    # x = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5])
    # y = np.array([-.4326, -.1656, 3.1253, 4.7877, 4.8535, 8.6909])

    # x = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5])
    # y = np.array([0, .25, 1.0, 2.25, 4.0, 6.25])

    x = np.array([0, 95, 10, 1.5, 2.0, 2.5])
    y = np.array([10, -0.9156, 4.2513, 9.0377, 8.9757, 7.9409])

    coefficients = curve_fit(x, y, 1)
    print(coefficients)

    plt.figure(1)
    plt.subplot(211)
    plt.plot(x, f(coefficients, x))
    plt.scatter(x, y)
    plt.title("Order = 1")
    plt.show()








