import numpy as np

def lagrange(x_values, y_values, x):
    result = 0.0
    for i in range(len(x_values)):
        term = y_values[i]
        for j in range(len(x_values)):
            if i != j:
                term *= (x - x_values[j]) / (x_values[i] - x_values[j])
        result += term
    return result

def spline(x_values, y_values, x):
    n = len(x_values)  # liczba punktów danych
    h = [x_values[i + 1] - x_values[i] for i in range(n - 1)]  # różnice między kolejnymi wartościami x
    alpha = [0.0] * n
    for i in range(1, n - 1):
        alpha[i] = (3 / h[i]) * (y_values[i + 1] - y_values[i]) - (3 / h[i - 1]) * (y_values[i] - y_values[i - 1])

    l = np.ones(n)
    mu = np.zeros(n)
    z = np.zeros(n)
    c = np.zeros(n)
    b = np.zeros(n)
    d = np.zeros(n)

    xi = 0.0
    yi = 0.0
    bi = 0.0
    ci = 0.0

    # Rozwiązujemy układ równań za pomocą trójdiagonalnej eliminacji Gaussa
    for i in range(1, n - 1):
        l[i] = 2 * (x_values[i + 1] - x_values[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

    l[n - 1] = 1.0
    z[n - 1] = 0.0
    c[n - 1] = 0.0

    # Na podstawie wyżej wymienionych wartości wyzaczamy b,c,d
    for j in range(n - 2, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = (y_values[j + 1] - y_values[j]) / h[j] - (h[j] / 3) * (c[j + 1] + 2 * c[j])
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])

    # Znajdowanie odpowiedniego przedziału
    for i in range(n - 1):
        if x_values[i] <= x <= x_values[i + 1]:
            xi = x_values[i]
            yi = y_values[i]
            hi = h[i]
            bi = b[i]
            ci = c[i]
            di = d[i]
            break

    # Obliczanie wartości interpolowanej dla x
    dx = x - xi
    interpolated_y = yi + bi * dx + ci * dx**2 + di * dx**3

    return interpolated_y