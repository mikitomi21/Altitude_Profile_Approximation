def lagrange(x_values, y_values, x):
    result = 0.0
    for i in range(len(x_values)):
        term = y_values[i]
        for j in range(len(x_values)):
            if i != j:
                term *= (x - x_values[j]) / (x_values[i] - x_values[j])
        result += term
    return result

def solve_tridiagonal_system(a, b, c, d):
    """
    Funkcja rozwiązująca układ trójdiagonalny Ax = d, gdzie:
    - a, b, c: listy odpowiadające diagonalom macierzy A
    - d: lista prawych stron

    Zwraca:
    - lista rozwiązań x
    """

    n = len(d)
    c_dash = [0.0] * n
    x = [0.0] * n

    c_dash[0] = c[0] / b[0]
    for i in range(1, n):
        c_dash[i] = c[i] / (b[i] - a[i] * c_dash[i - 1])

    d_dash = [0.0] * n
    d_dash[0] = d[0] / b[0]
    for i in range(1, n):
        d_dash[i] = (d[i] - a[i] * d_dash[i - 1]) / (b[i] - a[i] * c_dash[i - 1])

    x[n - 1] = d_dash[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = d_dash[i] - c_dash[i] * x[i + 1]

    return x


def spline(x_values, y_values, x):
    n = len(x_values)  # liczba punktów danych
    h = [x_values[i + 1] - x_values[i] for i in range(n - 1)]  # różnice między kolejnymi wartościami x
    alpha = [0.0] * n
    for i in range(1, n - 1):
        alpha[i] = (3 / h[i]) * (y_values[i + 1] - y_values[i]) - (3 / h[i - 1]) * (y_values[i] - y_values[i - 1])

    l = [1.0] * n
    mu = [0.0] * n
    z = [0.0] * n
    c = [0.0] * n
    b = [0.0] * n
    d = [0.0] * n

    for i in range(1, n - 1):
        l[i] = 2 * (x_values[i + 1] - x_values[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

    l[n - 1] = 1.0
    z[n - 1] = 0.0
    c[n - 1] = 0.0

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