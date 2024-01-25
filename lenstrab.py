import random
import math

# Define the Elliptic Curve
class EllipticCurve:
    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.p = p

    def is_on_curve(self, x, y):
        return (y * y) % self.p == (x * x * x + self.a * x + self.b) % self.p

class Point:
    def __init__(self, curve, x, y):
        self.curve = curve
        self.x = x
        self.y = y

        if not curve.is_on_curve(x, y):
            raise ValueError("The point is not on the curve.")

def point_addition(p1, p2):
    if p1.x == p2.x and p1.y == -p2.y:
        return None

    if p1 == p2:
        lam = (3 * p1.x * p1.x + p1.curve.a) * pow(2 * p1.y, -1, p1.curve.p)
    else:
        lam = (p2.y - p1.y) * pow(p2.x - p1.x, -1, p1.curve.p)

    x3 = (lam * lam - p1.x - p2.x) % p1.curve.p
    y3 = (lam * (p1.x - x3) - p1.y) % p1.curve.p

    return Point(p1.curve, x3, y3)

def point_multiplication(p, n):
    result = p
    for _ in range(1, n):
        result = point_addition(result, p)
        if result is None:
            break
    return result

def lenstras_algorithm(N, max_iterations=1000):
    for _ in range(max_iterations):
        a = random.randint(0, N-1)
        x = random.randint(0, N-1)
        y = random.randint(0, N-1)
        b = (y*y - x*x*x - a*x) % N

        curve = EllipticCurve(a, b, N)
        point = Point(curve, x, y)

        for q in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
            try:
                point = point_multiplication(point, q)
            except ValueError:
                p = math.gcd(curve.p, 2 * point.y)
                if p > 1 and p < N:
                    return p

    return None

# Example usage
N = 8051  # Change this to the number you want to factor
factor = lenstras_algorithm(N)
print(f"A non-trivial factor of {N} is: {factor}")
