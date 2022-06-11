import argparse
import sys
import math
import random

def inverse(a,n):
    return pow(a,-1,n)

class Point:
    """
    Point:
        Coordinates for point p=(x,y) on the plane.
    """
    def __init__(self, x: int, y: int) -> None:
        self.x = x
        self.y = y

    def __eq__(self, other) -> bool:
        return type(self) == type(other) and (self.x == other.x) and (self.y == other.y)

    def __repr__(self):
        return f"(x,y) = ({self.x}, {self.y})"

class CannotAdd(Exception):
    pass

class EllipticCurve:
    '''
    Elliptic curve:
        Takes the parameters a,b to describe Weierstrass equation y^2 = x^3 + ax + b over Z/NZ.
        Implements explicit sum of points on the elliptic curve itself and fast multiplication.
    '''
    def __init__(self, a: int, b: int, N: int) -> None:
        self.a = a
        self.b = b
        self.N = N
        self.div = 1

    def sum(self, p: Point, q: Point) -> Point:
        '''
        Calculates the sum of two points p, q in the elliptic curve.
        '''
        try:
            if type(p) == type(None):
                return q 
            if type(q) == type(None):
                return p
            if p == q:
                # Raise exception if division by 0 in Z/NZ
                if math.gcd(2*p.y, self.N) != 1:
                    self.div = math.gcd(2*p.y, self.N)
                    raise CannotAdd
                
                slope = (3*p.x*p.x + self.a) * inverse(2*p.y, self.N)
                x = (slope*slope - 2*p.x) % self.N
                y = (slope*(p.x-x) - p.y) % self.N
                return Point(x,y)
            
            else:
                # Raise exception if division by 0 in Z/NZ
                if math.gcd(q.x-p.x, self.N) != 1:
                    self.div = math.gcd(q.x-p.x, self.N)
                    raise CannotAdd

                slope = (q.y-p.y) * inverse(q.x-p.x, self.N)
                x3 = (slope*slope - p.x - q.x) % self.N
                y3 = (slope * (p.x-x3) - p.y) % self.N
                return Point(x3,y3)
        except CannotAdd:
            return None

    def mult(self, p: Point, k: int) -> Point:
        '''
        Fast calculation of k*p by doing 2^b1 p + 2^b2 p + ... looking at the binary expansion of k instead of adding p successively k times. 
        Reduces complexity from O(k) to O(log(k)).
        '''
        if p == None:
            return None

        binary_expansion = bin(k)[2:]
        m = len(binary_expansion)
        
        try:
            powers_of_two = [p]
            for _ in range(m-1):
                new_point = self.sum(powers_of_two[-1], powers_of_two[-1])
                powers_of_two.append(new_point)

            q = None
            for i, v in enumerate(binary_expansion):
                if v == "1":
                    q = self.sum(powers_of_two[i], q)
            return q
        
        except CannotAdd:
            return None
    
    def check(self, p: Point, q: Point) -> int:
        '''
        Returns a common factor between N and slope coordinates.
        '''
        if p == q:
            return math.gcd(2*p.y, self.N)
        else:
            return math.gcd(q.x-p.x, self.N)

    
def Lenstra(N):
    x0, y0, a = [random.randint(1,N) for _ in range(3)]
    b = (y0*y0 - x0*x0*x0 - a*x0) % N
    E = EllipticCurve(a,b,N)
    point = Point(x0,y0)
    i = 1
    while E.div == 1:
        try:
            if point == None:
                break
            point = E.mult(point,i)
            i += 1
        except CannotAdd:
            break
    return E.div

if __name__ == "__main__":
    '''
    Factor numbers of the form N = pq with p,q large primes using Lenstra's algorithm.
    If N is not of that form, the program returns just a divisor (not necessarily prime) of N (including maybe N itself).
    '''
    parser = argparse.ArgumentParser(description="Factor numbers of the form N = pq with p,q large primes using Lenstra's algorithm.")
    parser.add_argument('N', metavar='N', type=int, help='the number to factor')
    args = parser.parse_args()
    N = args.N
    print(Lenstra(N))