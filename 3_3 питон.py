import math

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# коэффициенты полиномов эйлера (для a)
euler = [
    [],  # 0
    [1],
    [1, 1],
    [1, 4, 1],
    [1, 11, 11, 1],
    [1, 26, 66, 26, 1],
    [1, 57, 302, 302, 57, 1],
    [1, 120, 1191, 2416, 1191, 120, 1],
    [1, 247, 4293, 15619, 15619, 4293, 247, 1],
    [1, 502, 14608, 88234, 156190, 88234, 14608, 502, 1],
    [1, 1013, 46834, 472063, 1561900, 1561900, 472063, 46834, 1013, 1]  # 10
]

def rational_sum(a, b):
    if a < 1 or a > 10:
        return ""
    
    coeffs = euler[a]
    deg = len(coeffs) - 1
    
    # вычисление числителя по схеме Горнера
    num = 0
    for i in range(deg, -1, -1):
        num = num * b + coeffs[i]
    num *= b
    
    # знаменатель (b-1)^(a+1)
    den = 1
    for i in range(a + 1):
        den *= (b - 1)
    
    g = gcd(num, den)
    num //= g
    den //= g
    
    return f"{num}/{den}"

def main():
    a, b = map(int, input().split())
    
    if b == 1:
        print("infinity")
        return
    
    res = rational_sum(a, b)
    if not res:
        print("irrational")
    else:
        print(res)

if __name__ == "__main__":
    main()