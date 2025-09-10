import random
import math
from typing import List, Tuple

# решето эратосфена
def sieve_of_eratosthenes(limit: int) -> List[int]:
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False

    for p in range(2, int(math.sqrt(limit)) + 1):
        if is_prime[p]:
            for i in range(p * p, limit + 1, p):
                is_prime[i] = False

    primes = [i for i, prime in enumerate(is_prime) if prime]
    return primes

# возведение в степень по модулю
def mod_pow(base: int, exp: int, mod: int) -> int:
    result = 1
    base = base % mod

    while exp > 0:
        if exp % 2 == 1:
            result = (result * base) % mod
        exp = exp >> 1
        base = (base * base) % mod

    return result

# тест миллера
def miller_test(n: int, factorization: List[Tuple[int, int]], t: int) -> bool:
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False

    # выбор t различных случайных чисел a_j
    a_values = []
    for _ in range(t):
        while True:
            a = random.randint(2, n - 2)
            if a not in a_values:
                a_values.append(a)
                break

    # проверка a_j^(n-1) mod n == 1
    for a in a_values:
        if mod_pow(a, n - 1, n) != 1:
            return False

    # для каждого простого делителя q_i
    for q, _ in factorization:
        all_one = True
        for a in a_values:
            exponent = (n - 1) // q
            if mod_pow(a, exponent, n) == 1:
                all_one = False
                break
        if all_one:
            return False

    return True

# генерация простого числа (миллер)
def generate_prime_miller(bit_length: int, primes: List[int], t: int) -> int:
    # границы для m
    min_m = (1 << (bit_length - 2)) + 1
    max_m = (1 << (bit_length - 1)) - 1

    while True:
        # строим число m в нужном диапазоне
        while True:
            m = 1
            factorization = []

            while m < max_m:
                q = random.choice(primes)
                alpha = random.randint(1, 3)

                # проверка на переполнение
                if m > max_m // (q ** alpha):
                    break

                m *= q ** alpha
                factorization.append((q, alpha))

            if min_m <= m <= max_m:
                break

        # вычисляем n = 2m + 1
        n = 2 * m + 1

        # доп проверка битовой длины
        if n.bit_length() != bit_length:
            continue

        # проверка n тестом миллера
        if miller_test(n, factorization, t):
            return n

# тест поклингтона (проверка простоты)
def pocklington_test(n: int, F: int, F_factorization: List[Tuple[int, int]], t: int) -> bool:
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False

    R = (n - 1) // F
    if F <= math.isqrt(n) - 1:
        return False

    # выбор t различных случайных чисел a_j
    a_values = []
    for _ in range(t):
        while True:
            a = random.randint(2, n - 2)
            if a not in a_values:
                a_values.append(a)
                break

    # проверка a_j^(n-1) mod n == 1
    for a in a_values:
        if mod_pow(a, n - 1, n) != 1:
            return False

    # для каждого a_j проверить условия
    for a in a_values:
        condition_met = False
        for q, _ in F_factorization:
            exponent = (n - 1) // q
            if mod_pow(a, exponent, n) != 1:
                condition_met = True
                break

        if not condition_met:
            return False

    return True

# генерация простого числа (поклингтона)
def generate_prime_pocklington(bit_length: int, primes: List[int], t: int) -> int:
    # границы для F (на 1 бит больше половины требуемого размера)
    min_F = (1 << (bit_length // 2)) + 1
    max_F = (1 << (bit_length // 2 + 1)) - 1

    while True:
        # построение числа F
        F = 1
        F_factorization = []

        while F < max_F:
            q = random.choice(primes)
            alpha = random.randint(1, 3)

            # проверка на переполнение
            if F > max_F // (q ** alpha):
                break

            F *= q ** alpha
            F_factorization.append((q, alpha))

        if F < min_F:
            continue

        # выбор случайного четного R
        R = 2 * random.randint(1, (1 << (bit_length // 2 - 1)) - 1)
        n = R * F + 1

        # проверка размера
        if n.bit_length() != bit_length:
            continue

        # проверка n тестом поклингтона
        if pocklington_test(n, F, F_factorization, t):
            return n

# генерация простого числа (гост)
def generate_prime_gost(t: int, q: int) -> int:
    max_p = (1 << t) - 1
    attempts = 0
    max_attempts = 1000

    while attempts < max_attempts:
        attempts += 1
        xi = random.random()

        # вычисление N
        N = math.ceil((2 ** (t - 1)) / q) + math.ceil((2 ** (t - 1)) * xi / q)
        if N % 2 != 0:
            N += 1

        # инициализация u
        for u in range(0, 4 * (q + 1), 2):
            # шаг 3: вычислить p
            p = (N + u) * q + 1

            # проверка размера p
            if p > max_p:
                break

            # проверка теоремы диемитко при a = 2
            if mod_pow(2, p - 1, p) == 1 and mod_pow(2, N + u, p) != 1:
                return p

    return -1  # если не удалось найти простое число

# функция для подбора q
def find_suitable_q(bit_length: int, primes: List[int]) -> int:
    q_bits = (bit_length + 1) // 2
    min_q = (1 << (q_bits - 1)) + 1
    max_q = (1 << q_bits) - 1

    # ищем q в таблице
    for q in primes:
        if min_q <= q <= max_q:
            return q

    # если нет в таблице == максимальное доступное простое
    return max(primes)

def print_results_table(primes: List[int], test_results: List[bool], k: int, method_name: str) -> None:
    print(f"\nРезультаты для метода {method_name}:")
    print("+----+------+---------------------+")
    print("| N  | P    | Результат проверки |")
    print("+----+------+---------------------+")

    for i, (prime, result) in enumerate(zip(primes, test_results), 1):
        print(f"| {i:2} | {prime:4} | {'+':19} |" if result else f"| {i:2} | {prime:4} | {'-':19} |")

    print("+----+------+---------------------+")
    print(f"| K  | {k:4} |                     |")
    print("+----+------+---------------------+")

def main():
    # таблица простых чисел < 500
    primes = sieve_of_eratosthenes(500)

    bit_length = 10
    t = 5
    total_numbers = 10

    # генерация миллера
    miller_numbers = []
    miller_test_results = []
    k_miller = 0

    for _ in range(total_numbers):
        prime = generate_prime_miller(bit_length, primes, t)
        result = miller_test(prime, [], t)
        miller_numbers.append(prime)
        miller_test_results.append(result)

        if not result:
            second_test = miller_test(prime, [], t * 2)
            if second_test:
                k_miller += 1

    # генерация поклингтона
    pocklington_numbers = []
    pocklington_test_results = []
    k_pocklington = 0

    for _ in range(total_numbers):
        prime = generate_prime_pocklington(bit_length, primes, t)
        result = miller_test(prime, [], t)
        pocklington_numbers.append(prime)
        pocklington_test_results.append(result)

        if not result:
            second_test = miller_test(prime, [], t * 2)
            if second_test:
                k_pocklington += 1

    # генерация гост
    gost_numbers = []
    gost_test_results = []
    k_gost = 0
    generated = 0

    while generated < total_numbers:
        q = find_suitable_q(bit_length, primes)
        prime = generate_prime_gost(bit_length, q)
        if prime != -1:
            result = miller_test(prime, [], t)
            gost_numbers.append(prime)
            gost_test_results.append(result)

            if not result:
                second_test = miller_test(prime, [], t * 2)
                if second_test:
                    k_gost += 1
            generated += 1

    # вывод
    print_results_table(miller_numbers, miller_test_results, k_miller, "Миллера")
    print_results_table(pocklington_numbers, pocklington_test_results, k_pocklington, "Поклингтона")
    print_results_table(gost_numbers, gost_test_results, k_gost, "ГОСТ Р 34.10-94")

if __name__ == "__main__":
    main()