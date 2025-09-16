#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <locale>
#include <iomanip>
#include <locale>

using namespace std;

// решето эратосфена
vector<int> sieve_of_eratosthenes(int limit) {
    vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;

    for (int p = 2; p * p <= limit; ++p) {
        if (is_prime[p]) {
            for (int i = p * p; i <= limit; i += p) {
                is_prime[i] = false;
            }
        }
    }

    vector<int> primes;
    for (int i = 2; i <= limit; ++i) {
        if (is_prime[i]) {
            primes.push_back(i);
        }
    }

    return primes;
}

// возведение в степень по модулю
long long mod_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    base = base % mod;

    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * base) % mod;
        }
        exp = exp >> 1;
        base = (base * base) % mod;
    }

    return result;
}

// тест миллера
bool miller_test(long long n, const vector<pair<int, int>>& factorization, int t) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<long long> dist(2, n - 2);

    // выбор t различных случайных чисел a_j
    vector<long long> a_values;
    for (int i = 0; i < t; ++i) {
        long long a;
        do {
            a = dist(gen);
        } while (find(a_values.begin(), a_values.end(), a) != a_values.end());
        a_values.push_back(a);
    }

    // проверка a_j^(n-1) mod n == 1
    for (long long a : a_values) {
        if (mod_pow(a, n - 1, n) != 1) {
            return false;
        }
    }

    // для каждого простого делителя q_i
    for (const auto& factor : factorization) {
        int q = factor.first;
        bool all_one = true;
        for (long long a : a_values) {
            long long exponent = (n - 1) / q;
            if (mod_pow(a, exponent, n) == 1) {
                all_one = false;
                break;
            }
        }
        if (all_one) return false;
    }

    return true;
}

// генерация простого числа (миллер)
long long generate_prime_miller(int bit_length, const vector<int>& primes, int t) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist_idx(0, primes.size() - 1);
    uniform_int_distribution<int> dist_alpha(1, 3);

    // границы для m
    const long long min_m = (1LL << (bit_length - 2)) + 1;
    const long long max_m = (1LL << (bit_length - 1)) - 1;

    while (true) {
        // строим число m в нужном диапазоне
        long long m;
        vector<pair<int, int>> factorization;

        do {
            m = 1;
            factorization.clear();

            while (m < max_m) {
                int idx = dist_idx(gen);
                int q = primes[idx];
                int alpha = dist_alpha(gen);

                // проверка на переполнение
                if (m > max_m / pow(q, alpha)) break;

                m *= pow(q, alpha);
                factorization.emplace_back(q, alpha);
            }
        } while (m < min_m || m > max_m);

        // вычисляем n = 2m + 1
        long long n = 2 * m + 1;

        // доп проверка битовой длины
        int bits = (int)log2(n) + 1;
        if (bits != bit_length) continue;

        // проверка n тестом Миллера
        if (miller_test(n, factorization, t)) {
            return n;
        }
    }
}

// тест поклингтона 
bool pocklington_test(long long n, long long F, const vector<pair<int, int>>& F_factorization, int t) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;

    long long R = (n - 1) / F;
    if (F <= sqrt(n) - 1) return false;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<long long> dist(2, n - 2);

    // выбор t различных случайных чисел a_j
    vector<long long> a_values;
    for (int i = 0; i < t; ++i) {
        long long a;
        do {
            a = dist(gen);
        } while (find(a_values.begin(), a_values.end(), a) != a_values.end());
        a_values.push_back(a);
    }

    // проверка a_j^(n-1) mod n == 1
    for (long long a : a_values) {
        if (mod_pow(a, n - 1, n) != 1) {
            return false;
        }
    }

    // для каждого a_j проверить условия
    for (long long a : a_values) {
        bool condition_met = false;
        for (const auto& factor : F_factorization) {
            int q = factor.first;
            long long exponent = (n - 1) / q;
            if (mod_pow(a, exponent, n) != 1) {
                condition_met = true;
                break;
            }
        }

        if (!condition_met) {
            return false;
        }
    }

    return true;
}

// генерация простого числа (поклингтона)
long long generate_prime_pocklington(int bit_length, const vector<int>& primes, int t) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist_idx(0, primes.size() - 1);
    uniform_int_distribution<int> dist_alpha(1, 3);
    uniform_int_distribution<long long> dist_R(1, (1LL << (bit_length / 2 - 1)) - 1);

    // границы для F (на 1 бит больше половины требуемого размера)
    const long long min_F = (1LL << (bit_length / 2)) + 1;
    const long long max_F = (1LL << (bit_length / 2 + 1)) - 1;

    while (true) {
        // построение числа F
        long long F = 1;
        vector<pair<int, int>> F_factorization;

        while (F < max_F) {
            int idx = dist_idx(gen);
            int q = primes[idx];
            int alpha = dist_alpha(gen);

            // проверка на переполнение
            if (F > max_F / pow(q, alpha)) break;

            F *= pow(q, alpha);
            F_factorization.emplace_back(q, alpha);
        }

        if (F < min_F) continue;

        // выбор случайного четного R
        long long R = 2 * dist_R(gen);
        long long n = R * F + 1;

        // проверка размера
        int bits = (int)log2(n) + 1;
        if (bits != bit_length) continue;

        // проверка n тестом поклингтона
        if (pocklington_test(n, F, F_factorization, t)) {
            return n;
        }
    }
}

// генерация простого числа (гост)
long long generate_prime_gost(int t, long long q) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist_xi(0.0, 1.0);

    const long long max_p = (1LL << t) - 1;
    int attempts = 0;
    const int max_attempts = 1000;

    while (attempts++ < max_attempts) {
        double xi = dist_xi(gen);

        // вычисление N
        long long N = ceil(pow(2, t - 1) / q) + ceil((pow(2, t - 1) * xi / q));
        if (N % 2 != 0) N++;

        // инициализация u
        for (long long u = 0; u < 4 * (q + 1); u += 2) {
            // Шаг 3: Вычислить p
            long long p = (N + u) * q + 1;

            // проверка размера p
            if (p > max_p) break;

            // проверка теоремы диемитко при a = 2
            if (mod_pow(2, p - 1, p) == 1 && mod_pow(2, N + u, p) != 1) {
                return p;
            }
        }
    }
    return -1; // если не удалось найти простое число
}

// функция для подбора q
long long find_suitable_q(int bit_length, const vector<int>& primes) {
    int q_bits = (bit_length + 1) / 2;
    long long min_q = (1LL << (q_bits - 1)) + 1;
    long long max_q = (1LL << q_bits) - 1;

    // ищем q в таблице
    for (int q : primes) {
        if (q >= min_q && q <= max_q) {
            return q;
        }
    }

    // если нет в таблице == максимальное доступное простое
    return *max_element(primes.begin(), primes.end());
}

void print_results_table(const vector<long long>& primes,
    const vector<bool>& test_results,
    int k,
    const string& method_name) {
    cout << "\nРезультаты для метода " << method_name << ":\n";
    cout << "+----+------+---------------------+\n";
    cout << "| N  | P    | Результат проверки |\n";
    cout << "+----+------+---------------------+\n";

    // данные
    for (size_t i = 0; i < primes.size(); ++i) {
        cout << "| " << setw(2) << i + 1 << " | "
            << setw(4) << primes[i] << " | "
            << setw(19) << (test_results[i] ? "+" : "-") << " |\n";
    }

    cout << "+----+------+---------------------+\n";
    cout << "| K  | " << setw(4) << k << " |                     |\n";
    cout << "+----+------+---------------------+\n";
}

int main() {
    setlocale(LC_ALL, "Russian");

    // таблица простых чисел < 500
    vector<int> primes = sieve_of_eratosthenes(500);

    const int bit_length = 10;
    const int t = 5;
    const int total_numbers = 10;

    // генерация миллера
    vector<long long> miller_numbers;
    vector<bool> miller_test_results;
    int k_miller = 0;

    for (int i = 0; i < total_numbers; ++i) {
        long long prime = generate_prime_miller(bit_length, primes, t);
        vector<pair<int, int>> empty_fact;
        bool result = miller_test(prime, empty_fact, t);
        miller_numbers.push_back(prime);
        miller_test_results.push_back(result);

        if (!result) {
            bool second_test = miller_test(prime, empty_fact, t * 2);
            if (second_test) k_miller++;
        }
    }

    // генерация поклингтона
    vector<long long> pocklington_numbers;
    vector<bool> pocklington_test_results;
    int k_pocklington = 0;

    for (int i = 0; i < total_numbers; ++i) {
        long long prime = generate_prime_pocklington(bit_length, primes, t);
        vector<pair<int, int>> empty_fact;
        bool result = miller_test(prime, empty_fact, t);
        pocklington_numbers.push_back(prime);
        pocklington_test_results.push_back(result);

        if (!result) {
            bool second_test = miller_test(prime, empty_fact, t * 2);
            if (second_test) k_pocklington++;
        }
    }

    // генерация гост
    vector<long long> gost_numbers;
    vector<bool> gost_test_results;
    int k_gost = 0;
    int generated = 0;

    while (generated < total_numbers) {
        long long q = find_suitable_q(bit_length, primes);
        long long prime = generate_prime_gost(bit_length, q);
        if (prime != -1) {
            vector<pair<int, int>> empty_fact;
            bool result = miller_test(prime, empty_fact, t);
            gost_numbers.push_back(prime);
            gost_test_results.push_back(result);

            if (!result) {
                bool second_test = miller_test(prime, empty_fact, t * 2);
                if (second_test) k_gost++;
            }
            generated++;
        }
    }

    // вывод
    print_results_table(miller_numbers, miller_test_results, k_miller, "Миллера");
    print_results_table(pocklington_numbers, pocklington_test_results, k_pocklington, "Поклингтона");
    print_results_table(gost_numbers, gost_test_results, k_gost, "ГОСТ Р 34.10-94");

    return 0;
}
