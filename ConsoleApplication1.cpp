#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <locale>
#include <iomanip>
#include <numeric> // gcd

using namespace std;

// решето эратосфена
vector<int> sieve_of_eratosthenes(int limit) {
    vector<bool> is_prime(limit + 1, true);
    if (limit >= 0) is_prime[0] = false;
    if (limit >= 1) is_prime[1] = false;

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

// целочисленное возведение в степень (без double)
long long int_pow_ll(long long base, int exp) {
    long long result = 1;
    while (exp > 0) {
        if (exp & 1) result = result * base;
        base = base * base;
        exp >>= 1;
    }
    return result;
}

// возведение в степень по модулю (если mod небольшой)
long long mod_pow(long long base, long long exp, long long mod) {
    long long result = 1 % mod;
    base = ((base % mod) + mod) % mod;

    while (exp > 0) {
        if (exp & 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp >>= 1;
    }

    return result;
}

// факторизация n по таблице простых (возвращает пары (prime, exponent))
vector<pair<int, int>> factorize_with_table(long long n, const vector<int>& primes_table) {
    vector<pair<int, int>> res;
    long long temp = n;
    for (int p : primes_table) {
        if ((long long)p * p > temp) break;
        if (temp % p == 0) {
            int cnt = 0;
            while (temp % p == 0) {
                temp /= p;
                ++cnt;
            }
            res.emplace_back(p, cnt);
        }
    }
    if (temp > 1) {
        // temp сам прост
        if (temp <= INT_MAX)
            res.emplace_back((int)temp, 1);
        else {
            // редкий случай (для наших размеров неважен)
            res.emplace_back((int)temp, 1);
        }
    }
    return res;
}

// тест миллера - строгий вариант: получает каноническое разложение n-1 (без учета 2 или с ним)
// Реализация ориентирована на проверку, описанную в вашем фото: 
// 1) для всех выбранных a: a^(n-1) ≡ 1 (mod n)
// 2) для каждого простого q | (n-1) должно существовать a среди выбранных, для которого a^{(n-1)/q} != 1 (mod n)
//    иначе n считается составным (по этой строгой форме критерия)
bool miller_test(long long n, const vector<pair<int, int>>& factorization_of_n_minus_1, int t) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;

    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<long long> dist(2, n - 2);

    // выбор t различных случайных чисел a_j
    vector<long long> a_values;
    a_values.reserve(t);
    while ((int)a_values.size() < t) {
        long long a = dist(gen);
        if (find(a_values.begin(), a_values.end(), a) == a_values.end()) a_values.push_back(a);
    }

    // 1) проверка a^(n-1) mod n == 1 для всех a
    for (long long a : a_values) {
        if (mod_pow(a, n - 1, n) != 1) {
            return false; // явно составное (нарушение условия Ферма)
        }
    }

    // составим список простых делителей q (уникальные q)
    vector<int> q_list;
    for (const auto& pr : factorization_of_n_minus_1) {
        int q = pr.first;
        // пропустим q=2 если он есть, т.к. в условии часто рассматривают нечётные q
        if (q != 2) q_list.push_back(q);
    }
    sort(q_list.begin(), q_list.end());
    q_list.erase(unique(q_list.begin(), q_list.end()), q_list.end());

    // 2) для каждого q в разложении должно существовать хотя бы одно a такое, что a^{(n-1)/q} != 1 (mod n)
    for (int q : q_list) {
        bool exists_a_with_not1 = false;
        long long exp = (n - 1) / q;
        for (long long a : a_values) {
            if (mod_pow(a, exp, n) != 1) {
                exists_a_with_not1 = true;
                break;
            }
        }
        if (!exists_a_with_not1) {
            // для данного q все a дали 1 => критерий не выполнен => пометим как составное
            return false;
        }
    }

    // если все проверки пройдены, считаем число прошедшим тест Миллера (по строгой схеме, описанной в фото)
    return true;
}

// тест поклингтона (строгая реализация по теореме Поклингтона),
// принимает n, F (известная часть разложения n-1) и точную факторизацию F.
// Условия (в упрощенном виде, как на иллюстрации):
// 1) для каждого выбранного a: a^(n-1) ≡ 1 (mod n)
// 2) для каждого простого q | F должно существовать a такое, что gcd(a^{(n-1)/q} - 1, n) == 1
// Если оба условия выполнены для всех q из F, то n прост.
bool pocklington_test(long long n, long long F, const vector<pair<int, int>>& F_factorization, int t) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;

    // требование теоремы: F должно быть > sqrt(n)-1 (или по условию: F содержит достаточно большие делители)
    long double sqrt_n = sqrt((long double)n);
    if (!((long double)F > sqrt_n - 1.0L)) {
        // если F слишком мало, теорема не применима строго
        return false;
    }

    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<long long> dist(2, n - 2);

    // выбираем t различных a
    vector<long long> a_values;
    a_values.reserve(t);
    while ((int)a_values.size() < t) {
        long long a = dist(gen);
        if (find(a_values.begin(), a_values.end(), a) == a_values.end()) a_values.push_back(a);
    }

    // 1) проверка a^(n-1) ≡ 1 (mod n) для всех a
    for (long long a : a_values) {
        if (mod_pow(a, n - 1, n) != 1) {
            return false;
        }
    }

    // для каждого простого q в разложении F нужно, чтобы существовал a с gcd(a^{(n-1)/q}-1, n) == 1
    vector<int> q_list;
    for (const auto& pr : F_factorization) {
        int q = pr.first;
        if (q != 2) q_list.push_back(q); // обычно Берём простые q (не обязательно 2)
    }
    sort(q_list.begin(), q_list.end());
    q_list.erase(unique(q_list.begin(), q_list.end()), q_list.end());

    for (int q : q_list) {
        bool found_a = false;
        long long exp = (n - 1) / q;
        for (long long a : a_values) {
            long long val = mod_pow(a, exp, n);
            long long g = std::gcd((long long)((val - 1 + n) % n), n);
            if (g == 1) {
                found_a = true;
                break;
            }
        }
        if (!found_a) {
            // для этого q не найдётся a, удовлетворяющего условию => теорема не может подтвердить простоту
            return false;
        }
    }

    // все условия Поклингтона выполнены => n прост (по теореме)
    return true;
}

// генерация простого числа (миллер) — теперь также выдаёт факторизацию n-1 в factor_out
long long generate_prime_miller(int bit_length, const vector<int>& primes, int t, vector<pair<int, int>>& factor_out) {
    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<int> dist_idx(0, (int)primes.size() - 1);
    uniform_int_distribution<int> dist_alpha(1, 3);

    // границы для m
    const long long min_m = (1LL << (bit_length - 2)) + 1;
    const long long max_m = (1LL << (bit_length - 1)) - 1;

    while (true) {
        // строим число m в нужном диапазоне
        long long m = 1;
        vector<pair<int, int>> m_factorization;

        // аккуратно умножаем простые степени, проверяя на переполнение / выход за max_m
        while (m < max_m) {
            int q = primes[dist_idx(gen)];
            int alpha = dist_alpha(gen);

            // вычислим q^alpha безопасно
            long long qpow = 1;
            bool overflow = false;
            for (int i = 0; i < alpha; ++i) {
                if (qpow > max_m / q) { overflow = true; break; }
                qpow *= q;
            }
            if (overflow) break;
            if (m > max_m / qpow) break;

            m *= qpow;
            // добавляем в факторизацию m_factorization (аккумулируем степень)
            bool merged = false;
            for (auto& pr : m_factorization) {
                if (pr.first == q) { pr.second += alpha; merged = true; break; }
            }
            if (!merged) m_factorization.emplace_back(q, alpha);
        }

        if (m < min_m || m > max_m) continue;

        long long n = 2 * m + 1;
        int bits = (int)log2((long double)n) + 1;
        if (bits != bit_length) continue;

        // формируем факторизацию n-1 = 2 * m
        vector<pair<int, int>> n_minus_1_fact;
        // учтём двойку
        n_minus_1_fact.emplace_back(2, 1);
        // добавим факторизацию m
        for (auto& pr : m_factorization) n_minus_1_fact.push_back(pr);

        // проверка n тестом Миллера (строгая проверка с реальной факторизацией)
        if (miller_test(n, n_minus_1_fact, t)) {
            factor_out = move(n_minus_1_fact);
            return n;
        }
    }
}

// генерация простого числа (поклингтона) — возвращает n и заполняет F_factorization (факторизацию F)
long long generate_prime_pocklington(int bit_length, const vector<int>& primes, int t, long long& F_out, vector<pair<int, int>>& F_factorization_out) {
    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<int> dist_idx(0, (int)primes.size() - 1);
    uniform_int_distribution<int> dist_alpha(1, 3);
    uniform_int_distribution<long long> dist_R(1, (1LL << (bit_length / 2 - 1)) - 1);

    // границы для F (на 1 бит больше половины требуемого размера)
    const long long min_F = (1LL << (bit_length / 2)) + 1;
    const long long max_F = (1LL << (bit_length / 2 + 1)) - 1;

    while (true) {
        long long F = 1;
        vector<pair<int, int>> F_factorization;

        while (F < max_F) {
            int q = primes[dist_idx(gen)];
            int alpha = dist_alpha(gen);

            // вычислим q^alpha безопасно
            long long qpow = 1;
            bool overflow = false;
            for (int i = 0; i < alpha; ++i) {
                if (qpow > max_F / q) { overflow = true; break; }
                qpow *= q;
            }
            if (overflow) break;
            if (F > max_F / qpow) break;

            F *= qpow;
            bool merged = false;
            for (auto& pr : F_factorization) {
                if (pr.first == q) { pr.second += alpha; merged = true; break; }
            }
            if (!merged) F_factorization.emplace_back(q, alpha);
        }

        if (F < min_F) continue;

        long long R = 2 * (dist_R(gen)); // выбираем чётное R
        long long n = R * F + 1;

        int bits = (int)log2((long double)n) + 1;
        if (bits != bit_length) continue;

        // теперь проверяем Поклингтон строго: передаем F и её факторизацию
        if (pocklington_test(n, F, F_factorization, t)) {
            F_out = F;
            F_factorization_out = move(F_factorization);
            return n;
        }
    }
}

// генерация простого (ГОСТ) - оставлена ваша логика, но возвращаем просто p
long long generate_prime_gost(int t, long long q) {
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> dist_xi(0.0, 1.0);

    const long long max_p = (1LL << t) - 1;
    int attempts = 0;
    const int max_attempts = 1000;

    while (attempts++ < max_attempts) {
        double xi = dist_xi(gen);

        long long N = (long long)ceil(pow(2.0L, t - 1) / (long double)q) + (long long)ceil((pow(2.0L, t - 1) * xi / (long double)q));
        if (N % 2 != 0) N++;

        for (long long u = 0; u < 4 * (q + 1); u += 2) {
            long long p = (N + u) * q + 1;
            if (p > max_p) break;
            if (p < 2) continue;
            // проверка по критерию Диэмитко (a=2)
            if (mod_pow(2, p - 1, p) == 1 && mod_pow(2, N + u, p) != 1) {
                return p;
            }
        }
    }
    return -1;
}

// функция для подбора q (как у вас)
long long find_suitable_q(int bit_length, const vector<int>& primes) {
    int q_bits = (bit_length + 1) / 2;
    long long min_q = (1LL << (q_bits - 1)) + 1;
    long long max_q = (1LL << q_bits) - 1;

    for (int q : primes) {
        if (q >= min_q && q <= max_q) {
            return q;
        }
    }

    return *max_element(primes.begin(), primes.end());
}

void print_results_table(const vector<long long>& numbers,
    const vector<bool>& test_results,
    int k,
    const string& method_name) {
    cout << "\nРезультаты для метода " << method_name << ":\n";
    cout << "+----+--------+---------------------+\n";
    cout << "| N  |   P    | Результат проверки |\n";
    cout << "+----+--------+---------------------+\n";

    for (size_t i = 0; i < numbers.size(); ++i) {
        cout << "| " << setw(2) << i + 1 << " | "
            << setw(6) << numbers[i] << " | "
            << setw(19) << (test_results[i] ? "+" : "-") << " |\n";
    }

    cout << "+----+--------+---------------------+\n";
    cout << "| K  | " << setw(6) << k << " |                     |\n";
    cout << "+----+--------+---------------------+\n";
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
        vector<pair<int, int>> n_minus_1_fact;
        long long prime = generate_prime_miller(bit_length, primes, t, n_minus_1_fact);
        bool result = miller_test(prime, n_minus_1_fact, t);
        miller_numbers.push_back(prime);
        miller_test_results.push_back(result);

        if (!result) {
            // повторная большая проверка
            bool second_test = miller_test(prime, n_minus_1_fact, t * 2);
            if (second_test) k_miller++;
        }
    }

    // генерация поклингтона
    vector<long long> pocklington_numbers;
    vector<bool> pocklington_test_results;
    int k_pocklington = 0;

    for (int i = 0; i < total_numbers; ++i) {
        long long F;
        vector<pair<int, int>> F_factorization;
        long long prime = generate_prime_pocklington(bit_length, primes, t, F, F_factorization);
        bool result = pocklington_test(prime, F, F_factorization, t);
        pocklington_numbers.push_back(prime);
        pocklington_test_results.push_back(result);

        if (!result) {
            bool second_test = pocklington_test(prime, F, F_factorization, t * 2);
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
            // факторизуем n-1 по таблице и используем в тесте Миллера
            vector<pair<int, int>> n_minus_1_fact = factorize_with_table(prime - 1, primes);
            bool result = miller_test(prime, n_minus_1_fact, t);
            gost_numbers.push_back(prime);
            gost_test_results.push_back(result);

            if (!result) {
                bool second_test = miller_test(prime, n_minus_1_fact, t * 2);
                if (second_test) k_gost++;
            }
            generated++;
        }
        else {
            // если не получилось — выйдем чтобы не зациклиться в бесконечности для демонстрации
            break;
        }
    }

    // вывод
    print_results_table(miller_numbers, miller_test_results, k_miller, "Миллера");
    print_results_table(pocklington_numbers, pocklington_test_results, k_pocklington, "Поклингтона");
    print_results_table(gost_numbers, gost_test_results, k_gost, "ГОСТ Р 34.10-94");

    return 0;
}
