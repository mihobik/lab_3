#include <iostream>
#include <vector>
#include <string>

using namespace std;

// нод
long long gcd(long long a, long long b) {
    return b == 0 ? a : gcd(b, a % b);
}

// коэффициенты полиномов эйлера (для a)
vector<vector<long long>> euler = {
    {}, // 0
    {1},
    {1, 1},
    {1, 4, 1},
    {1, 11, 11, 1},
    {1, 26, 66, 26, 1},
    {1, 57, 302, 302, 57, 1},
    {1, 120, 1191, 2416, 1191, 120, 1},
    {1, 247, 4293, 15619, 15619, 4293, 247, 1},
    {1, 502, 14608, 88234, 156190, 88234, 14608, 502, 1},
    {1, 1013, 46834, 472063, 1561900, 1561900, 472063, 46834, 1013, 1} // 10
};

// вычисление суммы как несократимой дроби
string rational_sum(int a, int b) {
    if (a < 1 || a > 10) return "";

    const vector<long long>& coeffs = euler[a];
    int deg = coeffs.size() - 1;

    // вычисление числителя
    long long num = 0;
    long long pow_b = 1;
    for (int i = deg; i >= 0; i--) {
        num = num * b + coeffs[i]; // схема горнера
    }
    num *= b;

    // знаменатель (b-1)^(a+1)
    long long den = 1;
    for (int i = 0; i < a + 1; i++) den *= (b - 1);

    long long g = gcd(num, den);
    num /= g;
    den /= g;

    return to_string(num) + "/" + to_string(den);
}

int main() {
    int a, b;
    cin >> a >> b;

    if (b == 1) {
        cout << "infinity" << endl;
        return 0;
    }

    string res = rational_sum(a, b);
    if (res.empty()) cout << "irrational" << endl;
    else cout << res << endl;

    return 0;
}
