#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>

using namespace std;

int main() {
    int n, m;
    cin >> n >> m;

    vector<int> nums(n);
    for (int i = 0; i < n; ++i) {
        cin >> nums[i];
    }

    vector<long long> prefix_sum(n + 1, 0);
    for (int i = 0; i < n; ++i) {
        prefix_sum[i + 1] = prefix_sum[i] + nums[i];
    }

    vector<long long> dp(n + 1, LLONG_MIN);
    dp[n] = 0;

    for (int i = n - 1; i >= 0; --i) {
        for (int k = 1; k <= m; ++k) {
            if (i + k > n) break;
            long long current_sum = prefix_sum[i + k] - prefix_sum[i];
            dp[i] = max(dp[i], current_sum - dp[i + k]);
        }
    }

    cout << (dp[0] > 0 ? 1 : 0) << endl;

    return 0;
}