import sys

def main():
    n, m = map(int, sys.stdin.readline().split())
    nums = list(map(int, sys.stdin.readline().split()))
    
    prefix_sum = [0] * (n + 1)
    for i in range(n):
        prefix_sum[i + 1] = prefix_sum[i] + nums[i]
    
    dp = [-float('inf')] * (n + 1)
    dp[n] = 0
    
    for i in range(n - 1, -1, -1):
        for k in range(1, m + 1):
            if i + k > n:
                break
            current_sum = prefix_sum[i + k] - prefix_sum[i]
            dp[i] = max(dp[i], current_sum - dp[i + k])
    
    print(1 if dp[0] > 0 else 0)

if __name__ == "__main__":
    main()