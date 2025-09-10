import math

def main():
    start = -4.0
    end = 6.0
    dx = 0.1

    print("   x\t|\ty")
    print("--------------------")

    x = start
    while x <= end:
        # прямая от x = -4 до x = 0
        if x < 0:
            y = -0.5 * x
        # четвертая четверть окружности
        elif x <= 2:
            under_sqrt = 4 - (x - 0) ** 2
            y = 2 - math.sqrt(under_sqrt) if under_sqrt >= 0 else float('nan')
        # первая четверть окружности
        elif x <= 4:
            under_sqrt = 4 - (x - 2) ** 2
            y = math.sqrt(under_sqrt) if under_sqrt >= 0 else float('nan')
        # прямая от x = 4 до x = 6
        else:
            y = -0.5 * x + 2

        print(f"{x:5.2f}\t|\t{y:.2f}")
        x += dx

if __name__ == "__main__":
    main()