def count_collisions(m1, m2, v1=0.0, v2=-1.0, ext=0.0):
    collisions = 0
    while True:
        if v1 < 0:
            v1 = -v1
            collisions += 1
            v1 *= (1 - ext)
            v2 *= (1 - ext)
        else:
            new_v1 = ((m1 - m2)/(m1 + m2))*v1 + (2*m2/(m1 + m2))*v2
            new_v2 = (2*m1/(m1 + m2))*v1 + ((m2 - m1)/(m1 + m2))*v2
            v1, v2 = new_v1, new_v2
            collisions += 1
            v1 *= (1 - ext)
            v2 *= (1 - ext)
        if 0 <= v1 < v2 and v2 >= 0:
            break
    return collisions

n = int(input("Введите n (степень 10 для отношения масс): "))
m1 = 1.0
m2 = 10.0**n * m1
result = count_collisions(m1, m2, v1=0.0, v2=-1.0)
print(f"При m2/m1=10^{n} число столкновений: {result}")