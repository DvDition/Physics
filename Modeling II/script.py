import numpy as np
import matplotlib.pyplot as plt

e = 1.602176634e-19
me = 9.10938356e-31

r = 8.5e-2
R = 18e-2
L = 26e-2
Vx = 1e6

gap = R - r
y0 = gap / 2.0
vy0 = 0.0

lnRR = np.log(R / r)
t_exit = L / Vx

def ay(y, U):
    return -e * U / (me * (r + y) * lnRR)

def simulate(U, dt=1e-11, save=False, max_steps=10_000_000):
    t = 0.0
    y = y0
    vy = vy0

    if save:
        X = [0.0]
        T = [0.0]
        Y = [y]
        VY = [vy]
        AY = [ay(y, U)]

    steps = 0
    while t < t_exit and steps < max_steps:
        y_old, vy_old, t_old = y, vy, t

        def f1(y, vy): return vy
        def f2(y, vy): return ay(y, U)

        k1y = f1(y, vy)
        k1v = f2(y, vy)

        k2y = f1(y + 0.5*dt*k1y, vy + 0.5*dt*k1v)
        k2v = f2(y + 0.5*dt*k1y, vy + 0.5*dt*k1v)

        k3y = f1(y + 0.5*dt*k2y, vy + 0.5*dt*k2v)
        k3v = f2(y + 0.5*dt*k2y, vy + 0.5*dt*k2v)

        k4y = f1(y + dt*k3y, vy + dt*k3v)
        k4v = f2(y + dt*k3y, vy + dt*k3v)

        y = y + dt*(k1y + 2*k2y + 2*k3y + k4y)/6.0
        vy = vy + dt*(k1v + 2*k2v + 2*k3v + k4v)/6.0
        t += dt
        steps += 1

        if y <= 0.0 or y >= gap:
            if y <= 0.0:
                boundary = 0.0
            else:
                boundary = gap
            denom = (y - y_old)
            if abs(denom) < 1e-16:
                alpha = 0.0
            else:
                alpha = (boundary - y_old) / denom
                alpha = max(0.0, min(1.0, alpha))

            t_coll = t_old + alpha * (t - t_old)
            vy_coll = vy_old + alpha * (vy - vy_old)

            if save:
                X.append(Vx * t_coll)
                T.append(t_coll)
                Y.append(boundary)
                VY.append(vy_coll)
                AY.append(ay(boundary, U))

                return True, (np.array(X), np.array(T), np.array(Y), np.array(VY), np.array(AY))
            else:
                return True, None

        if save:
            X.append(Vx * t)
            T.append(t)
            Y.append(y)
            VY.append(vy)
            AY.append(ay(y, U))
    if save:
        return False, (np.array(X), np.array(T), np.array(Y), np.array(VY), np.array(AY))
    else:
        return False, None

def find_min_U(U_start=1.0, dt=1e-11):
    U_low = 0.0
    U_high = U_start
    for _ in range(60):
        collided, _ = simulate(U_high, dt=dt, save=False)
        if collided:
            break
        U_high *= 2.0
    else:
        raise RuntimeError("Не удалось найти U_high при котором происходит столкновение")

    for _ in range(40):
        U_mid = 0.5 * (U_low + U_high)
        collided, _ = simulate(U_mid, dt=dt, save=False)
        if collided:
            U_high = U_mid
        else:
            U_low = U_mid

    return U_high

dt = 1e-11
U_min = find_min_U(U_start=5.0, dt=dt)
collided, data = simulate(U_min, dt=dt, save=True)
X, T, Y, VY, AY = data

t_sim = T[-1]
vy_final = VY[-1]
V_final = np.sqrt(Vx**2 + vy_final**2)

print(f"U_min = {U_min:.6f} В")
print(f"t_exit (теоретич.) = {t_exit:.6e} с")
print(f"t_sim (реально, по симуляции) = {t_sim:.6e} с")
print(f"V_кон = {V_final:.6e} м/с")

plt.figure(figsize=(10, 8))
plt.subplot(2,2,1)
plt.plot(X, Y)
plt.title("y(x)")
plt.xlabel("x, м")
plt.ylabel("y, м")
plt.grid()

plt.subplot(2,2,2)
plt.plot(T, VY)
plt.title("v_y(t)")
plt.xlabel("t, с")
plt.ylabel("v_y, м/с")
plt.grid()

plt.subplot(2,2,3)
plt.plot(T, AY)
plt.title("a_y(t)")
plt.xlabel("t, с")
plt.ylabel("a_y, м/с²")
plt.grid()

plt.subplot(2,2,4)
plt.plot(T, Y)
plt.title("y(t)")
plt.xlabel("t, с")
plt.ylabel("y, м")
plt.grid()

plt.tight_layout()
plt.show()