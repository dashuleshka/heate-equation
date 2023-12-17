import numpy as np
import matplotlib.pyplot as plt
# Заданные параметры
R = 2.0 # Радиус цилиндра
a_square = R **2 / np.pi **2 # Квадрат коэффициента теплопроводности a ^ 2
N_r = 10 # Количество шагов по радиусу
t_max = 1.0 # Максимальное время
# Шаги по радиусу и времени
h = R / N_r

# Вычисление шага по времени с учетом условия Куранта
courant_factor = 0.5 # Задаем коэффициент для условия Куранта
tau = courant_factor * h **2 / a_square # Вычисляем шаг по времени
N_t = int(t_max / tau) # Пересчитываем количество шагов по времени

# Начальное распределение температуры
def u0(r) :
    if r == 0 :
        return 0
    return np.sin(np.pi * r / R) / np.sqrt(r)

 # Инициализация массива для решения
nu = np.zeros((N_r + 1, N_t + 1))
r = np.linspace(0, R, N_r + 1)

 # Установка начального распределения температуры
for i in range(N_r + 1) :
    nu[i, 0] = np.sqrt(r[i]) * u0(r[i])

 # Определение функции погрешности f(r, t)
def f(r, t) :
    return (R **2 / (4 * np.pi **2 * r **2)) * np.sin(np.pi * r / R) * np.exp(-t)
 
 # Применение разностной схемы
for j in range(N_t) : # По времени
    for i in range(1, N_r) : # По радиусу(внутренние точки)
        nu[i, j + 1] = (nu[i, j] + (a_square * (nu[i + 1, j] - 2 * nu[i, j] + nu[i - 1, j]) / h **2 + a_square * nu[i, j] / (4 * r[i] **2)
                    - f(r[i], j * tau)) * tau)

 # Визуализация численного решения
R_grid, T_grid = np.meshgrid(r, np.linspace(0, t_max, N_t + 1))

def analytical_solution(r, t) :
    return np.sin(np.pi * r / R) * np.exp(-t)

difference = nu.T - analytical_solution(R_grid, T_grid)
fig = plt.figure()
ax = fig.add_subplot(111, projection = "3d")
p = ax.plot_surface(R_grid, T_grid, difference, cmap = "viridis")
ax.set_xlabel("r")
ax.set_ylabel("t")
ax.set_zlabel("nu(r, t)")
ax.set_title("График погрешностей")
plt.show()

 # Аналитическое решение
def analytical_solution(r, t) :
    return np.sin(np.pi * r / R) * np.exp(-t)

 # Расчет разности между численным и аналитическим решениями
difference = nu.T - analytical_solution(R_grid, T_grid)
norm_difference = np.max(np.abs(difference))
print(f"Норма разности между численным и аналитическим решениями: {norm_difference}")
