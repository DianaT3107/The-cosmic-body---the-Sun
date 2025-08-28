# подключение необходимых библитек
from scipy import integrate
import numpy as np
from matplotlib import pyplot as plt

#константы и "дано"
G = 6.67E-11 # гравитационная постоянная
M = 1.989E30 # масса Солнца (кг)
r = 1.49e13 # расстояние от Солнца до космического тела (м) (100а.е.)
v = 3000 # начальная скорость космического тела (м/с)
n = 100 # число точек интегрирования
t_start = 0 #начальное значение массива точек интегрирования
t_stop = 3.19e10 #конечное значение массива точек интегрирования
m = 2.2e14 # масса космического объекта (в кг)
B = G * M * m #константа введённая для упрощения записи
a = - G * M #константа введённая для упрощения записи
#нулевые массивы для записи данных
r_an = np.zeros(n)
x_an = np.zeros(n)
y_an = np.zeros(n)

##########   численное решение ##############

# полагается, что f зависит от (ode_vec, t), причём ode_vec - это [x, y, vx, vy] 
def f(ode_vec, t):
    x = ode_vec[0]
    y = ode_vec[1]
    vx = ode_vec[2]
    vy = ode_vec[3]
    
    #переменные введённые для упрощения записи
    R = np.sqrt(x**2 + y**2) 
    
    # система дифф. уравнений
    f0 = vx
    f1 = vy
    f2 = a * x / (R**3)
    f3 = a * y / (R**3)
    return [f0, f1, f2 , f3]

#время в секундах
t = np.linspace(t_start, t_stop, n) # массив точек интегрирования 

# начальные условия
x0 = r
y0 = 0
vx0 = 0
vy0 = v
ic = [x0, y0, vx0, vy0]

#интегрирование с помощью библеотеки SciPy
sol = integrate.odeint(f, ic, t)

#переобозначение для удобства
x = sol[:, 0] 
y = sol[:, 1] 

################### аналитическое решение ##########
L = np.outer(r, m * v) #момент импульса 
E = L ** 2 / (2 * m * r ** 2) - B / r  + m / (2 * v ** 2 )
# создание функции для полярного радиса из уравнения траектории движения
def r_exact(phi):
   #полная энергия
    
    #замена переменных
    p = L ** 2 / (G * M * m **2) #параметр
    e = np.sqrt(1 + (2 * E * L ** 2) / (G ** 2 * M ** 2 * m ** 3)) #эксцентриситет
    
    return p / (1 + e * np.cos(phi)) # уравнение траектории движения в виде конического сечения 

fi = np.linspace(0.0, 2.0 * np.pi, n) #массив точек

#цикл записи значений полярного радиуса
for i in range(n):
    r_an[i] = r_exact(fi[i])

#переход из полярных координат в декартовые 
for i in range(n) :
    x_an[i] = r_an[i] * np.cos(fi[i])
    y_an[i] = r_an[i] * np.sin(fi[i])
    
#построение орбиты в декартовых координатах
fig = plt.figure(figsize=(20, 20)) 

ax1 = fig.add_subplot(111)
plt.axis([-2e13, 2e13, -2e13, 2e13])
ax1.plot(x, y, 'k--' , label='Численное решение', linewidth = 2)   
ax1.plot(x_an, y_an, color = 'darkorange', label='аналитическое решение', linewidth = 1)
ax1.set_xlabel('x (m)')
ax1.set_ylabel('y (m)')
ax1.scatter(0, 0, color='orange', label='Sun',  marker='*', s = 200)
ax1.scatter(x0, y0, color='darkred', label='space body', marker='o', s = 40)
ax1.grid(True)
ax1.legend(loc='lower left')


fig.savefig("Orbita.png", orientation = 'landscape', dpi=500) 

