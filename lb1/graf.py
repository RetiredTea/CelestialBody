import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def load_data(method):
    df = pd.read_csv(f'trajectories{method}.csv')
    df['Time(days)'] = df['Time(h)'] / 24
    return df.pivot_table(index='Time(days)', columns='Body',
                        values=['X(m)', 'Y(m)', 'Z(m)'], aggfunc='first')

rk4 = load_data('')
rk5 = load_data('3')

# Расчет отклонений
bodies = [b for b in rk4.columns.get_level_values(1).unique() if b != 'Sun']
time = rk4.index.values

instant_error = np.zeros(len(time))
cumulative_error = np.zeros(len(time))

for body in bodies:
    dx = rk4[('X(m)', body)] - rk5[('X(m)', body)]
    dy = rk4[('Y(m)', body)] - rk5[('Y(m)', body)]
    dz = rk4[('Z(m)', body)] - rk5[('Z(m)', body)]
    body_error = np.sqrt(dx**2 + dy**2 + dz**2)
    instant_error += body_error
    cumulative_error += np.cumsum(body_error)

instant_error /= len(bodies)
cumulative_error /= len(bodies)

# Нормализация ошибок
instant_error_au = instant_error / 1.496e11
cumulative_error_au = cumulative_error / 1.496e11

# Визуализация
plt.style.use('dark_background')
fig = plt.figure(figsize=(20, 12), facecolor='#0f0f0f')
gs = GridSpec(2, 1, figure=fig, height_ratios=[3, 1])

ax0 = fig.add_subplot(gs[0], projection='3d')
ax0.set_facecolor('#0f0f0f')

for planet in ['Earth', 'Mars', 'Jupiter']:
    ax0.plot(rk4[('X(m)', planet)], rk4[('Y(m)', planet)], rk4[('Z(m)', planet)],
            label=f'{planet} (RK4)', alpha=0.8)
    ax0.plot(rk5[('X(m)', planet)], rk5[('Y(m)', planet)], rk5[('Z(m)', planet)],
            '--', alpha=0.6, label=f'{planet} (RK5)')

ax0.set_title("Сравнение орбитальных траекторий", fontsize=16, pad=20)
ax0.set_xlabel('X (AU)', fontsize=12)
ax0.set_ylabel('Y (AU)', fontsize=12)
ax0.set_zlabel('Z (AU)', fontsize=12)
ax0.xaxis.pane.fill = False
ax0.yaxis.pane.fill = False
ax0.zaxis.pane.fill = False
ax0.grid(color='gray', alpha=0.2)
ax0.legend(loc='upper right')

# График ошибок
ax1 = fig.add_subplot(gs[1])
ax1.plot(time, instant_error_au, 'lime', label='Мгновенная ошибка')
ax1.plot(time, cumulative_error_au, 'cyan', linewidth=2, label='Накопленная ошибка')
ax1.set_title("Динамика ошибок", fontsize=14)
ax1.set_xlabel("Время (дни)", fontsize=12)
ax1.set_ylabel("Ошибка (AU)", fontsize=12)
ax1.grid(alpha=0.2)
ax1.legend()
ax1.set_facecolor('#0f0f0f')

# Добавляем вторую ось для накопленной ошибки
ax2 = ax1.twinx()
cumulative_error_km = cumulative_error / 1e3  # переводим в тысячи км
ax2.plot(time, cumulative_error_km, 'white', linestyle=':', alpha=0.5)
ax2.set_ylabel("Накопленная ошибка (тыс. км)", color='white', fontsize=12)

plt.tight_layout()
plt.subplots_adjust(hspace=0.25)
plt.show()