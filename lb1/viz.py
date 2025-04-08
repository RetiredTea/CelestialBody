import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from collections import deque
import matplotlib.patches as mpatches

# Чтение данных
df = pd.read_csv('trajectories3.csv')
times = df['Time(h)'].unique()

# Параметры визуализации
PLANET_SCALE_FACTOR = 1e4
STARFIELD_SIZE = 500
TRAIL_LENGTH = 500 #200

# Настройка параметров тел
body_properties = {
    'Sun': {'color': '#FFD700', 'size': 139270 / PLANET_SCALE_FACTOR, 'glow': True}, #1392700
    'Mercury': {'color': '#A0522D', 'size': 4879 / PLANET_SCALE_FACTOR}, #4879
    'Venus': {'color': '#DEB887', 'size': 12104 / PLANET_SCALE_FACTOR},
    'Earth': {'color': '#4169E1', 'size': 12756 / PLANET_SCALE_FACTOR},
    'Moon': {'color': '#D3D3D3', 'size': 3475 / PLANET_SCALE_FACTOR},
    'Mars': {'color': '#CD5C5C', 'size': 6792 / PLANET_SCALE_FACTOR},
    'Jupiter': {'color': '#DAA520', 'size': 142984 / PLANET_SCALE_FACTOR},
    'Saturn': {'color': '#F4A460', 'size': 120536 / PLANET_SCALE_FACTOR, 'ring': True},
    'Uranus': {'color': '#87CEEB', 'size': 51118 / PLANET_SCALE_FACTOR},
    'Neptune': {'color': '#1E90FF', 'size': 49528 / PLANET_SCALE_FACTOR},
    'Pluto': {'color': '#A9A9A9', 'size': 2376 / PLANET_SCALE_FACTOR}
}

np.random.seed(42)
stars = np.random.uniform(-5e12, 5e12, (STARFIELD_SIZE, 3))

fig = plt.figure(figsize=(16, 12), facecolor='black')
ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor('black')

ax.grid(False)
ax.xaxis.pane.fill = ax.yaxis.pane.fill = ax.zaxis.pane.fill = False
ax.xaxis.set_major_formatter(lambda x, pos: f"{x / 1e12:.1f} AU")
ax.yaxis.set_major_formatter(lambda x, pos: f"{x / 1e12:.1f} AU")
ax.zaxis.set_major_formatter(lambda x, pos: f"{x / 1e12:.1f} AU")

scats = {}
lines = {}
rings = {}
history = {body: deque(maxlen=TRAIL_LENGTH) for body in body_properties}
annotations = {}

for body, props in body_properties.items():
    scats[body] = ax.plot([], [], [], 'o',
                          markersize=props['size'],
                          color=props['color'],
                          markeredgecolor='white' if props.get('glow') else props['color'],
                          zorder=10)[0]

    # След траектории
    lines[body], = ax.plot([], [], [],
                           color=props['color'],
                           alpha=0.3,
                           linewidth=1.5)
    if props.get('ring'):
        theta = np.linspace(0, 2 * np.pi, 100)
        x_ring = np.cos(theta) * props['size'] * 2.5
        y_ring = np.sin(theta) * props['size'] * 2.5
        rings[body] = ax.plot(x_ring, y_ring, 0,
                              color=props['color'],
                              alpha=0.5,
                              linewidth=1)[0]

    # Аннотации
    annotations[body] = ax.text(0, 0, 0, body,
                                color=props['color'],
                                fontsize=8,
                                visible=False)

# Добавление звездного фона
ax.scatter(stars[:, 0], stars[:, 1], stars[:, 2],
           s=0.1,
           color='white',
           alpha=0.3,
           marker='*')
speed_multiplier = 20
current_frame = 0
paused = False
show_labels = True


def update(frame):
    global current_frame, speed_multiplier, paused, show_labels

    if not paused:
        current_frame = (current_frame + speed_multiplier) % len(times)
        current_time = times[int(current_frame)]

        data = df[df['Time(h)'] == current_time]

        for body in body_properties:
            body_data = data[data['Body'] == body]
            if not body_data.empty:
                x = body_data['X(m)'].values[0]
                y = body_data['Y(m)'].values[0]
                z = body_data['Z(m)'].values[0]
                scats[body].set_data([x], [y])
                scats[body].set_3d_properties([z])
                history[body].append((x, y, z))
                if len(history[body]) > 1:
                    x_hist = [p[0] for p in history[body]]
                    y_hist = [p[1] for p in history[body]]
                    z_hist = [p[2] for p in history[body]]
                    lines[body].set_data(x_hist, y_hist)
                    lines[body].set_3d_properties(z_hist)
                if show_labels:
                    annotations[body].set_position((x, y))
                    annotations[body].set_3d_properties(z, zdir='z')
                    annotations[body].set_visible(True)
                else:
                    annotations[body].set_visible(False)
        ax.set_title(
            f'BILLYARD\nTime: {current_time / 24:.1f} hours | Speed: {speed_multiplier}x | [SPACE] Pause | [±] Speed | [L] Labels',
            color='white',
            fontsize=12,
            y=1.0)

    return list(scats.values()) + list(lines.values()) + list(annotations.values())


def on_key(event):
    global speed_multiplier, paused, show_labels
    if event.key == '=':
        speed_multiplier += 5
    elif event.key == '-':
        speed_multiplier -= 5
    elif event.key == ' ':
        paused = not paused
    elif event.key == 'l':
        show_labels = not show_labels


fig.canvas.mpl_connect('key_press_event', on_key)

# Легенда
legend_elements = [mpatches.Patch(color=props['color'], label=body)
                   for body, props in body_properties.items()]
ax.legend(handles=legend_elements,
          loc='upper left',
          fontsize=8,
          facecolor='black',
          labelcolor='white')

# Создание анимации
ani = animation.FuncAnimation(
    fig,
    update,
    frames=len(times),
    interval=30,
    blit=True,
    repeat=True,
    cache_frame_data=False
)
plt.axis('off')
plt.show()