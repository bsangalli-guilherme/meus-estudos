# -*- coding: utf-8 -*-
"""
Simulador de Fluido Euleriano 2D Interativo
Baseado nas equações de Navier-Stokes.

Interação:
- Clique e arraste o mouse para "empurrar" o fluido e adicionar densidade.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- Parâmetros da Simulação ---
GRID_SIZE = 128 # Tamanho da grade (N x N). Use potências de 2 (64, 128) para melhor performance.
DT = 0.1              # Passo de tempo (delta t). Controla a estabilidade.
DIFFUSION = 0.0001    # Taxa de difusão (viscosidade). Um valor pequeno torna o fluido "espesso".
VISCOSITY = 0.0001    # Viscosidade do campo de velocidade.
SIM_STEPS = 50        # Passos de iteração do solver (Gauss-Seidel). Mais passos = mais precisão.


class FluidSim:
    def __init__(self, size, dt, diff, visc):
        self.size = size
        self.dt = dt
        self.diff = diff
        self.visc = visc

        # Campos de velocidade (u, v) e seus estados anteriores (u0, v0)
        self.u = np.zeros((size, size))
        self.v = np.zeros((size, size))
        self.u0 = np.zeros((size, size))
        self.v0 = np.zeros((size, size))

        # Densidade do fluido (ex: corante) e seu estado anterior
        self.density = np.zeros((size, size))
        self.density0 = np.zeros((size, size))
        
        # --- NOVO: Definição do Obstáculo ---
        Y, X = np.ogrid[0:size, 0:size]
        center_x, center_y = size // 2, size // 2
        radius = size // 8
        dist_sq = (X - center_x)**2 + (Y - center_y)**2
        self.obstacle = dist_sq < radius**2


    def add_source(self, x, y, amount_density, amount_u, amount_v):
        """ Adiciona uma fonte de densidade e velocidade em uma posição. """
        if 0 <= x < self.size and 0 <= y < self.size:
            self.density[x, y] += amount_density
            self.u[x, y] += amount_u
            self.v[x, y] += amount_v

    def step(self):
        """ Avança a simulação em um passo de tempo. """
        # 1. Difusão da velocidade (viscosidade)
        self.diffuse(1, self.u0, self.u, self.visc)
        self.diffuse(2, self.v0, self.v, self.visc)
        
        # 2. Projeção: garante que o fluido seja incompressível
        self.project(self.u0, self.v0, self.u, self.v)
        
        # 3. Advecção: move a velocidade através do próprio campo de velocidade
        self.advect(1, self.u, self.u0, self.u0, self.v0)
        self.advect(2, self.v, self.v0, self.u0, self.v0)

        # 4. Projeção novamente para limpar divergências da advecção
        self.project(self.u, self.v, self.u0, self.v0)

        # 5. Move a densidade (corante) pelo campo de velocidade final
        self.diffuse(0, self.density0, self.density, self.diff)
        self.advect(0, self.density, self.density0, self.u, self.v)

        # --- NOVO: Aplicar condição de contorno do obstáculo ---
        # Força a velocidade e a densidade a serem zero dentro do obstáculo.
        self.u[self.obstacle] = 0
        self.v[self.obstacle] = 0
        self.density[self.obstacle] = 0


    def set_bnd(self, b, x):
        """ Define as condições de contorno (paredes da caixa). """
        # Bordas esquerda/direita
        x[0, :] = -x[1, :] if b == 1 else x[1, :]
        x[self.size - 1, :] = -x[self.size - 2, :] if b == 1 else x[self.size - 2, :]
        
        # Bordas superior/inferior
        x[:, 0] = -x[:, 1] if b == 2 else x[:, 1]
        x[:, self.size - 1] = -x[:, self.size - 2] if b == 2 else x[:, self.size - 2]

        # Cantos
        x[0, 0] = 0.5 * (x[1, 0] + x[0, 1])
        x[0, self.size - 1] = 0.5 * (x[1, self.size - 1] + x[0, self.size - 2])
        x[self.size - 1, 0] = 0.5 * (x[self.size - 2, 0] + x[self.size - 1, 1])
        x[self.size - 1, self.size - 1] = 0.5 * (x[self.size - 2, self.size - 1] + x[self.size - 1, self.size - 2])

    def diffuse(self, b, x, x0, rate):
        """ Passo de Difusão: espalha a propriedade (velocidade ou densidade). """
        a = self.dt * rate * (self.size - 2) * (self.size - 2)
        for _ in range(SIM_STEPS):
            x[1:-1, 1:-1] = (x0[1:-1, 1:-1] + a * (x[:-2, 1:-1] + x[2:, 1:-1] + x[1:-1, :-2] + x[1:-1, 2:])) / (1 + 4 * a)
            self.set_bnd(b, x)

    def project(self, velocX, velocY, p, div):
        """ Passo de Projeção: torna o campo de velocidade 'mass-conserving'. """
        div[1:-1, 1:-1] = -0.5 * (velocX[2:, 1:-1] - velocX[:-2, 1:-1] + velocY[1:-1, 2:] - velocY[1:-1, :-2]) / self.size
        p.fill(0)
        self.set_bnd(0, div)
        self.set_bnd(0, p)

        for _ in range(SIM_STEPS):
            p[1:-1, 1:-1] = (div[1:-1, 1:-1] + p[:-2, 1:-1] + p[2:, 1:-1] + p[1:-1, :-2] + p[1:-1, 2:]) / 4
            self.set_bnd(0, p)

        velocX[1:-1, 1:-1] -= 0.5 * (p[2:, 1:-1] - p[:-2, 1:-1]) * self.size
        velocY[1:-1, 1:-1] -= 0.5 * (p[1:-1, 2:] - p[1:-1, :-2]) * self.size
        self.set_bnd(1, velocX)
        self.set_bnd(2, velocY)

    def advect(self, b, d, d0, velocX, velocY):
        """ Passo de Advecção: move a propriedade ao longo do campo de velocidade. """
        dtx = self.dt * (self.size - 2)
        dty = self.dt * (self.size - 2)

        i, j = np.meshgrid(np.arange(1, self.size - 1), np.arange(1, self.size - 1))
        
        tmp1 = dtx * velocX[i, j]
        tmp2 = dty * velocY[i, j]
        x = i - tmp1
        y = j - tmp2

        x = np.clip(x, 0.5, self.size - 1.5)
        y = np.clip(y, 0.5, self.size - 1.5)

        i0 = x.astype(int)
        i1 = i0 + 1
        j0 = y.astype(int)
        j1 = j0 + 1

        s1 = x - i0
        s0 = 1 - s1
        t1 = y - j0
        t0 = 1 - t1

        d[i, j] = (s0 * (t0 * d0[i0, j0] + t1 * d0[i0, j1]) +
                   s1 * (t0 * d0[i1, j0] + t1 * d0[i1, j1]))
        
        self.set_bnd(b, d)


# --- Configuração da Visualização e Interação ---
fig, ax = plt.subplots()
ax.axis('off')

fluid = FluidSim(GRID_SIZE, DT, DIFFUSION, VISCOSITY)
# Use um colormap diferente para melhor visualização
im = ax.imshow(fluid.density, cmap='magma', origin='lower', interpolation='bilinear')

prev_mouse_pos = [0, 0]

def on_mouse_move(event):
    global prev_mouse_pos
    if event.inaxes is None or not event.button:
        return
    
    ix, iy = int(event.xdata), int(event.ydata)
    
    # Adiciona densidade e velocidade com base no movimento do mouse
    dx = event.xdata - prev_mouse_pos[0]
    dy = event.ydata - prev_mouse_pos[1]
    
    # Adiciona uma "mancha" de densidade
    for x in range(ix-2, ix+3):
        for y in range(iy-2, iy+3):
            fluid.add_source(x, y, 100, dx * 5.0, dy * 5.0)

    prev_mouse_pos = [event.xdata, event.ydata]

def on_mouse_press(event):
    global prev_mouse_pos
    if event.inaxes is None: return
    prev_mouse_pos = [event.xdata, event.ydata]

fig.canvas.mpl_connect('motion_notify_event', on_mouse_move)
fig.canvas.mpl_connect('button_press_event', on_mouse_press)

def update(frame):
    fluid.step()
    # --- NOVO: Mascarar o obstáculo para a visualização ---
    masked_density = np.ma.masked_array(fluid.density, mask=fluid.obstacle)
    im.set_data(masked_density.T) # Transpor para corresponder à orientação do imshow
    im.set_clim(0, 255) # Define um limite de cor fixo para evitar cintilação
    return [im]

ani = animation.FuncAnimation(fig, update, frames=200, interval=1, blit=True)
plt.show()

