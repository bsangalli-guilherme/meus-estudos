# -*- coding: utf-8 -*-
"""
Simulador de Fluido 2D com Lattice Boltzmann (LBM)
Versão otimizada para performance máxima com Numba,
reutilizando arrays para evitar alocação de memória em loop.
"""
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import numba

# --- Parâmetros da Simulação ---
Nx = 400
Ny = 100
tau = 0.53
Nt = 10000
plot_every = 5

# --- Modelo D2Q9 ---
NL = 9
OPPOSITE_INDICES = np.array([0, 5, 6, 7, 8, 1, 2, 3, 4], dtype=np.int64)
cxs = np.array([0, 0, 1, 1, 1, 0, -1, -1, -1], dtype=np.int64)
cys = np.array([0, 1, 1, 0, -1, -1, -1, 0, 1], dtype=np.int64)
weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36], dtype=np.float64)

# --- Obstáculo (cilindro) ---
Y, X = np.ogrid[0:Ny, 0:Nx]
cylinder = ((X - Nx // 4)**2 + (Y - Ny // 2)**2) < (13**2)
cylinder = cylinder.astype(np.bool_)  # assegura booleano

@numba.njit()
def run_simulation_steps_jit(steps, F, cylinder, cxs, cys, weights, OPPOSITE_INDICES, NL, tau, F_old, rho, ux, uy, usq):
    """
    Função principal da simulação, otimizada com Numba.
    Reutiliza os arrays F_old, rho, ux, uy, usq para evitar alocação de memória.
    """
    Ny_sim, Nx_sim = F.shape[0], F.shape[1]
    eps = 1e-12
    for _ in range(steps):
        # --- Streaming (lendo de cópia) ---
        F_old[:, :, :] = F[:, :, :] # Copia para o buffer pré-alocado
        for i in range(NL):
            cx = cxs[i]
            cy = cys[i]
            for y in range(Ny_sim):
                prev_y = (y - cy) % Ny_sim
                for x in range(Nx_sim):
                    prev_x = (x - cx) % Nx_sim
                    F[y, x, i] = F_old[prev_y, prev_x, i]

        # --- Bounce-back no obstáculo (SEM indexação booleana) ---
        for y in range(Ny_sim):
            for x in range(Nx_sim):
                if cylinder[y, x]:
                    tmp = np.empty(NL, dtype=np.float64)
                    for i in range(NL):
                        tmp[i] = F[y, x, i]
                    for i in range(NL):
                        opp = OPPOSITE_INDICES[i]
                        F[y, x, i] = tmp[opp]

        # --- Colisão (usando buffers pré-alocados) ---
        for y in range(Ny_sim):
            for x in range(Nx_sim):
                s = 0.0
                for i in range(NL):
                    s += F[y, x, i]
                rho[y, x] = s + eps

        for y in range(Ny_sim):
            for x in range(Nx_sim):
                s_x, s_y = 0.0, 0.0
                for i in range(NL):
                    fi = F[y, x, i]
                    s_x += fi * cxs[i]
                    s_y += fi * cys[i]
                ux[y, x] = s_x / rho[y, x]
                uy[y, x] = s_y / rho[y, x]

        for y in range(Ny_sim):
            for x in range(Nx_sim):
                if cylinder[y, x]:
                    ux[y, x] = 0.0
                    uy[y, x] = 0.0

        for y in range(Ny_sim):
            for x in range(Nx_sim):
                usq[y, x] = ux[y, x]*ux[y, x] + uy[y, x]*uy[y, x]

        for i in range(NL):
            w, cx, cy = weights[i], cxs[i], cys[i]
            for y in range(Ny_sim):
                for x in range(Nx_sim):
                    cu = 3.0 * (cx * ux[y, x] + cy * uy[y, x])
                    Feq = w * rho[y, x] * (1.0 + cu + 0.5 * cu * cu - 1.5 * usq[y, x])
                    F[y, x, i] += -(1.0 / tau) * (F[y, x, i] - Feq)

    # calcula macroscópicos finais (usando os mesmos buffers)
    for y in range(Ny_sim):
        for x in range(Nx_sim):
            s = 0.0
            for i in range(NL):
                s += F[y, x, i]
            rho[y, x] = s + eps

    for y in range(Ny_sim):
        for x in range(Nx_sim):
            s_x, s_y = 0.0, 0.0
            for i in range(NL):
                fi = F[y, x, i]
                s_x += fi * cxs[i]
                s_y += fi * cys[i]
            ux[y, x] = s_x / rho[y, x]
            uy[y, x] = s_y / rho[y, x]

    for y in range(Ny_sim):
        for x in range(Nx_sim):
            usq[y, x] = ux[y, x]*ux[y, x] + uy[y, x]*uy[y, x]

    return F, ux, uy, usq

# --- Condições Iniciais ---
F = np.ones((Ny, Nx, NL), dtype=np.float64) + 0.1 * np.random.randn(Ny, Nx, NL)
F[:, :, 3] = 2.8

# --- Pré-alocação dos Buffers de Memória ---
F_old = np.empty_like(F)
rho = np.empty((Ny, Nx), dtype=np.float64)
ux = np.empty((Ny, Nx), dtype=np.float64)
uy = np.empty((Ny, Nx), dtype=np.float64)
usq = np.empty((Ny, Nx), dtype=np.float64)

# --- Visualização ---
fig, ax = plt.subplots(figsize=(10, 3))
ax.axis('off')
im = ax.imshow(np.zeros((Ny, Nx)), origin='lower', cmap='plasma',
              interpolation='bilinear', vmin=0, vmax=0.2)
fig.colorbar(im, label='Velocidade')
fig.tight_layout()

def update(frame):
    global F
    # Passa os buffers pré-alocados para a função Numba
    new_F, _, _, usq_final = run_simulation_steps_jit(plot_every, F, cylinder, cxs, cys, weights, OPPOSITE_INDICES, NL, tau, F_old, rho, ux, uy, usq)
    F = new_F
    speed = np.sqrt(usq_final)
    speed_masked = np.ma.masked_array(speed, mask=cylinder)
    im.set_data(speed_masked)
    ax.set_title(f"Iteração = {frame * plot_every}")
    if frame % 10 == 0:
        print(f"Desenhando quadro {frame} (Iteração {frame * plot_every})")
    return [im]

num_frames = Nt // plot_every
ani = animation.FuncAnimation(fig, update, frames=num_frames, interval=1, blit=False)

plt.show()

