"""
FLIP Fluid Simulator em Python
Baseado no código original de Matthias Müller - Ten Minute Physics
Convertido para Python com PyGame e NumPy
"""

import pygame
import numpy as np
import math

# Constantes
U_FIELD = 0
V_FIELD = 1

FLUID_CELL = 0
AIR_CELL = 1
SOLID_CELL = 2

def clamp(x, min_val, max_val):
    return max(min_val, min(x, max_val))


class FlipFluid:
    def __init__(self, density, width, height, spacing, particle_radius, max_particles):
        # Configuração da grade de fluido
        self.density = density
        self.f_num_x = int(width / spacing) + 1
        self.f_num_y = int(height / spacing) + 1
        self.h = max(width / self.f_num_x, height / self.f_num_y)
        self.f_inv_spacing = 1.0 / self.h
        self.f_num_cells = self.f_num_x * self.f_num_y

        # Arrays da grade
        self.u = np.zeros(self.f_num_cells, dtype=np.float32)
        self.v = np.zeros(self.f_num_cells, dtype=np.float32)
        self.du = np.zeros(self.f_num_cells, dtype=np.float32)
        self.dv = np.zeros(self.f_num_cells, dtype=np.float32)
        self.prev_u = np.zeros(self.f_num_cells, dtype=np.float32)
        self.prev_v = np.zeros(self.f_num_cells, dtype=np.float32)
        self.p = np.zeros(self.f_num_cells, dtype=np.float32)
        self.s = np.zeros(self.f_num_cells, dtype=np.float32)
        self.cell_type = np.zeros(self.f_num_cells, dtype=np.int32)
        self.cell_color = np.zeros(3 * self.f_num_cells, dtype=np.float32)

        # Configuração de partículas
        self.max_particles = max_particles
        self.particle_pos = np.zeros(2 * max_particles, dtype=np.float32)
        self.particle_color = np.zeros(3 * max_particles, dtype=np.float32)
        self.particle_color[2::3] = 1.0  # Define canal azul como 1.0
        
        self.particle_vel = np.zeros(2 * max_particles, dtype=np.float32)
        self.particle_density = np.zeros(self.f_num_cells, dtype=np.float32)
        self.particle_rest_density = 0.0

        # Grade de partículas para detecção de colisão
        self.particle_radius = particle_radius
        self.p_inv_spacing = 1.0 / (2.2 * particle_radius)
        self.p_num_x = int(width * self.p_inv_spacing) + 1
        self.p_num_y = int(height * self.p_inv_spacing) + 1
        self.p_num_cells = self.p_num_x * self.p_num_y

        self.num_cell_particles = np.zeros(self.p_num_cells, dtype=np.int32)
        self.first_cell_particle = np.zeros(self.p_num_cells + 1, dtype=np.int32)
        self.cell_particle_ids = np.zeros(max_particles, dtype=np.int32)

        self.num_particles = 0

    def integrate_particles(self, dt, gravity):
        """Integra as velocidades e posições das partículas"""
        for i in range(self.num_particles):
            self.particle_vel[2*i + 1] += dt * gravity
            self.particle_pos[2*i] += self.particle_vel[2*i] * dt
            self.particle_pos[2*i + 1] += self.particle_vel[2*i + 1] * dt

    def push_particles_apart(self, num_iters):
        """Separa partículas que estão muito próximas"""
        color_diffusion_coeff = 0.001

        # Conta partículas por célula
        self.num_cell_particles.fill(0)
        
        for i in range(self.num_particles):
            x = self.particle_pos[2*i]
            y = self.particle_pos[2*i + 1]
            
            xi = clamp(int(x * self.p_inv_spacing), 0, self.p_num_x - 1)
            yi = clamp(int(y * self.p_inv_spacing), 0, self.p_num_y - 1)
            cell_nr = xi * self.p_num_y + yi
            self.num_cell_particles[cell_nr] += 1

        # Somas parciais
        first = 0
        for i in range(self.p_num_cells):
            first += self.num_cell_particles[i]
            self.first_cell_particle[i] = first
        self.first_cell_particle[self.p_num_cells] = first

        # Preenche partículas nas células
        for i in range(self.num_particles):
            x = self.particle_pos[2*i]
            y = self.particle_pos[2*i + 1]
            
            xi = clamp(int(x * self.p_inv_spacing), 0, self.p_num_x - 1)
            yi = clamp(int(y * self.p_inv_spacing), 0, self.p_num_y - 1)
            cell_nr = xi * self.p_num_y + yi
            self.first_cell_particle[cell_nr] -= 1
            self.cell_particle_ids[self.first_cell_particle[cell_nr]] = i

        # Separa partículas
        min_dist = 2.0 * self.particle_radius
        min_dist2 = min_dist * min_dist

        for _ in range(num_iters):
            for i in range(self.num_particles):
                px = self.particle_pos[2*i]
                py = self.particle_pos[2*i + 1]

                pxi = int(px * self.p_inv_spacing)
                pyi = int(py * self.p_inv_spacing)
                x0 = max(pxi - 1, 0)
                y0 = max(pyi - 1, 0)
                x1 = min(pxi + 1, self.p_num_x - 1)
                y1 = min(pyi + 1, self.p_num_y - 1)

                for xi in range(x0, x1 + 1):
                    for yi in range(y0, y1 + 1):
                        cell_nr = xi * self.p_num_y + yi
                        first = self.first_cell_particle[cell_nr]
                        last = self.first_cell_particle[cell_nr + 1]
                        
                        for j in range(first, last):
                            id_p = self.cell_particle_ids[j]
                            if id_p == i:
                                continue
                            
                            qx = self.particle_pos[2*id_p]
                            qy = self.particle_pos[2*id_p + 1]
                            
                            dx = qx - px
                            dy = qy - py
                            d2 = dx*dx + dy*dy
                            
                            if d2 > min_dist2 or d2 == 0.0:
                                continue
                            
                            d = math.sqrt(d2)
                            s = 0.5 * (min_dist - d) / d
                            dx *= s
                            dy *= s
                            
                            self.particle_pos[2*i] -= dx
                            self.particle_pos[2*i + 1] -= dy
                            self.particle_pos[2*id_p] += dx
                            self.particle_pos[2*id_p + 1] += dy

                            # Difunde cores
                            for k in range(3):
                                color0 = self.particle_color[3*i + k]
                                color1 = self.particle_color[3*id_p + k]
                                color = (color0 + color1) * 0.5
                                self.particle_color[3*i + k] = color0 + (color - color0) * color_diffusion_coeff
                                self.particle_color[3*id_p + k] = color1 + (color - color1) * color_diffusion_coeff

    def handle_particle_collisions(self, obstacle_x, obstacle_y, obstacle_w, obstacle_h):
        """Trata colisões de partículas com obstáculo retangular e paredes"""
        h = 1.0 / self.f_inv_spacing
        r = self.particle_radius

        # limites do tanque (centros das partículas não atravessam)
        min_x = h + r
        max_x = (self.f_num_x - 1) * h - r
        min_y = h + r
        max_y = (self.f_num_y - 1) * h - r

        # limites do retângulo (obstacle_x, obstacle_y são as coordenadas do centro do retângulo)
        rect_min_x = obstacle_x - obstacle_w * 0.5
        rect_max_x = obstacle_x + obstacle_w * 0.5
        rect_min_y = obstacle_y - obstacle_h * 0.5
        rect_max_y = obstacle_y + obstacle_h * 0.5

        for i in range(self.num_particles):
            x = self.particle_pos[2*i]
            y = self.particle_pos[2*i + 1]

            # Colisões com paredes (mantidos)
            if x < min_x:
                x = min_x
                self.particle_vel[2*i] = 0.0
            if x > max_x:
                x = max_x
                self.particle_vel[2*i] = 0.0
            if y < min_y:
                y = min_y
                self.particle_vel[2*i + 1] = 0.0
            if y > max_y:
                y = max_y
                self.particle_vel[2*i + 1] = 0.0

            # Colisão com o retângulo: se centro da partícula estiver dentro do retângulo,
            # empurra pela borda mais próxima.
            if rect_min_x <= x <= rect_max_x and rect_min_y <= y <= rect_max_y:
                left_dist   = abs(x - rect_min_x)
                right_dist  = abs(rect_max_x - x)
                top_dist    = abs(y - rect_min_y)
                bottom_dist = abs(rect_max_y - y)
                min_dist = min(left_dist, right_dist, top_dist, bottom_dist)

                if min_dist == left_dist:
                    x = rect_min_x - r
                    self.particle_vel[2*i] = 0.0
                elif min_dist == right_dist:
                    x = rect_max_x + r
                    self.particle_vel[2*i] = 0.0
                elif min_dist == top_dist:
                    y = rect_min_y - r
                    self.particle_vel[2*i + 1] = 0.0
                else:
                    y = rect_max_y + r
                    self.particle_vel[2*i + 1] = 0.0

            self.particle_pos[2*i] = x
            self.particle_pos[2*i + 1] = y

    def update_particle_density(self):
        """Atualiza a densidade de partículas na grade"""
        n = self.f_num_y
        h = self.h
        h1 = self.f_inv_spacing
        h2 = 0.5 * h

        self.particle_density.fill(0.0)

        for i in range(self.num_particles):
            x = clamp(self.particle_pos[2*i], h, (self.f_num_x - 1) * h)
            y = clamp(self.particle_pos[2*i + 1], h, (self.f_num_y - 1) * h)

            x0 = int((x - h2) * h1)
            tx = ((x - h2) - x0 * h) * h1
            x1 = min(x0 + 1, self.f_num_x - 2)
            
            y0 = int((y - h2) * h1)
            ty = ((y - h2) - y0 * h) * h1
            y1 = min(y0 + 1, self.f_num_y - 2)

            sx = 1.0 - tx
            sy = 1.0 - ty

            if x0 < self.f_num_x and y0 < self.f_num_y:
                self.particle_density[x0 * n + y0] += sx * sy
            if x1 < self.f_num_x and y0 < self.f_num_y:
                self.particle_density[x1 * n + y0] += tx * sy
            if x1 < self.f_num_x and y1 < self.f_num_y:
                self.particle_density[x1 * n + y1] += tx * ty
            if x0 < self.f_num_x and y1 < self.f_num_y:
                self.particle_density[x0 * n + y1] += sx * ty

        if self.particle_rest_density == 0.0:
            total = 0.0
            num_fluid_cells = 0

            for i in range(self.f_num_cells):
                if self.cell_type[i] == FLUID_CELL:
                    total += self.particle_density[i]
                    num_fluid_cells += 1

            if num_fluid_cells > 0:
                self.particle_rest_density = total / num_fluid_cells

    def transfer_velocities(self, to_grid, flip_ratio=0.9):
        """Transfere velocidades entre partículas e grade"""
        n = self.f_num_y
        h = self.h
        h1 = self.f_inv_spacing
        h2 = 0.5 * h

        if to_grid:
            self.prev_u[:] = self.u
            self.prev_v[:] = self.v

            self.du.fill(0.0)
            self.dv.fill(0.0)
            self.u.fill(0.0)
            self.v.fill(0.0)

            self.cell_type[:] = np.where(self.s == 0.0, SOLID_CELL, AIR_CELL)

            for i in range(self.num_particles):
                x = self.particle_pos[2*i]
                y = self.particle_pos[2*i + 1]
                xi = clamp(int(x * h1), 0, self.f_num_x - 1)
                yi = clamp(int(y * h1), 0, self.f_num_y - 1)
                cell_nr = xi * n + yi
                if self.cell_type[cell_nr] == AIR_CELL:
                    self.cell_type[cell_nr] = FLUID_CELL

        for component in range(2):
            dx = 0.0 if component == 0 else h2
            dy = h2 if component == 0 else 0.0

            f = self.u if component == 0 else self.v
            prev_f = self.prev_u if component == 0 else self.prev_v
            d = self.du if component == 0 else self.dv

            for i in range(self.num_particles):
                x = clamp(self.particle_pos[2*i], h, (self.f_num_x - 1) * h)
                y = clamp(self.particle_pos[2*i + 1], h, (self.f_num_y - 1) * h)

                x0 = min(int((x - dx) * h1), self.f_num_x - 2)
                tx = ((x - dx) - x0 * h) * h1
                x1 = min(x0 + 1, self.f_num_x - 2)
                
                y0 = min(int((y - dy) * h1), self.f_num_y - 2)
                ty = ((y - dy) - y0 * h) * h1
                y1 = min(y0 + 1, self.f_num_y - 2)

                sx = 1.0 - tx
                sy = 1.0 - ty

                d0 = sx * sy
                d1 = tx * sy
                d2 = tx * ty
                d3 = sx * ty

                nr0 = x0 * n + y0
                nr1 = x1 * n + y0
                nr2 = x1 * n + y1
                nr3 = x0 * n + y1

                if to_grid:
                    pv = self.particle_vel[2*i + component]
                    f[nr0] += pv * d0
                    d[nr0] += d0
                    f[nr1] += pv * d1
                    d[nr1] += d1
                    f[nr2] += pv * d2
                    d[nr2] += d2
                    f[nr3] += pv * d3
                    d[nr3] += d3
                else:
                    offset = n if component == 0 else 1
                    valid0 = 1.0 if (self.cell_type[nr0] != AIR_CELL or self.cell_type[nr0 - offset] != AIR_CELL) else 0.0
                    valid1 = 1.0 if (self.cell_type[nr1] != AIR_CELL or self.cell_type[nr1 - offset] != AIR_CELL) else 0.0
                    valid2 = 1.0 if (self.cell_type[nr2] != AIR_CELL or self.cell_type[nr2 - offset] != AIR_CELL) else 0.0
                    valid3 = 1.0 if (self.cell_type[nr3] != AIR_CELL or self.cell_type[nr3 - offset] != AIR_CELL) else 0.0

                    v = self.particle_vel[2*i + component]
                    d_sum = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3

                    if d_sum > 0.0:
                        pic_v = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + 
                                valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d_sum
                        corr = (valid0 * d0 * (f[nr0] - prev_f[nr0]) + 
                               valid1 * d1 * (f[nr1] - prev_f[nr1]) +
                               valid2 * d2 * (f[nr2] - prev_f[nr2]) + 
                               valid3 * d3 * (f[nr3] - prev_f[nr3])) / d_sum
                        flip_v = v + corr

                        self.particle_vel[2*i + component] = (1.0 - flip_ratio) * pic_v + flip_ratio * flip_v

            if to_grid:
                mask = d > 0.0
                f[mask] /= d[mask]

                # Restaura células sólidas
                for i in range(self.f_num_x):
                    for j in range(self.f_num_y):
                        solid = self.cell_type[i * n + j] == SOLID_CELL
                        if solid or (i > 0 and self.cell_type[(i-1) * n + j] == SOLID_CELL):
                            self.u[i * n + j] = self.prev_u[i * n + j]
                        if solid or (j > 0 and self.cell_type[i * n + j - 1] == SOLID_CELL):
                            self.v[i * n + j] = self.prev_v[i * n + j]

    def solve_incompressibility(self, num_iters, dt, over_relaxation, compensate_drift=True):
        """Resolve a incompressibilidade do fluido"""
        self.p.fill(0.0)
        self.prev_u[:] = self.u
        self.prev_v[:] = self.v

        n = self.f_num_y
        cp = self.density * self.h / dt

        for _ in range(num_iters):
            for i in range(1, self.f_num_x - 1):
                for j in range(1, self.f_num_y - 1):
                    if self.cell_type[i*n + j] != FLUID_CELL:
                        continue

                    center = i * n + j
                    left = (i - 1) * n + j
                    right = (i + 1) * n + j
                    bottom = i * n + j - 1
                    top = i * n + j + 1

                    sx0 = self.s[left]
                    sx1 = self.s[right]
                    sy0 = self.s[bottom]
                    sy1 = self.s[top]
                    s = sx0 + sx1 + sy0 + sy1
                    
                    if s == 0.0:
                        continue

                    div = self.u[right] - self.u[center] + self.v[top] - self.v[center]

                    if self.particle_rest_density > 0.0 and compensate_drift:
                        k = 1.0
                        compression = self.particle_density[i*n + j] - self.particle_rest_density
                        if compression > 0.0:
                            div = div - k * compression

                    p = -div / s
                    p *= over_relaxation
                    self.p[center] += cp * p

                    self.u[center] -= sx0 * p
                    self.u[right] += sx1 * p
                    self.v[center] -= sy0 * p
                    self.v[top] += sy1 * p

    def update_particle_colors(self):
        """Atualiza as cores das partículas baseado na densidade"""
        h1 = self.f_inv_spacing

        for i in range(self.num_particles):
            s = 0.01
            self.particle_color[3*i] = clamp(self.particle_color[3*i] - s, 0.0, 1.0)
            self.particle_color[3*i + 1] = clamp(self.particle_color[3*i + 1] - s, 0.0, 1.0)
            self.particle_color[3*i + 2] = clamp(self.particle_color[3*i + 2] + s, 0.0, 1.0)

            x = self.particle_pos[2*i]
            y = self.particle_pos[2*i + 1]
            xi = clamp(int(x * h1), 1, self.f_num_x - 1)
            yi = clamp(int(y * h1), 1, self.f_num_y - 1)
            cell_nr = xi * self.f_num_y + yi

            d0 = self.particle_rest_density

            if d0 > 0.0:
                rel_density = self.particle_density[cell_nr] / d0
                if rel_density < 0.7:
                    s = 0.8
                    self.particle_color[3*i] = s
                    self.particle_color[3*i + 1] = s
                    self.particle_color[3*i + 2] = 1.0

    def simulate(self, dt, gravity, flip_ratio, num_pressure_iters, num_particle_iters,
                over_relaxation, compensate_drift, separate_particles, 
                obstacle_x, obstacle_y, obstacle_w, obstacle_h):
        """Executa um passo da simulação"""
        num_sub_steps = 1
        sdt = dt / num_sub_steps

        for _ in range(num_sub_steps):
            self.integrate_particles(sdt, gravity)
            if separate_particles:
                self.push_particles_apart(num_particle_iters)
            # agora passa retângulo completo (x,y,w,h)
            self.handle_particle_collisions(obstacle_x, obstacle_y, obstacle_w, obstacle_h)
            self.transfer_velocities(True)
            self.update_particle_density()
            self.solve_incompressibility(num_pressure_iters, sdt, over_relaxation, compensate_drift)
            self.transfer_velocities(False, flip_ratio)

        self.update_particle_colors()


class Scene:
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.gravity = -9.81
        self.dt = 1.0 / 60.0
        self.flip_ratio = 0.9
        self.num_pressure_iters = 50
        self.num_particle_iters = 2
        self.over_relaxation = 1.9
        self.compensate_drift = True
        self.separate_particles = True

        # obstáculo retangular (centro + largura/altura)
        self.obstacle_x = width / 2.0
        self.obstacle_y = height / 2.0
        self.obstacle_w = 0.4    # largura (em unidades de simulação)
        self.obstacle_h = 0.25   # altura
        self.obstacle_dragging = False

        # comece despausado para ver movimento
        self.paused = False
        self.show_particles = True
        self.fluid = None

        self.setup_scene()

    def setup_scene(self):
        """Configura a cena inicial"""
        # resolução / grade da simulação (um equilíbrio: não muito alto para PCs fracos)
        res = 60
        tank_height = 1.0 * self.height
        tank_width = 1.0 * self.width
        h = tank_height / res
        density = 1000.0

        # quantidade de água (fração do tanque)
        rel_water_height = 0.5  # metade da altura
        rel_water_width = 0.4   # 40% da largura

        # raio das partículas (menor raio => mais partículas)
        r = 0.5 * h
        dx = 2.0 * r
        dy = math.sqrt(3.0) / 2.0 * dx

        # calcula quantos encaixam horizontal/verticalmente dentro da região de água
        num_x = int((rel_water_width * tank_width - 2.0 * h - 2.0 * r) / dx)
        num_y = int((rel_water_height * tank_height - 2.0 * h - 2.0 * r) / dy)
        if num_x < 1:
            num_x = 1
        if num_y < 1:
            num_y = 1

        max_particles = num_x * num_y

        # Cria o fluido
        self.fluid = FlipFluid(density, tank_width, tank_height, h, r, max_particles)

        # Colocando as partículas de modo que a primeira linha fique encostada no fundo
        # base_y é a altura do primeiro centro de partícula (encostado no fundo)
        base_y = h + r   # borda inferior definida como h + r em várias partes do código
        start_x = h + r   # margem esquerda
        p = 0

        for i in range(num_x):
            for j in range(num_y):
                # deslocamento em x e y; j=0 -> primeira linha (no fundo)
                x = start_x + dx * i + (0.0 if j % 2 == 0 else r)
                y = base_y + dy * j
                self.fluid.particle_pos[p] = x
                self.fluid.particle_pos[p + 1] = y
                # assegura velocidade inicial zero
                self.fluid.particle_vel[p] = 0.0
                self.fluid.particle_vel[p + 1] = 0.0
                p += 2

        self.fluid.num_particles = num_x * num_y

        # Configura células da grade para o tanque
        n = self.fluid.f_num_y
        for i in range(self.fluid.f_num_x):
            for j in range(self.fluid.f_num_y):
                s = 1.0  # fluido
                if i == 0 or i == self.fluid.f_num_x - 1 or j == 0:
                    s = 0.0  # sólido
                self.fluid.s[i * n + j] = s

        # posição inicial do obstáculo: center already set in __init__

    def set_obstacle(self, x, y, w=None, h=None):
        self.obstacle_x = x
        self.obstacle_y = y
        if w is not None:
            self.obstacle_w = w
        if h is not None:
            self.obstacle_h = h


def main():
    pygame.init()
    
    # Configurações da janela (menor para visualização leve)
    sim_height = 1.5
    window_height = 360
    c_scale = window_height / sim_height
    sim_width = 2.0
    
    # janela px (largura, altura)
    screen = pygame.display.set_mode((480, window_height))
    pygame.display.set_caption("FLIP Fluid Simulator")
    clock = pygame.time.Clock()

    # Cria a cena
    scene = Scene(sim_width, sim_height)

    running = True
    mouse_down = False
    drag_offset_x = 0.0
    drag_offset_y = 0.0

    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    scene.paused = not scene.paused
            elif event.type == pygame.MOUSEBUTTONDOWN:
                mouse_down = True
                mx, my = pygame.mouse.get_pos()
                # converte coordenadas da tela para coordenadas de simulação (repara no flip Y)
                sim_x = mx / c_scale
                sim_y = (window_height - my) / c_scale
                # verifica se clicou dentro do retângulo (usando centro + w/h)
                if (scene.obstacle_x - scene.obstacle_w * 0.5 <= sim_x <= scene.obstacle_x + scene.obstacle_w * 0.5 and
                    scene.obstacle_y - scene.obstacle_h * 0.5 <= sim_y <= scene.obstacle_y + scene.obstacle_h * 0.5):
                    scene.obstacle_dragging = True
                    drag_offset_x = scene.obstacle_x - sim_x
                    drag_offset_y = scene.obstacle_y - sim_y
                else:
                    # se não clicou no sólido, move o obstáculo para o ponto (comportamento anterior)
                    scene.obstacle_x = sim_x
                    scene.obstacle_y = sim_y
                    scene.paused = False
            elif event.type == pygame.MOUSEBUTTONUP:
                mouse_down = False
                scene.obstacle_dragging = False
            elif event.type == pygame.MOUSEMOTION and scene.obstacle_dragging:
                mx, my = pygame.mouse.get_pos()
                sim_x = mx / c_scale
                sim_y = (window_height - my) / c_scale
                scene.obstacle_x = sim_x + drag_offset_x
                scene.obstacle_y = sim_y + drag_offset_y

        # Simula
        if not scene.paused:
            scene.fluid.simulate(
                scene.dt, scene.gravity, scene.flip_ratio,
                scene.num_pressure_iters, scene.num_particle_iters,
                scene.over_relaxation, scene.compensate_drift,
                scene.separate_particles, scene.obstacle_x,
                scene.obstacle_y, scene.obstacle_w, scene.obstacle_h
            )

        # Desenha
        screen.fill((0, 0, 0))

        if scene.show_particles:
            for i in range(scene.fluid.num_particles):
                x = int(scene.fluid.particle_pos[2*i] * c_scale)
                y = int(window_height - scene.fluid.particle_pos[2*i + 1] * c_scale)
                rr = int(scene.fluid.particle_color[3*i] * 255)
                gg = int(scene.fluid.particle_color[3*i + 1] * 255)
                bb = int(scene.fluid.particle_color[3*i + 2] * 255)
                
                radius = max(2, int(scene.fluid.particle_radius * c_scale))
                pygame.draw.circle(screen, (rr, gg, bb), (x, y), radius)

        # Desenha obstáculo retangular (preenchido)
        ox = int((scene.obstacle_x - scene.obstacle_w * 0.5) * c_scale)
        oy = int(window_height - (scene.obstacle_y + scene.obstacle_h * 0.5) * c_scale)
        ow = int(scene.obstacle_w * c_scale)
        oh = int(scene.obstacle_h * c_scale)
        pygame.draw.rect(screen, (200, 0, 0), (ox, oy, ow, oh), 0)

        # Informações na tela
        font = pygame.font.Font(None, 24)
        info_text = f"Partículas: {scene.fluid.num_particles} | Espaço: Pausar/Retomar"
        text = font.render(info_text, True, (255, 255, 255))
        screen.blit(text, (10, 10))

        if scene.paused:
            pause_text = font.render("PAUSADO", True, (255, 255, 0))
            screen.blit(pause_text, (10, 40))

        pygame.display.flip()
        clock.tick(60)

    pygame.quit()


if __name__ == "__main__":
    main()
