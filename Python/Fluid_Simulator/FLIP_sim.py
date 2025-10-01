"""
FLIP Fluid (adaptado)
"""
import sys
import math
import time
import numpy as np
import pygame
import traceback
import os


try:
    from numba import njit
    NUMBA_AVAILABLE = True
except Exception:
    def njit(func=None, **kwargs):
        def _dec(f):
            return f
        if func is None:
            return _dec
        return _dec(func)
    NUMBA_AVAILABLE = False

def clamp(a, lo, hi):
    return max(lo, min(hi, a))


# ----- core simulation functions -----
@njit
def integrate_particles_nb(particle_pos, particle_vel, num_particles, dt, gravity):
    for i in range(num_particles):
        particle_vel[2*i+1] += dt * gravity
    for i in range(num_particles):
        particle_pos[2*i] += particle_vel[2*i] * dt
        particle_pos[2*i+1] += particle_vel[2*i+1] * dt

@njit
def build_cell_lists_nb(particle_pos, pInvSpacing, pNumX, pNumY, num_particles, numCellParticles, firstCellParticle, cellParticleIds):
    for i in range(len(numCellParticles)):
        numCellParticles[i] = 0
    for i in range(num_particles):
        x = particle_pos[2*i]
        y = particle_pos[2*i+1]
        xi = int(x * pInvSpacing)
        if xi < 0: xi = 0
        if xi >= pNumX: xi = pNumX - 1
        yi = int(y * pInvSpacing)
        if yi < 0: yi = 0
        if yi >= pNumY: yi = pNumY - 1
        cellNr = xi * pNumY + yi
        numCellParticles[cellNr] += 1
    s = 0
    for i in range(pNumX * pNumY):
        s += numCellParticles[i]
        firstCellParticle[i] = s
    firstCellParticle[pNumX * pNumY] = s
    for i in range(num_particles):
        x = particle_pos[2*i]
        y = particle_pos[2*i+1]
        xi = int(x * pInvSpacing)
        if xi < 0: xi = 0
        if xi >= pNumX: xi = pNumX - 1
        yi = int(y * pInvSpacing)
        if yi < 0: yi = 0
        if yi >= pNumY: yi = pNumY - 1
        cellNr = xi * pNumY + yi
        firstCellParticle[cellNr] -= 1
        idx = firstCellParticle[cellNr]
        cellParticleIds[idx] = i

@njit
def push_particles_apart_nb(particle_pos, particle_vel, particle_color,
                            num_particles, pInvSpacing, pNumX, pNumY,
                            numCellParticles, firstCellParticle, cellParticleIds,
                            particleRadius, numIters, colorDiffusionCoeff):
    minDist = 2.0 * particleRadius
    minDist2 = minDist * minDist
    for it in range(numIters):
        for i in range(num_particles):
            px = particle_pos[2*i]; py = particle_pos[2*i+1]
            pxi = int(px * pInvSpacing)
            pyi = int(py * pInvSpacing)
            x0 = pxi - 1
            if x0 < 0: x0 = 0
            y0 = pyi - 1
            if y0 < 0: y0 = 0
            x1 = pxi + 1
            if x1 >= pNumX: x1 = pNumX - 1
            y1 = pyi + 1
            if y1 >= pNumY: y1 = pNumY - 1
            for xi in range(x0, x1+1):
                for yi in range(y0, y1+1):
                    cellNr = xi * pNumY + yi
                    first = firstCellParticle[cellNr]
                    last = firstCellParticle[cellNr+1]
                    for j in range(first, last):
                        id = cellParticleIds[j]
                        if id == i:
                            continue
                        qx = particle_pos[2*id]; qy = particle_pos[2*id+1]
                        dx = qx - px; dy = qy - py
                        d2 = dx*dx + dy*dy
                        if d2 > minDist2 or d2 == 0.0:
                            continue
                        d = math.sqrt(d2)
                        s = 0.5 * (minDist - d) / d
                        dxs = dx * s; dys = dy * s
                        particle_pos[2*i]   -= dxs
                        particle_pos[2*i+1] -= dys
                        particle_pos[2*id]  += dxs
                        particle_pos[2*id+1]+= dys
                        for k in range(3):
                            c0 = particle_color[3*i + k]
                            c1 = particle_color[3*id + k]
                            c = 0.5*(c0 + c1)
                            particle_color[3*i + k] = c0 + (c - c0) * colorDiffusionCoeff
                            particle_color[3*id + k] = c1 + (c - c1) * colorDiffusionCoeff

@njit
def handle_particle_collisions_nb(particle_pos, particle_vel, num_particles,
                                  fInvSpacing, fNumX, fNumY, particleRadius,
                                  obstacleX, obstacleY, obstacleRadius,
                                  scene_obstacleVelX, scene_obstacleVelY):
    h = 1.0 / fInvSpacing
    r = particleRadius
    orad = obstacleRadius
    minDist = orad + r
    minDist2 = minDist * minDist
    minX = h + r
    maxX = (fNumX - 1) * h - r
    minY = h + r
    maxY = (fNumY - 1) * h - r
    
    for i in range(num_particles):
        x = particle_pos[2*i]; y = particle_pos[2*i+1]
        
        # COLISÃO COM OBSTÁCULO - SISTEMA DE MÁSCARA APRIMORADO
        dx = x - obstacleX; dy = y - obstacleY
        d2 = dx*dx + dy*dy
        d = math.sqrt(d2 + 1e-12)
        
        # Aplica força de repulsão progressiva quando está próximo
        if d < minDist:
            # Normaliza o vetor de direção
            nx = dx / d; ny = dy / d
            
            # Calcula penetração
            penetration = minDist - d
            
            # Força de repulsão exponencial (mais forte quanto mais penetra)
            pushStrength = 0.5 + (penetration / minDist) * 0.5
            
            # Empurra partícula para fora
            particle_pos[2*i] += nx * penetration * pushStrength
            particle_pos[2*i+1] += ny * penetration * pushStrength
            
            # Componente normal da velocidade
            vn = particle_vel[2*i] * nx + particle_vel[2*i+1] * ny
            
            # Se está se aproximando do obstáculo
            if vn < 0.0:
                # Remove completamente a velocidade normal
                particle_vel[2*i] -= vn * nx
                particle_vel[2*i+1] -= vn * ny
                
                # Aplica fricção tangencial
                friction = 0.3
                particle_vel[2*i] *= (1.0 - friction)
                particle_vel[2*i+1] *= (1.0 - friction)
                
                # Adiciona pequena repulsão
                bounceStrength = 0.05
                particle_vel[2*i] += nx * bounceStrength
                particle_vel[2*i+1] += ny * bounceStrength
        
        # COLISÃO COM PAREDES
        x = particle_pos[2*i]
        y = particle_pos[2*i+1]
        
        if x < minX:
            x = minX
            if particle_vel[2*i] < 0.0:
                particle_vel[2*i] *= -0.3
        if x > maxX:
            x = maxX
            if particle_vel[2*i] > 0.0:
                particle_vel[2*i] *= -0.3
        if y < minY:
            y = minY
            if particle_vel[2*i+1] < 0.0:
                particle_vel[2*i+1] *= -0.3
        if y > maxY:
            y = maxY
            if particle_vel[2*i+1] > 0.0:
                particle_vel[2*i+1] *= -0.3
        particle_pos[2*i] = x; particle_pos[2*i+1] = y
@njit
def update_particle_density_nb(particle_pos, particleDensity, num_particles,
                               fNumX, fNumY, h, fInvSpacing):
    n = fNumY
    h2 = 0.5 * h
    for i in range(len(particleDensity)):
        particleDensity[i] = 0.0
    for i in range(num_particles):
        x = particle_pos[2*i]; y = particle_pos[2*i+1]
        if x < h: x = h
        if x > (fNumX - 1) * h: x = (fNumX - 1) * h
        if y < h: y = h
        if y > (fNumY - 1) * h: y = (fNumY - 1) * h
        x0 = int((x - h2) * fInvSpacing)
        tx = ((x - h2) - x0 * h) * fInvSpacing
        x1 = x0 + 1
        if x1 > fNumX-2:
            x1 = fNumX-2
        y0 = int((y - h2) * fInvSpacing)
        ty = ((y - h2) - y0 * h) * fInvSpacing
        y1 = y0 + 1
        if y1 > fNumY-2:
            y1 = fNumY-2
        sx = 1.0 - tx; sy = 1.0 - ty
        if 0 <= x0 < fNumX and 0 <= y0 < fNumY:
            particleDensity[x0*n + y0] += sx*sy
        if 0 <= x1 < fNumX and 0 <= y0 < fNumY:
            particleDensity[x1*n + y0] += tx*sy
        if 0 <= x1 < fNumX and 0 <= y1 < fNumY:
            particleDensity[x1*n + y1] += tx*ty
        if 0 <= x0 < fNumX and 0 <= y1 < fNumY:
            particleDensity[x0*n + y1] += sx*ty

@njit
def transfer_velocities_to_grid_nb(particle_pos, particle_vel, u, v, du, dv,
                                   num_particles, fNumX, fNumY, h, fInvSpacing, component_zero_is_u):
    n = fNumY
    h2 = 0.5 * h
    dx = 0.0 if component_zero_is_u else h2
    dy = h2 if component_zero_is_u else 0.0

    for p in range(num_particles):
        x = particle_pos[2*p]
        y = particle_pos[2*p+1]
        if x < h: x = h
        if x > (fNumX - 1) * h: x = (fNumX - 1) * h
        if y < h: y = h
        if y > (fNumY - 1) * h: y = (fNumY - 1) * h

        x0 = int((x - dx) * fInvSpacing)
        if x0 > fNumX - 2: x0 = fNumX - 2
        if x0 < 0: x0 = 0
        tx = ((x - dx) - x0 * h) * fInvSpacing
        x1 = x0 + 1

        y0 = int((y - dy) * fInvSpacing)
        if y0 > fNumY - 2: y0 = fNumY - 2
        if y0 < 0: y0 = 0
        ty = ((y - dy) - y0 * h) * fInvSpacing
        y1 = y0 + 1

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

        pv = particle_vel[2*p + (0 if component_zero_is_u else 1)]

        if component_zero_is_u:
            u[nr0] += pv * d0; du[nr0] += d0
            u[nr1] += pv * d1; du[nr1] += d1
            u[nr2] += pv * d2; du[nr2] += d2
            u[nr3] += pv * d3; du[nr3] += d3
        else:
            v[nr0] += pv * d0; dv[nr0] += d0
            v[nr1] += pv * d1; dv[nr1] += d1
            v[nr2] += pv * d2; dv[nr2] += d2
            v[nr3] += pv * d3; dv[nr3] += d3

@njit
def apply_obstacle_to_grid_nb(u, v, cellType, fNumX, fNumY, h, 
                              obstacleX, obstacleY, obstacleRadius):
    """Aplica influência do obstáculo no grid de velocidades"""
    n = fNumY
    orad2 = obstacleRadius * obstacleRadius
    
    for i in range(1, fNumX-1):
        for j in range(1, fNumY-1):
            idx = i*n + j
            
            # Posição do centro da célula
            x = (i + 0.5) * h
            y = (j + 0.5) * h
            
            dx = x - obstacleX
            dy = y - obstacleY
            d2 = dx*dx + dy*dy
            
            # Se a célula está próxima ou dentro do obstáculo
            if d2 < orad2 * 2.25:  # margem de 1.5x o raio
                if d2 > 0.0:
                    d = math.sqrt(d2)
                    # Marca como sólido se está muito próximo
                    if d < obstacleRadius * 1.1:
                        cellType[idx] = 2  # SOLID
                        # Zera velocidades nesta célula
                        u[idx] = 0.0
                        v[idx] = 0.0
                        u[(i+1)*n + j] = 0.0
                        v[i*n + (j+1)] = 0.0
                    # Se está na região de influência, reduz velocidades radiais
                    elif d < obstacleRadius * 1.5:
                        nx = dx / d
                        ny = dy / d
                        factor = (obstacleRadius * 1.5 - d) / (obstacleRadius * 0.5)
                        factor = max(0.0, min(1.0, factor)) * 0.7
                        
                        # Reduz componentes de velocidade que apontam para o obstáculo
                        vn_u = u[idx] * nx
                        vn_v = v[idx] * ny
                        if vn_u < 0: u[idx] *= (1.0 - factor)
                        if vn_v < 0: v[idx] *= (1.0 - factor)

@njit
def solve_incompressibility_nb(u, v, p, cellType, particleDensity,
                               fNumX, fNumY, density, h, dt,
                               numIters, overRelaxation, compensateDrift, particleRestDensity):
    n = fNumY
    cp = density * h / dt
    for i in range(len(p)):
        p[i] = 0.0
    
    for it in range(numIters):
        for i in range(1, fNumX-1):
            for j in range(1, fNumY-1):
                idx = i*n + j
                if cellType[idx] != 0:
                    continue
                left = (i-1)*n + j
                right= (i+1)*n + j
                bottom= i*n + j - 1
                top = i*n + j + 1
                sx0 = 1.0 if cellType[left] != 2 else 0.0
                sx1 = 1.0 if cellType[right]!= 2 else 0.0
                sy0 = 1.0 if cellType[bottom]!=2 else 0.0
                sy1 = 1.0 if cellType[top] !=2 else 0.0
                s = sx0 + sx1 + sy0 + sy1
                if s == 0.0:
                    continue
                div = u[right] - u[idx] + v[top] - v[idx]
                if particleRestDensity > 0.0 and compensateDrift:
                    compression = particleDensity[idx] - particleRestDensity
                    if compression > 0.0:
                        div = div - 0.8 * compression
                P = -div / s
                P = max(-0.5, min(0.5, P))
                P = P * overRelaxation
                p[idx] += cp * P
                u[idx] -= sx0 * P
                u[right] += sx1 * P
                v[idx] -= sy0 * P
                v[top] += sy1 * P


class FlipFluid:
    def __init__(self, density, width, height, spacing, particleRadius, maxParticles):
        self.density = density
        self.fNumX = int(math.floor(width / spacing)) + 1
        self.fNumY = int(math.floor(height / spacing)) + 1
        self.h = max(width / self.fNumX, height / self.fNumY)
        self.fInvSpacing = 1.0 / self.h
        self.fNumCells = self.fNumX * self.fNumY

        self.u = np.zeros(self.fNumCells, dtype=np.float32)
        self.v = np.zeros(self.fNumCells, dtype=np.float32)
        self.du = np.zeros(self.fNumCells, dtype=np.float32)
        self.dv = np.zeros(self.fNumCells, dtype=np.float32)
        self.prevU = np.zeros(self.fNumCells, dtype=np.float32)
        self.prevV = np.zeros(self.fNumCells, dtype=np.float32)
        self.p = np.zeros(self.fNumCells, dtype=np.float32)
        self.s = np.zeros(self.fNumCells, dtype=np.float32)
        self.cellType = np.ones(self.fNumCells, dtype=np.int32) * 1
        self.cellColor = np.zeros(3 * self.fNumCells, dtype=np.float32)

        self.maxParticles = maxParticles
        self.particlePos = np.zeros(2 * self.maxParticles, dtype=np.float32)
        self.particleVel = np.zeros(2 * self.maxParticles, dtype=np.float32)
        self.particleColor = np.zeros(3 * self.maxParticles, dtype=np.float32)
        self.particleColor[2::3] = 1.0
        self.particleDensity = np.zeros(self.fNumCells, dtype=np.float32)
        self.particleRestDensity = 0.0

        self.particleRadius = particleRadius
        self.pInvSpacing = 1.0 / (2.2 * particleRadius)
        self.pNumX = int(math.floor(width * self.pInvSpacing)) + 1
        self.pNumY = int(math.floor(height * self.pInvSpacing)) + 1
        self.pNumCells = self.pNumX * self.pNumY

        self.numCellParticles = np.zeros(self.pNumCells, dtype=np.int32)
        self.firstCellParticle = np.zeros(self.pNumCells + 1, dtype=np.int32)
        self.cellParticleIds = np.zeros(self.maxParticles, dtype=np.int32)

        self.numParticles = 0

    def integrate_particles(self, dt, gravity):
        if NUMBA_AVAILABLE:
            integrate_particles_nb(self.particlePos, self.particleVel, self.numParticles, dt, gravity)
        else:
            self.particleVel[1::2] += dt * gravity
            self.particlePos[0::2] += self.particleVel[0::2] * dt
            self.particlePos[1::2] += self.particleVel[1::2] * dt

        viscosity = 0.05
        damp = math.exp(-viscosity * dt)
        self.particleVel *= damp

        maxVel = 3.0
        vx = self.particleVel[0::2]
        vy = self.particleVel[1::2]
        speed = np.sqrt(vx * vx + vy * vy)
        if np.any(speed > maxVel):
            mask = speed > maxVel
            scale = maxVel / (speed[mask] + 1e-12)
            vx[mask] *= scale
            vy[mask] *= scale

    def push_particles_apart(self, numIters=2):
        if NUMBA_AVAILABLE:
            build_cell_lists_nb(self.particlePos, self.pInvSpacing, self.pNumX, self.pNumY,
                                self.numParticles, self.numCellParticles, self.firstCellParticle, self.cellParticleIds)
            push_particles_apart_nb(self.particlePos, self.particleVel, self.particleColor,
                                    self.numParticles, self.pInvSpacing, self.pNumX, self.pNumY,
                                    self.numCellParticles, self.firstCellParticle, self.cellParticleIds,
                                    self.particleRadius, numIters, 0.001)
        else:
            # fallback simplificado
            pass

    def handle_particle_collisions(self, obstacleX, obstacleY, obstacleRadius, scene_obstacleVelX, scene_obstacleVelY):
        handle_particle_collisions_nb(self.particlePos, self.particleVel, self.numParticles,
                                    self.fInvSpacing, self.fNumX, self.fNumY, self.particleRadius,
                                    obstacleX, obstacleY, obstacleRadius, scene_obstacleVelX, scene_obstacleVelY)
                
    def update_particle_density(self):
        update_particle_density_nb(self.particlePos, self.particleDensity, self.numParticles,
                                   self.fNumX, self.fNumY, self.h, self.fInvSpacing)
        if self.particleRestDensity == 0.0:
            s = 0.0; num = 0
            for i in range(self.fNumCells):
                if self.cellType[i] == 0:
                    s += self.particleDensity[i]; num += 1
            if num > 0:
                self.particleRestDensity = s / num

    def transfer_velocities_to_grid(self, obstacleX, obstacleY, obstacleRadius):
        self.prevU[:] = self.u
        self.prevV[:] = self.v
        self.du.fill(0.0)
        self.dv.fill(0.0)
        self.u.fill(0.0)
        self.v.fill(0.0)

        n = self.fNumY
        for i in range(self.fNumCells):
            self.cellType[i] = 2 if self.s[i] == 0.0 else 1

        for i in range(self.numParticles):
            x = self.particlePos[2*i]
            y = self.particlePos[2*i+1]
            xi = clamp(int(x * self.fInvSpacing), 0, self.fNumX - 1)
            yi = clamp(int(y * self.fInvSpacing), 0, self.fNumY - 1)
            cellNr = xi * n + yi
            if self.cellType[cellNr] == 1:
                self.cellType[cellNr] = 0

        transfer_velocities_to_grid_nb(self.particlePos, self.particleVel, self.u, self.v, self.du, self.dv,
                                    self.numParticles, self.fNumX, self.fNumY, self.h, self.fInvSpacing, True)
        transfer_velocities_to_grid_nb(self.particlePos, self.particleVel, self.u, self.v, self.du, self.dv,
                                    self.numParticles, self.fNumX, self.fNumY, self.h, self.fInvSpacing, False)

        for k in range(len(self.u)):
            if self.du[k] > 0.0:
                self.u[k] /= self.du[k]
        for k in range(len(self.v)):
            if self.dv[k] > 0.0:
                self.v[k] /= self.dv[k]

        # APLICA INFLUÊNCIA DO OBSTÁCULO NO GRID
        apply_obstacle_to_grid_nb(self.u, self.v, self.cellType, self.fNumX, self.fNumY, 
                                 self.h, obstacleX, obstacleY, obstacleRadius)

        for i in range(self.fNumX):
            for j in range(self.fNumY):
                idx = i * n + j
                solid = (self.cellType[idx] == 2)
                if solid or (i > 0 and self.cellType[(i-1)*n + j] == 2):
                    self.u[idx] = self.prevU[idx]
                if solid or (j > 0 and self.cellType[i*n + j - 1] == 2):
                    self.v[idx] = self.prevV[idx]

    def transfer_velocities_from_grid(self, flipRatio):
        n = self.fNumY
        h2 = 0.5 * self.h
        for i in range(self.numParticles):
            x = self.particlePos[2*i]; y = self.particlePos[2*i+1]
            x = clamp(x, self.h, (self.fNumX - 1) * self.h)
            y = clamp(y, self.h, (self.fNumY - 1) * self.h)
            for component in (0, 1):
                dx = 0.0 if component == 0 else h2
                dy = h2 if component == 0 else 0.0
                x0 = clamp(int((x - dx) * self.fInvSpacing), 0, self.fNumX-2)
                tx = ((x - dx) - x0*self.h) * self.fInvSpacing; x1 = x0 + 1
                y0 = clamp(int((y - dy) * self.fInvSpacing), 0, self.fNumY-2)
                ty = ((y - dy) - y0*self.h) * self.fInvSpacing; y1 = y0 + 1
                sx = 1.0 - tx; sy = 1.0 - ty
                d0 = sx*sy; d1 = tx*sy; d2 = tx*ty; d3 = sx*ty
                nr0 = x0*n + y0; nr1 = x1*n + y0; nr2 = x1*n + y1; nr3 = x0*n + y1

                valid0 = 1.0 if (self.cellType[nr0] != 1) else 0.0
                valid1 = 1.0 if (self.cellType[nr1] != 1) else 0.0
                valid2 = 1.0 if (self.cellType[nr2] != 1) else 0.0
                valid3 = 1.0 if (self.cellType[nr3] != 1) else 0.0

                d = valid0*d0 + valid1*d1 + valid2*d2 + valid3*d3
                if d > 0.0:
                    fvals = prevvals = 0.0
                    if valid0: 
                        fvals += d0 * (self.u[nr0] if component==0 else self.v[nr0])
                        prevvals += d0 * (self.prevU[nr0] if component==0 else self.prevV[nr0])
                    if valid1: 
                        fvals += d1 * (self.u[nr1] if component==0 else self.v[nr1])
                        prevvals += d1 * (self.prevU[nr1] if component==0 else self.prevV[nr1])
                    if valid2: 
                        fvals += d2 * (self.u[nr2] if component==0 else self.v[nr2])
                        prevvals += d2 * (self.prevU[nr2] if component==0 else self.prevV[nr2])
                    if valid3: 
                        fvals += d3 * (self.u[nr3] if component==0 else self.v[nr3])
                        prevvals += d3 * (self.prevU[nr3] if component==0 else self.prevV[nr3])

                    picV = fvals / d
                    corr = (fvals - prevvals) / d
                    flipV = self.particleVel[2*i + component] + corr
                    newV = (1.0 - flipRatio) * picV + flipRatio * flipV
                    self.particleVel[2*i + component] = newV

        maxVel = 8.0
        vx = self.particleVel[0::2]; vy = self.particleVel[1::2]
        speed = np.sqrt(vx*vx + vy*vy)
        if np.any(speed > maxVel):
            mask = speed > maxVel
            scale = maxVel / (speed[mask] + 1e-12)
            vx[mask] *= scale; vy[mask] *= scale

    def solve_incompressibility(self, numIters, dt, overRelaxation, compensateDrift=True):
        solve_incompressibility_nb(self.u, self.v, self.p, self.cellType, self.particleDensity,
                                   self.fNumX, self.fNumY, self.density, self.h, dt,
                                   numIters, overRelaxation, compensateDrift, self.particleRestDensity)

    def update_particle_colors(self):
        for i in range(self.numParticles):
            s = 0.01
            self.particleColor[3*i] = clamp(self.particleColor[3*i] - s, 0.0, 1.0)
            self.particleColor[3*i+1] = clamp(self.particleColor[3*i+1] - s, 0.0, 1.0)
            self.particleColor[3*i+2] = clamp(self.particleColor[3*i+2] + s, 0.0, 1.0)
            h1 = self.fInvSpacing
            x = self.particlePos[2*i]; y = self.particlePos[2*i+1]
            xi = clamp(int(x * h1), 1, self.fNumX-1)
            yi = clamp(int(y * h1), 1, self.fNumY-1)
            cellNr = xi * self.fNumY + yi
            d0 = self.particleRestDensity
            if d0 > 0.0:
                rel = self.particleDensity[cellNr] / d0
                if rel < 0.7:
                    s2 = 0.8
                    self.particleColor[3*i] = s2
                    self.particleColor[3*i+1] = s2
                    self.particleColor[3*i+2] = 1.0

    def update_cell_colors(self):
        self.cellColor.fill(0.0)
        for i in range(self.fNumCells):
            if self.cellType[i] == 2:
                self.cellColor[3*i] = 0.5; self.cellColor[3*i+1] = 0.5; self.cellColor[3*i+2] = 0.5
            elif self.cellType[i] == 0:
                d = self.particleDensity[i]
                if self.particleRestDensity > 0.0:
                    d /= self.particleRestDensity
                val = min(max((d - 0.0)/(2.0-0.0), 0.0), 0.9999)
                r = min(max(2.0*(val-0.5), 0.0),1.0)
                g = min(max(1.0 - abs(val-0.5)*2.0, 0.0),1.0)
                b = min(max(1.0 - val*2.0, 0.0),1.0)
                self.cellColor[3*i]=r; self.cellColor[3*i+1]=g; self.cellColor[3*i+2]=b

    def simulate(self, dt, gravity, flipRatio, numPressureIters, numParticleIters,
                 overRelaxation, compensateDrift, separateParticles,
                 obstacleX, obstacleY, obstacleRadius):
        self.integrate_particles(dt, gravity)
        if separateParticles and self.numParticles > 0:
            self.push_particles_apart(numParticleIters)
        self.handle_particle_collisions(obstacleX, obstacleY, obstacleRadius, 0.0, 0.0)
        self.transfer_velocities_to_grid(obstacleX, obstacleY, obstacleRadius)
        self.update_particle_density()
        self.solve_incompressibility(numPressureIters, dt, overRelaxation, compensateDrift)
        self.transfer_velocities_from_grid(flipRatio)
        self.update_particle_colors()
        self.update_cell_colors()


# ----- pygame renderer & scene setup -----

def world_to_screen(x, y, sim_w, sim_h, screen_w, screen_h):
    sx = int(x / sim_w * screen_w)
    sy = int(screen_h - (y / sim_h * screen_h))
    return sx, sy

def run_simulation(resolution=100, particle_scale=0.3, show_particles=True, max_fps=30):
    pygame.init()
    screen_w, screen_h = 800, 600
    screen = pygame.display.set_mode((screen_w, screen_h))
    clock = pygame.time.Clock()

    simHeight = 3.0
    cScale = screen_h / simHeight
    simWidth = screen_w / cScale

    tankHeight = 1.0 * simHeight
    tankWidth = 1.0 * simWidth
    h = tankHeight / resolution
    density = 500.0

    relWaterHeight = 0.6
    relWaterWidth = 0.5
    r = particle_scale * h
    dx = 2.5 * r
    dy = math.sqrt(3.0) / 2.0 * dx
    numX = int(math.floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx))
    numY = int(math.floor((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy))
    maxParticles = max(1, numX * numY)
    
    print(f"Número de partículas: {maxParticles} (numX={numX}, numY={numY})")

    fluid = FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles)
    fluid.numParticles = numX * numY
    p = 0
    for i in range(numX):
        for j in range(numY):
            fluid.particlePos[p] = h + r + dx * i + (0.0 if (j % 2 == 0) else r); p += 1
            fluid.particlePos[p] = h + r + dy * j; p += 1

    n = fluid.fNumY
    for i in range(fluid.fNumX):
        for j in range(fluid.fNumY):
            s = 1.0
            if i == 0 or i == fluid.fNumX-1 or j == 0:
                s = 0.0
            fluid.s[i*n + j] = s

    gravity = -9.81
    dt = 1.0 / 60.0
    flipRatio = 0.3
    numPressureIters = 60
    numParticleIters = 2
    overRelaxation = 1.5
    compensateDrift = True
    separateParticles = True
    obstacleX = 3.0
    obstacleY = 2.0
    obstacleRadius = 0.15

    running = True
    paused = False
    dragging = False
    font = pygame.font.SysFont("Arial", 16)

    while running:
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                running = False
            elif ev.type == pygame.KEYDOWN:
                if ev.key == pygame.K_p:
                    paused = not paused
                if ev.key == pygame.K_ESCAPE:
                    running = False
            elif ev.type == pygame.MOUSEBUTTONDOWN:
                if ev.button == 1:
                    mx, my = pygame.mouse.get_pos()
                    ox = mx / screen_w * simWidth
                    oy = (screen_h - my) / screen_h * simHeight
                    dx = ox - obstacleX
                    dy = oy - obstacleY
                    dist = math.sqrt(dx*dx + dy*dy)
                    if dist <= obstacleRadius * 1.5:
                        dragging = True
                    else:
                        obstacleX, obstacleY = ox, oy
            elif ev.type == pygame.MOUSEBUTTONUP:
                if ev.button == 1:
                    dragging = False
            elif ev.type == pygame.MOUSEMOTION:
                if dragging:
                    mx, my = pygame.mouse.get_pos()
                    obstacleX = mx / screen_w * simWidth
                    obstacleY = (screen_h - my) / screen_h * simHeight
                    obstacleX = clamp(obstacleX, obstacleRadius, tankWidth - obstacleRadius)
                    obstacleY = clamp(obstacleY, obstacleRadius, tankHeight - obstacleRadius)

        if not paused:
            fluid.simulate(dt, gravity, flipRatio, numPressureIters, numParticleIters,
                           overRelaxation, compensateDrift, separateParticles,
                           obstacleX, obstacleY, obstacleRadius)
        
        screen.fill((0,0,0))
        if show_particles:
            for i in range(fluid.numParticles):
                x = fluid.particlePos[2*i]; y = fluid.particlePos[2*i+1]
                sx, sy = world_to_screen(x, y, simWidth, simHeight, screen_w, screen_h)
                c = (int(255*fluid.particleColor[3*i]), int(255*fluid.particleColor[3*i+1]), int(255*fluid.particleColor[3*i+2]))
                pygame.draw.circle(screen, c, (sx, sy), max(1, int(1 + fluid.particleRadius / simWidth * screen_w)))
        
        ox_s, oy_s = world_to_screen(obstacleX, obstacleY, simWidth, simHeight, screen_w, screen_h)
        color = (255, 255, 0) if dragging else (255, 0, 0)
        pygame.draw.circle(screen, color, (ox_s, oy_s), max(2, int((obstacleRadius)/simWidth*screen_w)), 2)
        
        drag_text = " [ARRASTANDO]" if dragging else ""
        info = f"particles: {fluid.numParticles}  numba: {NUMBA_AVAILABLE}  fps: {int(clock.get_fps())}{drag_text}"
        txt = font.render(info, True, (255,255,255))
        screen.blit(txt, (10,10))
        
        pygame.display.flip()
        clock.tick(max_fps)

    pygame.quit()


if __name__ == "__main__":
    print("Iniciando FLIP Fluid - Colisão Melhorada")
    print(f"NumPy versão: {np.__version__}")
    print(f"Numba disponível: {NUMBA_AVAILABLE}")

    try:
        run_simulation(resolution=90, particle_scale=0.25, show_particles=True, max_fps=60)
    except Exception as e:
        tb = traceback.format_exc()
        print("Erro durante a execução:\n")
        print(tb)
        with open("error.log", "w", encoding="utf-8") as f:
            f.write("Erro durante a execução do flip_fast.py\n")
            f.write(tb)
        try:
            pygame.init()
            screen = pygame.display.set_mode((800, 200))
            pygame.display.set_caption("Erro - veja error.log")
            font = pygame.font.SysFont("Arial", 14)
            lines = tb.splitlines()[-8:]
            running = True
            while running:
                for ev in pygame.event.get():
                    if ev.type == pygame.QUIT or (ev.type == pygame.KEYDOWN and ev.key == pygame.K_ESCAPE):
                        running = False
                screen.fill((20, 20, 20))
                y = 10
                screen.blit(font.render("Erro durante a execução. Veja error.log", True, (255, 100, 100)), (10, y)); y += 22
                for ln in lines:
                    screen.blit(font.render(ln[:120], True, (220, 220, 220)), (10, y))
                    y += 18
                pygame.display.flip()
                pygame.time.wait(100)
            pygame.quit()
        except Exception:
            print("Falha ao exibir janela de erro.")