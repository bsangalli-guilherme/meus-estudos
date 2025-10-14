"""
FLIP Fluid - VERSÃO ULTRA OTIMIZADA + INTERFACE INTERATIVA
Otimizações adicionais + controle em tempo real de propriedades
"""
import sys
import math
import time
import numpy as np
import pygame
import traceback

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


# ----- CORE SIMULATION - ULTRA OTIMIZADO -----

@njit
def integrate_particles_soa(pos_x, pos_y, vel_x, vel_y, num_particles, dt, gravity, viscosity, maxVel):
    """Integração ultra-otimizada"""
    damp = math.exp(-viscosity * dt)
    maxVel2 = maxVel * maxVel
    
    for i in range(num_particles):
        vel_y[i] += dt * gravity
        vel_x[i] *= damp
        vel_y[i] *= damp
        
        vx = vel_x[i]
        vy = vel_y[i]
        speed2 = vx * vx + vy * vy
        
        if speed2 > maxVel2:
            scale = maxVel / math.sqrt(speed2 + 1e-12)
            vel_x[i] = vx * scale
            vel_y[i] = vy * scale
        
        pos_x[i] += vel_x[i] * dt
        pos_y[i] += vel_y[i] * dt


@njit
def build_cell_lists_soa(pos_x, pos_y, pInvSpacing, pNumX, pNumY, num_particles, 
                         numCellParticles, firstCellParticle, cellParticleIds):
    """Construção de listas otimizada"""
    for i in range(len(numCellParticles)):
        numCellParticles[i] = 0
    
    for i in range(num_particles):
        xi = int(pos_x[i] * pInvSpacing)
        xi = max(0, min(xi, pNumX - 1))
        yi = int(pos_y[i] * pInvSpacing)
        yi = max(0, min(yi, pNumY - 1))
        numCellParticles[xi * pNumY + yi] += 1
    
    s = 0
    for i in range(pNumX * pNumY):
        s += numCellParticles[i]
        firstCellParticle[i] = s
    firstCellParticle[pNumX * pNumY] = s
    
    for i in range(num_particles):
        xi = int(pos_x[i] * pInvSpacing)
        xi = max(0, min(xi, pNumX - 1))
        yi = int(pos_y[i] * pInvSpacing)
        yi = max(0, min(yi, pNumY - 1))
        cellNr = xi * pNumY + yi
        firstCellParticle[cellNr] -= 1
        cellParticleIds[firstCellParticle[cellNr]] = i


@njit
def push_particles_apart_soa(pos_x, pos_y, vel_x, vel_y, color_r, color_g, color_b,
                             num_particles, pInvSpacing, pNumX, pNumY,
                             numCellParticles, firstCellParticle, cellParticleIds,
                             particleRadius, numIters, colorDiffusionCoeff):
    """Separação de partículas ultra-otimizada"""
    minDist = 2.0 * particleRadius
    minDist2 = minDist * minDist
    
    for it in range(numIters):
        for i in range(num_particles):
            px = pos_x[i]
            py = pos_y[i]
            pxi = int(px * pInvSpacing)
            pyi = int(py * pInvSpacing)
            
            x0 = max(0, pxi - 1)
            y0 = max(0, pyi - 1)
            x1 = min(pNumX - 1, pxi + 1)
            y1 = min(pNumY - 1, pyi + 1)
            
            for xi in range(x0, x1+1):
                for yi in range(y0, y1+1):
                    cellNr = xi * pNumY + yi
                    first = firstCellParticle[cellNr]
                    last = firstCellParticle[cellNr+1]
                    
                    for j in range(first, last):
                        id = cellParticleIds[j]
                        if id == i:
                            continue
                        
                        dx = pos_x[id] - px
                        dy = pos_y[id] - py
                        d2 = dx*dx + dy*dy
                        
                        if d2 > minDist2 or d2 < 1e-12:
                            continue
                        
                        d = math.sqrt(d2)
                        s = 0.5 * (minDist - d) / d
                        
                        pos_x[i] -= dx * s
                        pos_y[i] -= dy * s
                        pos_x[id] += dx * s
                        pos_y[id] += dy * s
                        
                        # Difusão de cor otimizada
                        c_r = 0.5*(color_r[i] + color_r[id])
                        color_r[i] += (c_r - color_r[i]) * colorDiffusionCoeff
                        color_r[id] += (c_r - color_r[id]) * colorDiffusionCoeff
                        
                        c_g = 0.5*(color_g[i] + color_g[id])
                        color_g[i] += (c_g - color_g[i]) * colorDiffusionCoeff
                        color_g[id] += (c_g - color_g[id]) * colorDiffusionCoeff
                        
                        c_b = 0.5*(color_b[i] + color_b[id])
                        color_b[i] += (c_b - color_b[i]) * colorDiffusionCoeff
                        color_b[id] += (c_b - color_b[id]) * colorDiffusionCoeff


@njit
def handle_particle_collisions_soa(pos_x, pos_y, vel_x, vel_y, num_particles,
                                   fInvSpacing, fNumX, fNumY, particleRadius,
                                   obstacleX, obstacleY, obstacleRadius,
                                   scene_obstacleVelX, scene_obstacleVelY,
                                   obstacle_friction, obstacle_restitution):
    """Colisões otimizadas com parâmetros ajustáveis"""
    h = 1.0 / fInvSpacing
    r = particleRadius
    minDist = obstacleRadius + r
    minDist2 = minDist * minDist
    
    minX = h + r
    maxX = (fNumX - 1) * h - r
    minY = h + r
    maxY = (fNumY - 1) * h - r
    
    for i in range(num_particles):
        x = pos_x[i]
        y = pos_y[i]
        
        # Colisão com obstáculo
        dx = x - obstacleX
        dy = y - obstacleY
        d2 = dx*dx + dy*dy
        
        if d2 < minDist2:
            d = math.sqrt(d2 + 1e-12)
            nx = dx / d
            ny = dy / d
            penetration = minDist - d
            pushStrength = 0.5 + (penetration / minDist) * 0.5
            
            pos_x[i] += nx * penetration * pushStrength
            pos_y[i] += ny * penetration * pushStrength
            
            vn = vel_x[i] * nx + vel_y[i] * ny
            
            if vn < 0.0:
                vel_x[i] -= vn * nx
                vel_y[i] -= vn * ny
                
                vel_x[i] *= (1.0 - obstacle_friction)
                vel_y[i] *= (1.0 - obstacle_friction)
                
                vel_x[i] += nx * obstacle_restitution
                vel_y[i] += ny * obstacle_restitution
        
        # Colisão com paredes
        x = pos_x[i]
        y = pos_y[i]
        
        if x < minX:
            x = minX
            if vel_x[i] < 0.0:
                vel_x[i] *= -0.3
        if x > maxX:
            x = maxX
            if vel_x[i] > 0.0:
                vel_x[i] *= -0.3
        if y < minY:
            y = minY
            if vel_y[i] < 0.0:
                vel_y[i] *= -0.3
        if y > maxY:
            y = maxY
            if vel_y[i] > 0.0:
                vel_y[i] *= -0.3
        
        pos_x[i] = x
        pos_y[i] = y


@njit
def update_particle_density_soa(pos_x, pos_y, particleDensity, num_particles,
                                fNumX, fNumY, h, fInvSpacing):
    """Densidade otimizada"""
    n = fNumY
    h2 = 0.5 * h
    
    for i in range(len(particleDensity)):
        particleDensity[i] = 0.0
    
    for i in range(num_particles):
        x = max(h, min(pos_x[i], (fNumX - 1) * h))
        y = max(h, min(pos_y[i], (fNumY - 1) * h))
        
        x0 = int((x - h2) * fInvSpacing)
        tx = ((x - h2) - x0 * h) * fInvSpacing
        x1 = min(x0 + 1, fNumX - 2)
        
        y0 = int((y - h2) * fInvSpacing)
        ty = ((y - h2) - y0 * h) * fInvSpacing
        y1 = min(y0 + 1, fNumY - 2)
        
        sx = 1.0 - tx
        sy = 1.0 - ty
        
        if 0 <= x0 < fNumX and 0 <= y0 < fNumY:
            particleDensity[x0*n + y0] += sx*sy
        if 0 <= x1 < fNumX and 0 <= y0 < fNumY:
            particleDensity[x1*n + y0] += tx*sy
        if 0 <= x1 < fNumX and 0 <= y1 < fNumY:
            particleDensity[x1*n + y1] += tx*ty
        if 0 <= x0 < fNumX and 0 <= y1 < fNumY:
            particleDensity[x0*n + y1] += sx*ty


@njit
def transfer_velocities_to_grid_soa(pos_x, pos_y, vel_x, vel_y, u, v, du, dv,
                                   num_particles, fNumX, fNumY, h, fInvSpacing, 
                                   component_is_u):
    """Transferência para grid otimizada"""
    n = fNumY
    h2 = 0.5 * h
    dx = 0.0 if component_is_u else h2
    dy = h2 if component_is_u else 0.0

    for p in range(num_particles):
        x = max(h, min(pos_x[p], (fNumX - 1) * h))
        y = max(h, min(pos_y[p], (fNumY - 1) * h))

        x0 = max(0, min(int((x - dx) * fInvSpacing), fNumX - 2))
        tx = ((x - dx) - x0 * h) * fInvSpacing
        x1 = x0 + 1

        y0 = max(0, min(int((y - dy) * fInvSpacing), fNumY - 2))
        ty = ((y - dy) - y0 * h) * fInvSpacing
        y1 = y0 + 1

        sx = 1.0 - tx
        sy = 1.0 - ty
        
        nr0 = x0 * n + y0
        nr1 = x1 * n + y0
        nr2 = x1 * n + y1
        nr3 = x0 * n + y1

        pv = vel_x[p] if component_is_u else vel_y[p]
        d0 = sx * sy
        d1 = tx * sy
        d2 = tx * ty
        d3 = sx * ty

        if component_is_u:
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
def transfer_velocities_from_grid_soa(pos_x, pos_y, vel_x, vel_y, u, v, prevU, prevV,
                                     cellType, num_particles, fNumX, fNumY, h, 
                                     fInvSpacing, flipRatio, maxVel):
    """Transferência do grid otimizada"""
    n = fNumY
    h2 = 0.5 * h
    maxVel2 = maxVel * maxVel
    
    for i in range(num_particles):
        x = max(h, min(pos_x[i], (fNumX - 1) * h))
        y = max(h, min(pos_y[i], (fNumY - 1) * h))
        
        # Componente U
        x0 = max(0, min(int((x - 0.0) * fInvSpacing), fNumX-2))
        tx = ((x - 0.0) - x0*h) * fInvSpacing
        x1 = x0 + 1
        
        y0 = max(0, min(int((y - h2) * fInvSpacing), fNumY-2))
        ty = ((y - h2) - y0*h) * fInvSpacing
        y1 = y0 + 1
        
        sx = 1.0 - tx
        sy = 1.0 - ty
        
        nr0 = x0*n + y0
        nr1 = x1*n + y0
        nr2 = x1*n + y1
        nr3 = x0*n + y1

        valid0 = 1.0 if (cellType[nr0] != 1) else 0.0
        valid1 = 1.0 if (cellType[nr1] != 1) else 0.0
        valid2 = 1.0 if (cellType[nr2] != 1) else 0.0
        valid3 = 1.0 if (cellType[nr3] != 1) else 0.0

        d0 = sx*sy
        d1 = tx*sy
        d2 = tx*ty
        d3 = sx*ty
        d = valid0*d0 + valid1*d1 + valid2*d2 + valid3*d3
        
        if d > 0.0:
            fvals = prevvals = 0.0
            
            if valid0 > 0:
                fvals += d0 * u[nr0]
                prevvals += d0 * prevU[nr0]
            if valid1 > 0:
                fvals += d1 * u[nr1]
                prevvals += d1 * prevU[nr1]
            if valid2 > 0:
                fvals += d2 * u[nr2]
                prevvals += d2 * prevU[nr2]
            if valid3 > 0:
                fvals += d3 * u[nr3]
                prevvals += d3 * prevU[nr3]

            picV = fvals / d
            corr = (fvals - prevvals) / d
            flipV = vel_x[i] + corr
            vel_x[i] = (1.0 - flipRatio) * picV + flipRatio * flipV
        
        # Componente V
        x0 = max(0, min(int((x - h2) * fInvSpacing), fNumX-2))
        tx = ((x - h2) - x0*h) * fInvSpacing
        x1 = x0 + 1
        
        y0 = max(0, min(int((y - 0.0) * fInvSpacing), fNumY-2))
        ty = ((y - 0.0) - y0*h) * fInvSpacing
        y1 = y0 + 1
        
        sx = 1.0 - tx
        sy = 1.0 - ty
        
        nr0 = x0*n + y0
        nr1 = x1*n + y0
        nr2 = x1*n + y1
        nr3 = x0*n + y1

        valid0 = 1.0 if (cellType[nr0] != 1) else 0.0
        valid1 = 1.0 if (cellType[nr1] != 1) else 0.0
        valid2 = 1.0 if (cellType[nr2] != 1) else 0.0
        valid3 = 1.0 if (cellType[nr3] != 1) else 0.0

        d0 = sx*sy
        d1 = tx*sy
        d2 = tx*ty
        d3 = sx*ty
        d = valid0*d0 + valid1*d1 + valid2*d2 + valid3*d3
        
        if d > 0.0:
            fvals = prevvals = 0.0
            
            if valid0 > 0:
                fvals += d0 * v[nr0]
                prevvals += d0 * prevV[nr0]
            if valid1 > 0:
                fvals += d1 * v[nr1]
                prevvals += d1 * prevV[nr1]
            if valid2 > 0:
                fvals += d2 * v[nr2]
                prevvals += d2 * prevV[nr2]
            if valid3 > 0:
                fvals += d3 * v[nr3]
                prevvals += d3 * prevV[nr3]

            picV = fvals / d
            corr = (fvals - prevvals) / d
            flipV = vel_y[i] + corr
            vel_y[i] = (1.0 - flipRatio) * picV + flipRatio * flipV
        
        # Limita velocidade
        vx = vel_x[i]
        vy = vel_y[i]
        speed2 = vx*vx + vy*vy
        if speed2 > maxVel2:
            scale = maxVel / math.sqrt(speed2 + 1e-12)
            vel_x[i] = vx * scale
            vel_y[i] = vy * scale


@njit
def apply_obstacle_to_grid_nb(u, v, cellType, fNumX, fNumY, h, 
                              obstacleX, obstacleY, obstacleRadius):
    """Aplica obstáculo no grid"""
    n = fNumY
    orad2 = obstacleRadius * obstacleRadius
    search_radius2 = orad2 * 2.25
    
    for i in range(1, fNumX-1):
        for j in range(1, fNumY-1):
            idx = i*n + j
            x = (i + 0.5) * h
            y = (j + 0.5) * h
            
            dx = x - obstacleX
            dy = y - obstacleY
            d2 = dx*dx + dy*dy
            
            if d2 < search_radius2:
                if d2 > 1e-12:
                    d = math.sqrt(d2)
                    if d < obstacleRadius * 1.1:
                        cellType[idx] = 2
                        u[idx] = 0.0
                        v[idx] = 0.0
                        u[(i+1)*n + j] = 0.0
                        v[i*n + (j+1)] = 0.0
                    elif d < obstacleRadius * 1.5:
                        nx = dx / d
                        ny = dy / d
                        factor = (obstacleRadius * 1.5 - d) / (obstacleRadius * 0.5)
                        factor = max(0.0, min(1.0, factor)) * 0.7
                        
                        vn_u = u[idx] * nx
                        vn_v = v[idx] * ny
                        if vn_u < 0: u[idx] *= (1.0 - factor)
                        if vn_v < 0: v[idx] *= (1.0 - factor)


@njit
def solve_incompressibility_nb(u, v, p, cellType, particleDensity,
                               fNumX, fNumY, density, h, dt,
                               numIters, overRelaxation, compensateDrift, 
                               particleRestDensity):
    """Resolve incompressibilidade"""
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
                right = (i+1)*n + j
                bottom = i*n + j - 1
                top = i*n + j + 1
                
                sx0 = 1.0 if cellType[left] != 2 else 0.0
                sx1 = 1.0 if cellType[right] != 2 else 0.0
                sy0 = 1.0 if cellType[bottom] != 2 else 0.0
                sy1 = 1.0 if cellType[top] != 2 else 0.0
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


# ----- CLASSE PRINCIPAL -----

class FlipFluidSoA:
    """FLIP Fluid Simulator - Ultra Otimizado"""
    
    def __init__(self, density, width, height, spacing, particleRadius, maxParticles):
        self.density = density
        self.fNumX = int(math.floor(width / spacing)) + 1
        self.fNumY = int(math.floor(height / spacing)) + 1
        self.h = max(width / self.fNumX, height / self.fNumY)
        self.fInvSpacing = 1.0 / self.h
        self.fNumCells = self.fNumX * self.fNumY

        # Grid arrays
        self.u = np.zeros(self.fNumCells, dtype=np.float32)
        self.v = np.zeros(self.fNumCells, dtype=np.float32)
        self.du = np.zeros(self.fNumCells, dtype=np.float32)
        self.dv = np.zeros(self.fNumCells, dtype=np.float32)
        self.prevU = np.zeros(self.fNumCells, dtype=np.float32)
        self.prevV = np.zeros(self.fNumCells, dtype=np.float32)
        self.p = np.zeros(self.fNumCells, dtype=np.float32)
        self.s = np.zeros(self.fNumCells, dtype=np.float32)
        self.cellType = np.ones(self.fNumCells, dtype=np.int32)
        self.particleDensity = np.zeros(self.fNumCells, dtype=np.float32)

        # Particle arrays (SoA)
        self.maxParticles = maxParticles
        self.particle_pos_x = np.zeros(maxParticles, dtype=np.float32)
        self.particle_pos_y = np.zeros(maxParticles, dtype=np.float32)
        self.particle_vel_x = np.zeros(maxParticles, dtype=np.float32)
        self.particle_vel_y = np.zeros(maxParticles, dtype=np.float32)
        self.particle_color_r = np.zeros(maxParticles, dtype=np.float32)
        self.particle_color_g = np.zeros(maxParticles, dtype=np.float32)
        self.particle_color_b = np.ones(maxParticles, dtype=np.float32)
        
        self.particleRestDensity = 0.0
        self.particleRadius = particleRadius
        
        # Spatial hashing
        self.pInvSpacing = 1.0 / (2.2 * particleRadius)
        self.pNumX = int(math.floor(width * self.pInvSpacing)) + 1
        self.pNumY = int(math.floor(height * self.pInvSpacing)) + 1
        self.pNumCells = self.pNumX * self.pNumY

        self.numCellParticles = np.zeros(self.pNumCells, dtype=np.int32)
        self.firstCellParticle = np.zeros(self.pNumCells + 1, dtype=np.int32)
        self.cellParticleIds = np.zeros(self.maxParticles, dtype=np.int32)

        self.numParticles = 0

    def integrate_particles(self, dt, gravity, viscosity, maxVel):
        """Integra movimento das partículas"""
        integrate_particles_soa(self.particle_pos_x, self.particle_pos_y,
                               self.particle_vel_x, self.particle_vel_y,
                               self.numParticles, dt, gravity, viscosity, maxVel)

    def push_particles_apart(self, numIters):
        """Separa partículas"""
        if NUMBA_AVAILABLE and self.numParticles > 0:
            build_cell_lists_soa(self.particle_pos_x, self.particle_pos_y,
                                self.pInvSpacing, self.pNumX, self.pNumY,
                                self.numParticles, self.numCellParticles,
                                self.firstCellParticle, self.cellParticleIds)
            
            push_particles_apart_soa(self.particle_pos_x, self.particle_pos_y,
                                    self.particle_vel_x, self.particle_vel_y,
                                    self.particle_color_r, self.particle_color_g,
                                    self.particle_color_b, self.numParticles,
                                    self.pInvSpacing, self.pNumX, self.pNumY,
                                    self.numCellParticles, self.firstCellParticle,
                                    self.cellParticleIds, self.particleRadius,
                                    numIters, 0.001)

    def handle_particle_collisions(self, obstacleX, obstacleY, obstacleRadius, 
                                   friction, restitution):
        """Trata colisões"""
        handle_particle_collisions_soa(self.particle_pos_x, self.particle_pos_y,
                                      self.particle_vel_x, self.particle_vel_y,
                                      self.numParticles, self.fInvSpacing,
                                      self.fNumX, self.fNumY, self.particleRadius,
                                      obstacleX, obstacleY, obstacleRadius,
                                      0.0, 0.0, friction, restitution)

    def update_particle_density(self):
        """Atualiza densidade de partículas"""
        update_particle_density_soa(self.particle_pos_x, self.particle_pos_y,
                                   self.particleDensity, self.numParticles,
                                   self.fNumX, self.fNumY, self.h, self.fInvSpacing)
        
        if self.particleRestDensity == 0.0:
            fluid_cells = self.cellType == 0
            if np.any(fluid_cells):
                self.particleRestDensity = np.mean(self.particleDensity[fluid_cells])

    def transfer_velocities_to_grid(self, obstacleX, obstacleY, obstacleRadius):
        """Transfere velocidades para o grid"""
        self.prevU[:] = self.u
        self.prevV[:] = self.v
        self.du.fill(0.0)
        self.dv.fill(0.0)
        self.u.fill(0.0)
        self.v.fill(0.0)

        n = self.fNumY
        self.cellType[:] = np.where(self.s == 0.0, 2, 1)

        if self.numParticles > 0:
            xi = np.clip((self.particle_pos_x[:self.numParticles] * self.fInvSpacing).astype(np.int32), 
                        0, self.fNumX - 1)
            yi = np.clip((self.particle_pos_y[:self.numParticles] * self.fInvSpacing).astype(np.int32), 
                        0, self.fNumY - 1)
            cellNrs = xi * n + yi
            self.cellType[cellNrs] = np.where(self.cellType[cellNrs] == 1, 0, self.cellType[cellNrs])

        transfer_velocities_to_grid_soa(self.particle_pos_x, self.particle_pos_y,
                                       self.particle_vel_x, self.particle_vel_y,
                                       self.u, self.v, self.du, self.dv,
                                       self.numParticles, self.fNumX, self.fNumY,
                                       self.h, self.fInvSpacing, True)
        
        transfer_velocities_to_grid_soa(self.particle_pos_x, self.particle_pos_y,
                                       self.particle_vel_x, self.particle_vel_y,
                                       self.u, self.v, self.du, self.dv,
                                       self.numParticles, self.fNumX, self.fNumY,
                                       self.h, self.fInvSpacing, False)

        mask_u = self.du > 0.0
        self.u[mask_u] /= self.du[mask_u]
        mask_v = self.dv > 0.0
        self.v[mask_v] /= self.dv[mask_v]

        apply_obstacle_to_grid_nb(self.u, self.v, self.cellType, self.fNumX,
                                 self.fNumY, self.h, obstacleX, obstacleY,
                                 obstacleRadius)

        # OTIMIZAÇÃO: Loop único para células sólidas
        for i in range(self.fNumX):
            for j in range(self.fNumY):
                idx = i * n + j
                if self.cellType[idx] == 2:
                    self.u[idx] = self.prevU[idx]
                    self.v[idx] = self.prevV[idx]
                elif i > 0 and self.cellType[(i-1)*n + j] == 2:
                    self.u[idx] = self.prevU[idx]
                if j > 0 and self.cellType[i*n + j - 1] == 2:
                    self.v[idx] = self.prevV[idx]

    def transfer_velocities_from_grid(self, flipRatio, maxVel):
        """Transfere velocidades do grid"""
        transfer_velocities_from_grid_soa(
            self.particle_pos_x, self.particle_pos_y,
            self.particle_vel_x, self.particle_vel_y,
            self.u, self.v, self.prevU, self.prevV,
            self.cellType, self.numParticles,
            self.fNumX, self.fNumY, self.h, self.fInvSpacing,
            flipRatio, maxVel
        )

    def solve_incompressibility(self, numIters, dt, overRelaxation, compensateDrift):
        """Resolve incompressibilidade"""
        solve_incompressibility_nb(self.u, self.v, self.p, self.cellType,
                                   self.particleDensity, self.fNumX, self.fNumY,
                                   self.density, self.h, dt, numIters,
                                   overRelaxation, compensateDrift,
                                   self.particleRestDensity)

    def update_particle_colors(self, base_color):
        """Atualiza cores das partículas"""
        n = self.numParticles
        if n == 0:
            return
        
        # Aplica cor base
        self.particle_color_r[:n] = base_color[0]
        self.particle_color_g[:n] = base_color[1]
        self.particle_color_b[:n] = base_color[2]

    def simulate(self, dt, gravity, flipRatio, numPressureIters, numParticleIters,
                 overRelaxation, compensateDrift, separateParticles,
                 obstacleX, obstacleY, obstacleRadius, viscosity, maxVel,
                 obstacle_friction, obstacle_restitution, base_color):
        """Executa passo de simulação"""
        self.integrate_particles(dt, gravity, viscosity, maxVel)
        
        if separateParticles and self.numParticles > 0:
            self.push_particles_apart(numParticleIters)
        
        self.handle_particle_collisions(obstacleX, obstacleY, obstacleRadius,
                                       obstacle_friction, obstacle_restitution)
        self.transfer_velocities_to_grid(obstacleX, obstacleY, obstacleRadius)
        self.update_particle_density()
        self.solve_incompressibility(numPressureIters, dt, overRelaxation, compensateDrift)
        self.transfer_velocities_from_grid(flipRatio, maxVel)
        self.update_particle_colors(base_color)


# ----- CONTROLE DE INTERFACE -----

class UISlider:
    """Slider interativo para ajustar parâmetros"""
    def __init__(self, x, y, width, min_val, max_val, initial_val, label, format_str="{:.2f}"):
        self.x = x
        self.y = y
        self.width = width
        self.height = 20
        self.min_val = min_val
        self.max_val = max_val
        self.value = initial_val
        self.label = label
        self.format_str = format_str
        self.dragging = False
        self.rect = pygame.Rect(x, y, width, self.height)
        
    def handle_event(self, event):
        """Processa eventos do mouse"""
        if event.type == pygame.MOUSEBUTTONDOWN:
            if self.rect.collidepoint(event.pos):
                self.dragging = True
                self._update_value(event.pos[0])
        elif event.type == pygame.MOUSEBUTTONUP:
            self.dragging = False
        elif event.type == pygame.MOUSEMOTION:
            if self.dragging:
                self._update_value(event.pos[0])
    
    def _update_value(self, mouse_x):
        """Atualiza valor baseado na posição do mouse"""
        rel_x = max(0, min(mouse_x - self.x, self.width))
        ratio = rel_x / self.width
        self.value = self.min_val + ratio * (self.max_val - self.min_val)
    
    def draw(self, screen, font):
        """Desenha o slider"""
        # Barra de fundo
        pygame.draw.rect(screen, (60, 60, 60), self.rect)
        
        # Barra de progresso
        progress_width = int((self.value - self.min_val) / (self.max_val - self.min_val) * self.width)
        progress_rect = pygame.Rect(self.x, self.y, progress_width, self.height)
        pygame.draw.rect(screen, (100, 150, 255), progress_rect)
        
        # Borda
        pygame.draw.rect(screen, (200, 200, 200), self.rect, 2)
        
        # Label e valor
        text = f"{self.label}: {self.format_str.format(self.value)}"
        txt_surf = font.render(text, True, (255, 255, 255))
        screen.blit(txt_surf, (self.x, self.y - 20))


class UIButton:
    """Botão interativo"""
    def __init__(self, x, y, width, height, label, callback):
        self.rect = pygame.Rect(x, y, width, height)
        self.label = label
        self.callback = callback
        self.hovered = False
        
    def handle_event(self, event):
        """Processa eventos"""
        if event.type == pygame.MOUSEMOTION:
            self.hovered = self.rect.collidepoint(event.pos)
        elif event.type == pygame.MOUSEBUTTONDOWN:
            if self.hovered:
                self.callback()
    
    def draw(self, screen, font):
        """Desenha o botão"""
        color = (100, 150, 255) if self.hovered else (60, 90, 180)
        pygame.draw.rect(screen, color, self.rect)
        pygame.draw.rect(screen, (200, 200, 200), self.rect, 2)
        
        txt = font.render(self.label, True, (255, 255, 255))
        txt_rect = txt.get_rect(center=self.rect.center)
        screen.blit(txt, txt_rect)


class ControlPanel:
    """Painel de controle interativo"""
    def __init__(self, x, y, width):
        self.x = x
        self.y = y
        self.width = width
        self.visible = True
        self.sliders = []
        self.buttons = []
        self.current_y = y + 10
        
        # Presets de fluidos
        self.presets = {
            "Água": {
                "gravity": -9.81, "viscosity": 0.05, "density": 1000,
                "color": (0.2, 0.5, 1.0)
            },
            "Mel": {
                "gravity": -9.81, "viscosity": 2.5, "density": 1400,
                "color": (0.9, 0.7, 0.1)
            },
            "Óleo": {
                "gravity": -9.81, "viscosity": 0.8, "density": 900,
                "color": (0.3, 0.25, 0.1)
            },
            "Mercúrio": {
                "gravity": -9.81, "viscosity": 0.02, "density": 13600,
                "color": (0.7, 0.7, 0.8)
            },
            "Lava": {
                "gravity": -9.81, "viscosity": 5.0, "density": 3100,
                "color": (1.0, 0.3, 0.0)
            }
        }
    
    def add_slider(self, min_val, max_val, initial_val, label, format_str="{:.2f}"):
        """Adiciona um slider"""
        slider = UISlider(self.x + 10, self.current_y, self.width - 20,
                         min_val, max_val, initial_val, label, format_str)
        self.sliders.append(slider)
        self.current_y += 45
        return slider
    
    def add_button(self, label, callback):
        """Adiciona um botão"""
        button = UIButton(self.x + 10, self.current_y, self.width - 20, 30,
                         label, callback)
        self.buttons.append(button)
        self.current_y += 40
        return button
    
    def handle_event(self, event):
        """Processa eventos"""
        if not self.visible:
            return
        for slider in self.sliders:
            slider.handle_event(event)
        for button in self.buttons:
            button.handle_event(event)
    
    def draw(self, screen, font):
        """Desenha o painel"""
        if not self.visible:
            return
        
        # Fundo do painel
        panel_height = self.current_y - self.y + 10
        pygame.draw.rect(screen, (20, 20, 20, 200), 
                        (self.x, self.y, self.width, panel_height))
        pygame.draw.rect(screen, (100, 100, 100), 
                        (self.x, self.y, self.width, panel_height), 2)
        
        # Título
        title = font.render("⚙️ CONTROLES", True, (255, 255, 255))
        screen.blit(title, (self.x + 10, self.y + 10))
        
        # Desenha sliders e botões
        for slider in self.sliders:
            slider.draw(screen, font)
        for button in self.buttons:
            button.draw(screen, font)


def world_to_screen(x, y, sim_w, sim_h, screen_w, screen_h):
    """Converte coordenadas do mundo para a tela"""
    sx = int(x / sim_w * screen_w)
    sy = int(screen_h - (y / sim_h * screen_h))
    return sx, sy


# ----- SIMULAÇÃO PRINCIPAL -----

def run_simulation(resolution=90, particle_scale=0.25, max_fps=60):
    """Executa simulação com interface interativa"""
    pygame.init()
    screen_w, screen_h = 1280, 720
    screen = pygame.display.set_mode((screen_w, screen_h))
    pygame.display.set_caption("FLIP Fluid - Ultra Otimizado + Controles")
    clock = pygame.time.Clock()
    font = pygame.font.SysFont("Arial", 14)
    
    simHeight = 3
    cScale = screen_h / simHeight
    simWidth = screen_w / cScale
    
    tankHeight = 1.0 * simHeight
    tankWidth = 1.0 * simWidth
    h = tankHeight / resolution
    density = 1000
    
    relWaterHeight = 0.7
    relWaterWidth = 0.6
    r = particle_scale * h
    dx = 3 * r
    dy = math.sqrt(3.0) / 2.0 * dx
    numX = int(math.floor((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx))
    numY = int(math.floor((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy))
    maxParticles = max(1, numX * numY)
    
    print(f"╔{'═'*58}╗")
    print(f"║ FLIP Fluid - ULTRA OTIMIZADO + CONTROLES INTERATIVOS    ║")
    print(f"╚{'═'*58}╝")
    print(f"\n📊 Partículas: {maxParticles} | Resolução: {resolution}")
    print(f"⚡ Numba: {'✅' if NUMBA_AVAILABLE else '⚠️'}")
    
    # Cria simulador
    fluid = FlipFluidSoA(density, tankWidth, tankHeight, h, r, maxParticles)
    fluid.numParticles = numX * numY
    
    # Inicializa partículas
    p = 0
    for i in range(numX):
        for j in range(numY):
            fluid.particle_pos_x[p] = h + r + dx * i + (0.0 if (j % 2 == 0) else r)
            fluid.particle_pos_y[p] = h + r + dy * j
            p += 1
    
    # Configura paredes
    n = fluid.fNumY
    for i in range(fluid.fNumX):
        for j in range(fluid.fNumY):
            s = 1.0 if not (i == 0 or i == fluid.fNumX-1 or j == 0) else 0.0
            fluid.s[i*n + j] = s
    
    # Parâmetros iniciais
    params = {
        "gravity": -9.81,
        "viscosity": 0.05,
        "flipRatio": 0.3,
        "overRelaxation": 1.5,
        "obstacleRadius": 0.5,
        "obstacle_friction": 0.3,
        "obstacle_restitution": 0.05,
        "maxVel": 3.0,
        "color": [0.2, 0.5, 1.0]
    }
    
    dt = 1.0 / 60.0
    numPressureIters = 50
    numParticleIters = 1
    obstacleX = 3.0
    obstacleY = 2.0
    
    # Cria painel de controle
    panel = ControlPanel(10, 40, 280)
    
    # Adiciona sliders
    gravity_slider = panel.add_slider(-20, 0, params["gravity"], "Gravidade")
    viscosity_slider = panel.add_slider(0.01, 5.0, params["viscosity"], "Viscosidade")
    flip_slider = panel.add_slider(0.0, 1.0, params["flipRatio"], "FLIP Ratio")
    obstacle_size_slider = panel.add_slider(0.2, 1.5, params["obstacleRadius"], "Tamanho Obstáculo")
    friction_slider = panel.add_slider(0.0, 1.0, params["obstacle_friction"], "Atrito")
    restitution_slider = panel.add_slider(0.0, 0.5, params["obstacle_restitution"], "Repulsão")
    
    # Sliders de cor
    panel.current_y += 10
    color_r_slider = panel.add_slider(0.0, 1.0, params["color"][0], "Cor R")
    color_g_slider = panel.add_slider(0.0, 1.0, params["color"][1], "Cor G")
    color_b_slider = panel.add_slider(0.0, 1.0, params["color"][2], "Cor B")
    
    # Botões de preset
    panel.current_y += 10
    
    def apply_preset(preset_name):
        preset = panel.presets[preset_name]
        gravity_slider.value = preset["gravity"]
        viscosity_slider.value = preset["viscosity"]
        color_r_slider.value = preset["color"][0]
        color_g_slider.value = preset["color"][1]
        color_b_slider.value = preset["color"][2]
        print(f"✅ Preset aplicado: {preset_name}")
    
    panel.add_button("Água 💧", lambda: apply_preset("Água"))
    panel.add_button("Mel 🍯", lambda: apply_preset("Mel"))
    panel.add_button("Óleo 🛢️", lambda: apply_preset("Óleo"))
    panel.add_button("Mercúrio ☿", lambda: apply_preset("Mercúrio"))
    panel.add_button("Lava 🌋", lambda: apply_preset("Lava"))
    
    running = True
    paused = False
    dragging = False
    show_stats = True
    frame_times = []
    sim_times = []
    
    print(f"\n🎮 Controles:")
    print(f"  H - Mostrar/Ocultar painel")
    print(f"  P - Pausar/Continuar")
    print(f"  S - Stats")
    print(f"  Mouse - Arrastar obstáculo")
    print(f"  ESC - Sair\n")
    
    while running:
        frame_start = time.time()
        
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                running = False
            elif ev.type == pygame.KEYDOWN:
                if ev.key == pygame.K_h:
                    panel.visible = not panel.visible
                elif ev.key == pygame.K_p:
                    paused = not paused
                elif ev.key == pygame.K_s:
                    show_stats = not show_stats
                elif ev.key == pygame.K_ESCAPE:
                    running = False
            elif ev.type == pygame.MOUSEBUTTONDOWN and ev.button == 1:
                # Verifica se clicou no painel
                if not (panel.visible and ev.pos[0] < panel.width + 20):
                    mx, my = ev.pos
                    ox = mx / screen_w * simWidth
                    oy = (screen_h - my) / screen_h * simHeight
                    dx_click = ox - obstacleX
                    dy_click = oy - obstacleY
                    dist = math.sqrt(dx_click*dx_click + dy_click*dy_click)
                    if dist <= params["obstacleRadius"] * 1.5:
                        dragging = True
                    else:
                        obstacleX, obstacleY = ox, oy
            elif ev.type == pygame.MOUSEBUTTONUP and ev.button == 1:
                dragging = False
            elif ev.type == pygame.MOUSEMOTION:
                if dragging:
                    mx, my = ev.pos
                    obstacleX = mx / screen_w * simWidth
                    obstacleY = (screen_h - my) / screen_h * simHeight
                    obstacleX = np.clip(obstacleX, params["obstacleRadius"], 
                                       tankWidth - params["obstacleRadius"])
                    obstacleY = np.clip(obstacleY, params["obstacleRadius"], 
                                       tankHeight - params["obstacleRadius"])
            
            panel.handle_event(ev)
        
        # Atualiza parâmetros dos sliders
        params["gravity"] = gravity_slider.value
        params["viscosity"] = viscosity_slider.value
        params["flipRatio"] = flip_slider.value
        params["obstacleRadius"] = obstacle_size_slider.value
        params["obstacle_friction"] = friction_slider.value
        params["obstacle_restitution"] = restitution_slider.value
        params["color"] = [color_r_slider.value, color_g_slider.value, color_b_slider.value]
        
        # Simulação
        sim_start = time.time()
        if not paused:
            fluid.simulate(dt, params["gravity"], params["flipRatio"], 
                          numPressureIters, numParticleIters,
                          params["overRelaxation"], True, True,
                          obstacleX, obstacleY, params["obstacleRadius"],
                          params["viscosity"], params["maxVel"],
                          params["obstacle_friction"], params["obstacle_restitution"],
                          params["color"])
        sim_time = time.time() - sim_start
        
        # Renderização
        screen.fill((0, 0, 0))
        
        # Desenha partículas
        for i in range(fluid.numParticles):
            x = fluid.particle_pos_x[i]
            y = fluid.particle_pos_y[i]
            sx, sy = world_to_screen(x, y, simWidth, simHeight, screen_w, screen_h)
            
            r_val = int(255 * fluid.particle_color_r[i])
            g_val = int(255 * fluid.particle_color_g[i])
            b_val = int(255 * fluid.particle_color_b[i])
            c = (r_val, g_val, b_val)
            
            radius = max(1, int(1 + fluid.particleRadius / simWidth * screen_w))
            pygame.draw.circle(screen, c, (sx, sy), radius)
        
        # Desenha obstáculo
        ox_s, oy_s = world_to_screen(obstacleX, obstacleY, simWidth, simHeight, 
                                     screen_w, screen_h)
        color = (255, 255, 0) if dragging else (255, 100, 100)
        obstacle_radius = max(2, int(params["obstacleRadius"]/simWidth*screen_w))
        pygame.draw.circle(screen, color, (ox_s, oy_s), obstacle_radius, 3)
        
        # Painel de controle
        panel.draw(screen, font)
        
        # Estatísticas
        if show_stats:
            frame_time = time.time() - frame_start
            frame_times.append(frame_time)
            sim_times.append(sim_time)
            
            if len(frame_times) > 30:
                frame_times.pop(0)
                sim_times.pop(0)
            
            avg_frame = sum(frame_times) / len(frame_times)
            avg_sim = sum(sim_times) / len(sim_times)
            fps = 1.0 / avg_frame if avg_frame > 0 else 0
            
            status = " ⏸️" if paused else (" 🖱️" if dragging else "")
            
            stats_lines = [
                f"FPS: {int(clock.get_fps())} (avg: {fps:.1f})  |  Partículas: {fluid.numParticles}{status}",
                f"Sim: {avg_sim*1000:.1f}ms  |  Frame: {avg_frame*1000:.1f}ms"
            ]
            
            y_pos = screen_h - 60
            for line in stats_lines:
                txt = font.render(line, True, (255, 255, 255))
                bg_rect = txt.get_rect()
                bg_rect.bottomleft = (10, y_pos)
                bg_rect.inflate_ip(10, 4)
                pygame.draw.rect(screen, (0, 0, 0, 180), bg_rect)
                screen.blit(txt, (10, y_pos - 18))
                y_pos += 22
        
        pygame.display.flip()
        clock.tick(max_fps)
    
    pygame.quit()
    print("\n✅ Simulação encerrada.")


if __name__ == "__main__":
    print("=" * 60)
    print("FLIP Fluid Simulator")
    print("=" * 60)
    print(f"NumPy: {np.__version__} | Numba: {NUMBA_AVAILABLE}\n")
    
    try:
        run_simulation(resolution=90, particle_scale=0.15, max_fps=60)
    except Exception as e:
        print(f"\n❌ Erro: {e}")
        traceback.print_exc()