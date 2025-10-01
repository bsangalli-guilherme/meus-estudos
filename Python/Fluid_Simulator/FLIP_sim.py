"""
FLIP Fluid (adaptado e otimizado para Python)
- Versão otimizada para CPU (i5). Usa NumPy e, opcionalmente, Numba.
- Render com pygame (pontos). Ajuste parâmetros no bloco `if __name__ == "__main__":`
"""
import sys
import math
import time
import numpy as np
import pygame
import traceback
import os

# Try to import numba for speedups. If unavailable, code still runs (slower).
try:
    from numba import njit
    NUMBA_AVAILABLE = True
except Exception:
    def njit(func=None, **kwargs):
        # dummy decorator
        def _dec(f):
            return f
        if func is None:
            return _dec
        return _dec(func)
    NUMBA_AVAILABLE = False

# small helpers
def clamp(a, lo, hi):
    return max(lo, min(hi, a))


# ----- core simulation functions (numba-friendly signatures) -----
# We'll implement heavy loops as @njit functions that operate on plain numpy arrays.
# They accept arrays and size parameters, so Numba can compile them.

@njit
def integrate_particles_nb(particle_pos, particle_vel, num_particles, dt, gravity):
    for i in range(num_particles):
        particle_vel[2*i+1] += dt * gravity
    for i in range(num_particles):
        particle_pos[2*i] += particle_vel[2*i] * dt
        particle_pos[2*i+1] += particle_vel[2*i+1] * dt

@njit
def build_cell_lists_nb(particle_pos, pInvSpacing, pNumX, pNumY, num_particles, numCellParticles, firstCellParticle, cellParticleIds):
    # count
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
    # prefix sums into firstCellParticle (we store end indices then use decrement trick)
    s = 0
    for i in range(pNumX * pNumY):
        s += numCellParticles[i]
        firstCellParticle[i] = s
    firstCellParticle[pNumX * pNumY] = s
    # fill in reverse order
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
                        # diffuse colors
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
        
        # COLISÃO COM OBSTÁCULO
        dx = x - obstacleX; dy = y - obstacleY
        d2 = dx*dx + dy*dy
        if d2 < minDist2 and d2 > 0.0:
            d = math.sqrt(d2)
            pushFactor = (minDist - d) / d
            particle_pos[2*i] += dx * pushFactor
            particle_pos[2*i+1] += dy * pushFactor
            # reflexão com quique
            nx = dx / d; ny = dy / d
            vn = particle_vel[2*i] * nx + particle_vel[2*i+1] * ny
            particle_vel[2*i] -= 1.8 * vn * nx
            particle_vel[2*i+1] -= 1.8 * vn * ny
        
        # COLISÃO COM PAREDES
        if x < minX:
            x = minX; particle_vel[2*i] = 0.0
        if x > maxX:
            x = maxX; particle_vel[2*i] = 0.0
        if y < minY:
            y = minY; particle_vel[2*i+1] = 0.0
        if y > maxY:
            y = maxY; particle_vel[2*i+1] = 0.0
        particle_pos[2*i] = x; particle_pos[2*i+1] = y

@njit
def update_particle_density_nb(particle_pos, particleDensity, num_particles,
                               fNumX, fNumY, h, fInvSpacing):
    n = fNumY
    h2 = 0.5 * h
    # zero densities
    for i in range(len(particleDensity)):
        particleDensity[i] = 0.0
    for i in range(num_particles):
        x = particle_pos[2*i]; y = particle_pos[2*i+1]
        # clamp
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
    """
    Acumula a componente especificada (u se component_zero_is_u=True, senão v)
    nos arrays u/v e nos denominadores du/dv. Não zera os arrays: o caller deve zerá-los.
    """
    n = fNumY
    h2 = 0.5 * h

    dx = 0.0 if component_zero_is_u else h2
    dy = h2 if component_zero_is_u else 0.0

    for p in range(num_particles):
        x = particle_pos[2*p]
        y = particle_pos[2*p+1]

        # clamp para evitar índices inválidos
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
def finalize_grid_from_accum_nb(arr, denom):
    for i in range(len(arr)):
        if denom[i] > 0.0:
            arr[i] = arr[i] / denom[i]

@njit
def transfer_velocities_from_grid_nb(particle_pos, particle_vel, u, v, prevU, prevV,
                                     cellType, num_particles, fNumX, fNumY, h, fInvSpacing, flipRatio):
    # simplified: PIC/FLIP blend ignoring air-cell validity checks for speed
    n = fNumY
    h2 = 0.5 * h
    for p in range(num_particles):
        x = particle_pos[2*p]; y = particle_pos[2*p+1]
        # clamp
        if x < h: x = h
        if x > (fNumX - 1) * h: x = (fNumX - 1) * h
        if y < h: y = h
        if y > (fNumY - 1) * h: y = (fNumY - 1) * h
        # u-component (dx=0, dy=h2)
        for component in (0, 1):
            dx = 0.0 if component == 0 else h2
            dy = h2 if component == 0 else 0.0
            x0 = int((x - dx) * fInvSpacing)
            if x0 > fNumX-2: x0 = fNumX-2
            if x0 < 0: x0 = 0
            tx = ((x - dx) - x0*h) * fInvSpacing
            x1 = x0 + 1
            y0 = int((y - dy) * fInvSpacing)
            if y0 > fNumY-2: y0 = fNumY-2
            if y0 < 0: y0 = 0
            ty = ((y - dy) - y0*h) * fInvSpacing
            y1 = y0 + 1
            sx = 1.0 - tx; sy = 1.0 - ty
            d0 = sx*sy; d1 = tx*sy; d2 = tx*ty; d3 = sx*ty
            nr0 = x0*n + y0
            nr1 = x1*n + y0
            nr2 = x1*n + y1
            nr3 = x0*n + y1
            fvals = 0.0
            prevvals = 0.0
            denom = 0.0
            # naive validity: if cellType != 1 (air) treat as valid. This is a simplification for speed.
            if cellType[nr0] != 1: 
                fvals += d0 * (u[nr0] if component==0 else v[nr0]); prevvals += d0 * (prevU[nr0] if component==0 else prevV[nr0]); denom += d0
            if cellType[nr1] != 1:
                fvals += d1 * (u[nr1] if component==0 else v[nr1]); prevvals += d1 * (prevU[nr1] if component==0 else prevV[nr1]); denom += d1
            if cellType[nr2] != 1:
                fvals += d2 * (u[nr2] if component==0 else v[nr2]); prevvals += d2 * (prevU[nr2] if component==0 else prevV[nr2]); denom += d2
            if cellType[nr3] != 1:
                fvals += d3 * (u[nr3] if component==0 else v[nr3]); prevvals += d3 * (prevU[nr3] if component==0 else prevV[nr3]); denom += d3
            if denom > 0.0:
                picV = fvals / denom
                corr = (fvals - prevvals) / denom
                flipV = particle_vel[2*p + component] + corr
                particle_vel[2*p + component] = (1.0 - flipRatio) * picV + flipRatio * flipV

@njit
def solve_incompressibility_nb(u, v, p, cellType, particleDensity,
                               fNumX, fNumY, density, h, dt,
                               numIters, overRelaxation, compensateDrift, particleRestDensity):
    n = fNumY
    cp = density * h / dt
    # zero pressures
    for i in range(len(p)):
        p[i] = 0.0
    pcap = 0.5  # limite por iteração no p (ajuste para estabilidade)
    for it in range(numIters):
        for i in range(1, fNumX-1):
            for j in range(1, fNumY-1):
                idx = i*n + j
                if cellType[idx] != 0:  # 0 == FLUID_CELL
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
                        div = div - 1.0 * compression
                P = -div / s
                # cap P to avoid grandes saltos
                if P > pcap: P = pcap
                if P < -pcap: P = -pcap
                P = P * overRelaxation
                p[idx] += cp * P
                u[idx] -= sx0 * P
                u[right] += sx1 * P
                v[idx] -= sy0 * P
                v[top] += sy1 * P


# ----- high level FlipFluid class (holds arrays and calls njit functions) -----

class FlipFluid:
    def __init__(self, density, width, height, spacing, particleRadius, maxParticles):
        self.density = density
        self.fNumX = int(math.floor(width / spacing)) + 1
        self.fNumY = int(math.floor(height / spacing)) + 1
        self.h = max(width / self.fNumX, height / self.fNumY)
        self.fInvSpacing = 1.0 / self.h
        self.fNumCells = self.fNumX * self.fNumY

        # fluid fields (u and v defined on same cell centers as a simplification)
        self.u = np.zeros(self.fNumCells, dtype=np.float32)
        self.v = np.zeros(self.fNumCells, dtype=np.float32)
        self.du = np.zeros(self.fNumCells, dtype=np.float32)
        self.dv = np.zeros(self.fNumCells, dtype=np.float32)
        self.prevU = np.zeros(self.fNumCells, dtype=np.float32)
        self.prevV = np.zeros(self.fNumCells, dtype=np.float32)
        self.p = np.zeros(self.fNumCells, dtype=np.float32)
        self.s = np.zeros(self.fNumCells, dtype=np.float32)  # solid mask (1 fluid / 0 solid)
        # cellType: 0=FLUID,1=AIR,2=SOLID (we'll set later)
        self.cellType = np.ones(self.fNumCells, dtype=np.int32) * 1
        self.cellColor = np.zeros(3 * self.fNumCells, dtype=np.float32)

        # particles
        self.maxParticles = maxParticles
        self.particlePos = np.zeros(2 * self.maxParticles, dtype=np.float32)
        self.particleVel = np.zeros(2 * self.maxParticles, dtype=np.float32)
        self.particleColor = np.zeros(3 * self.maxParticles, dtype=np.float32)
        self.particleColor[2::3] = 1.0  # blue initial
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
        """
        Integra posição/velocidade das partículas, aplica amortecimento (viscosity),
        clamp de velocidade e colisão suave com paredes — método da classe FlipFluid.
        """
        # integração (com Numba quando disponível)
        if NUMBA_AVAILABLE:
            integrate_particles_nb(self.particlePos, self.particleVel, self.numParticles, dt, gravity)
        else:
            # fallback (vectorizado)
            self.particleVel[1::2] += dt * gravity
            self.particlePos[0::2] += self.particleVel[0::2] * dt
            self.particlePos[1::2] += self.particleVel[1::2] * dt

        # --- damping (viscosity) e clamp de velocidade para estabilidade ---
        viscosity = 2.0   # ajuste entre 0.0 (sem amortecimento) e ~5.0 (muito amortecido)
        damp = math.exp(-viscosity * dt)
        self.particleVel[0::2] *= damp
        self.particleVel[1::2] *= damp

        # clamp de velocidade (magnitude)
        maxVel = 6.0      # limite razoável; diminua se ainda explodir
        vx = self.particleVel[0::2]
        vy = self.particleVel[1::2]
        speed = np.sqrt(vx * vx + vy * vy)
        if np.any(speed > maxVel):
            mask = speed > maxVel
            scale = maxVel / (speed[mask] + 1e-12)
            vx[mask] *= scale
            vy[mask] *= scale

        # garante partículas dentro do tanque (colisão suave com paredes)
        h = 1.0 / self.fInvSpacing
        minX = h + self.particleRadius
        maxX = (self.fNumX - 1) * h - self.particleRadius
        minY = h + self.particleRadius
        maxY = (self.fNumY - 1) * h - self.particleRadius

        for i in range(self.numParticles):
            x = self.particlePos[2*i]; y = self.particlePos[2*i+1]
            if x < minX:
                self.particlePos[2*i] = minX
                if self.particleVel[2*i] < 0.0:
                    self.particleVel[2*i] = 0.0
            if x > maxX:
                self.particlePos[2*i] = maxX
                if self.particleVel[2*i] > 0.0:
                    self.particleVel[2*i] = 0.0
            if y < minY:
                self.particlePos[2*i+1] = minY
                if self.particleVel[2*i+1] < 0.0:
                    self.particleVel[2*i+1] = 0.0
            if y > maxY:
                self.particlePos[2*i+1] = maxY
                if self.particleVel[2*i+1] > 0.0:
                    self.particleVel[2*i+1] = 0.0


    def push_particles_apart(self, numIters=2):
        # build cell lists then call JIT loop
        if NUMBA_AVAILABLE:
            build_cell_lists_nb(self.particlePos, self.pInvSpacing, self.pNumX, self.pNumY,
                                self.numParticles, self.numCellParticles, self.firstCellParticle, self.cellParticleIds)
            push_particles_apart_nb(self.particlePos, self.particleVel, self.particleColor,
                                    self.numParticles, self.pInvSpacing, self.pNumX, self.pNumY,
                                    self.numCellParticles, self.firstCellParticle, self.cellParticleIds,
                                    self.particleRadius, numIters, 0.001)
        else:
            # fallback (less optimized but functional)
            self.numCellParticles.fill(0)
            for i in range(self.numParticles):
                x = self.particlePos[2*i]; y = self.particlePos[2*i+1]
                xi = int(x * self.pInvSpacing); yi = int(y * self.pInvSpacing)
                xi = max(0, min(self.pNumX-1, xi))
                yi = max(0, min(self.pNumY-1, yi))
                self.numCellParticles[xi * self.pNumY + yi] += 1
            first = 0
            for i in range(self.pNumCells):
                first += self.numCellParticles[i]
                self.firstCellParticle[i] = first
            self.firstCellParticle[self.pNumCells] = first
            for i in range(self.numParticles):
                x = self.particlePos[2*i]; y = self.particlePos[2*i+1]
                xi = int(x * self.pInvSpacing); yi = int(y * self.pInvSpacing)
                xi = max(0, min(self.pNumX-1, xi))
                yi = max(0, min(self.pNumY-1, yi))
                cellNr = xi * self.pNumY + yi
                self.firstCellParticle[cellNr] -= 1
                self.cellParticleIds[self.firstCellParticle[cellNr]] = i
            # pairwise push (naive)
            minDist = 2.0 * self.particleRadius
            minDist2 = minDist*minDist
            for _ in range(numIters):
                for i in range(self.numParticles):
                    px = self.particlePos[2*i]; py = self.particlePos[2*i+1]
                    pxi = int(px * self.pInvSpacing); pyi = int(py * self.pInvSpacing)
                    x0 = max(pxi-1, 0); y0 = max(pyi-1, 0)
                    x1 = min(pxi+1, self.pNumX-1); y1 = min(pyi+1, self.pNumY-1)
                    for xi in range(x0, x1+1):
                        for yi in range(y0, y1+1):
                            cellNr = xi * self.pNumY + yi
                            first = self.firstCellParticle[cellNr]
                            last = self.firstCellParticle[cellNr+1]
                            for j in range(first, last):
                                id = self.cellParticleIds[j]
                                if id == i: continue
                                qx = self.particlePos[2*id]; qy = self.particlePos[2*id+1]
                                dx = qx - px; dy = qy - py
                                d2 = dx*dx + dy*dy
                                if d2 > minDist2 or d2 == 0.0: continue
                                d = math.sqrt(d2)
                                s = 0.5 * (minDist - d) / d
                                dxs = dx * s; dys = dy * s
                                self.particlePos[2*i]   -= dxs
                                self.particlePos[2*i+1] -= dys
                                self.particlePos[2*id]  += dxs
                                self.particlePos[2*id+1]+= dys
                                # diffuse color slightly
                                for k in range(3):
                                    c0 = self.particleColor[3*i + k]
                                    c1 = self.particleColor[3*id + k]
                                    c = 0.5*(c0 + c1)
                                    self.particleColor[3*i + k] = c0 + (c - c0) * 0.001
                                    self.particleColor[3*id + k] = c1 + (c - c1) * 0.001

    def handle_particle_collisions(self, obstacleX, obstacleY, obstacleRadius, scene_obstacleVelX, scene_obstacleVelY):
        if NUMBA_AVAILABLE:
            handle_particle_collisions_nb(self.particlePos, self.particleVel, self.numParticles,
                                        self.fInvSpacing, self.fNumX, self.fNumY, self.particleRadius,
                                        obstacleX, obstacleY, obstacleRadius, scene_obstacleVelX, scene_obstacleVelY)
        else:
            h = 1.0 / self.fInvSpacing
            r = self.particleRadius
            orad = obstacleRadius
            minDist = orad + r
            minDist2 = minDist * minDist
            minX = h + r; maxX = (self.fNumX - 1) * h - r
            minY = h + r; maxY = (self.fNumY - 1) * h - r
            
            for i in range(self.numParticles):
                x = self.particlePos[2*i]; y = self.particlePos[2*i+1]
                
                # COLISÃO COM OBSTÁCULO (círculo vermelho)
                dx = x - obstacleX; dy = y - obstacleY
                d2 = dx*dx + dy*dy
                if d2 < minDist2 and d2 > 0.0:  # adicionei d2 > 0.0 para evitar divisão por zero
                    d = math.sqrt(d2)
                    # empurra partícula para fora do obstáculo
                    pushFactor = (minDist - d) / d
                    self.particlePos[2*i] += dx * pushFactor
                    self.particlePos[2*i+1] += dy * pushFactor
                    # reflexão da velocidade (quique)
                    nx = dx / d; ny = dy / d
                    vn = self.particleVel[2*i] * nx + self.particleVel[2*i+1] * ny
                    self.particleVel[2*i] -= 1.8 * vn * nx  # 1.8 para quique mais forte
                    self.particleVel[2*i+1] -= 1.8 * vn * ny
                
                # COLISÃO COM PAREDES DO TANQUE
                if x < minX:
                    x = minX; self.particleVel[2*i] = 0.0
                if x > maxX:
                    x = maxX; self.particleVel[2*i] = 0.0
                if y < minY:
                    y = minY; self.particleVel[2*i+1] = 0.0
                if y > maxY:
                    y = maxY; self.particleVel[2*i+1] = 0.0
                self.particlePos[2*i] = x; self.particlePos[2*i+1] = y
                
    def update_particle_density(self):
        if NUMBA_AVAILABLE:
            update_particle_density_nb(self.particlePos, self.particleDensity, self.numParticles,
                                       self.fNumX, self.fNumY, self.h, self.fInvSpacing)
        else:
            n = self.fNumY
            h2 = 0.5 * self.h
            self.particleDensity.fill(0.0)
            for i in range(self.numParticles):
                x = self.particlePos[2*i]; y = self.particlePos[2*i+1]
                x = clamp(x, self.h, (self.fNumX - 1) * self.h)
                y = clamp(y, self.h, (self.fNumY - 1) * self.h)
                x0 = int((x - h2) * self.fInvSpacing)
                tx = ((x - h2) - x0 * self.h) * self.fInvSpacing
                x1 = min(x0 + 1, self.fNumX-2)
                y0 = int((y - h2) * self.fInvSpacing)
                ty = ((y - h2) - y0 * self.h) * self.fInvSpacing
                y1 = min(y0 + 1, self.fNumY-2)
                sx = 1.0 - tx; sy = 1.0 - ty
                if 0 <= x0 < self.fNumX and 0 <= y0 < self.fNumY: self.particleDensity[x0*n + y0] += sx*sy
                if 0 <= x1 < self.fNumX and 0 <= y0 < self.fNumY: self.particleDensity[x1*n + y0] += tx*sy
                if 0 <= x1 < self.fNumX and 0 <= y1 < self.fNumY: self.particleDensity[x1*n + y1] += tx*ty
                if 0 <= x0 < self.fNumX and 0 <= y1 < self.fNumY: self.particleDensity[x0*n + y1] += sx*ty
        # set rest density once
        if self.particleRestDensity == 0.0:
            s = 0.0; num = 0
            for i in range(self.fNumCells):
                if self.cellType[i] == 0:
                    s += self.particleDensity[i]; num += 1
            if num > 0:
                self.particleRestDensity = s / num

    def transfer_velocities_to_grid(self):
        # prepare arrays (zero antes de chamar o JIT)
        self.prevU[:] = self.u
        self.prevV[:] = self.v

        self.du.fill(0.0)
        self.dv.fill(0.0)
        self.u.fill(0.0)
        self.v.fill(0.0)

        # inicializa cellType (2=SOLID,1=AIR)
        n = self.fNumY
        for i in range(self.fNumCells):
            self.cellType[i] = 2 if self.s[i] == 0.0 else 1

        # marque células que contém partículas como FLUID (0)
        for i in range(self.numParticles):
            x = self.particlePos[2*i]
            y = self.particlePos[2*i+1]
            xi = int(x * self.fInvSpacing); yi = int(y * self.fInvSpacing)
            if xi < 0: xi = 0
            if xi > self.fNumX - 1: xi = self.fNumX - 1
            if yi < 0: yi = 0
            if yi > self.fNumY - 1: yi = self.fNumY - 1
            cellNr = xi * n + yi
            if self.cellType[cellNr] == 1:
                self.cellType[cellNr] = 0  # FLUID_CELL == 0

        # chame o JIT para acumular u e v separadamente (caller já zerou arrays)
        if NUMBA_AVAILABLE:
            transfer_velocities_to_grid_nb(self.particlePos, self.particleVel, self.u, self.v, self.du, self.dv,
                                        self.numParticles, self.fNumX, self.fNumY, self.h, self.fInvSpacing, True)
            transfer_velocities_to_grid_nb(self.particlePos, self.particleVel, self.u, self.v, self.du, self.dv,
                                        self.numParticles, self.fNumX, self.fNumY, self.h, self.fInvSpacing, False)
        else:
            # fallback em Python (igual ao anterior, mas usando a mesma lógica)
            n = self.fNumY
            h2 = 0.5 * self.h
            # acumular componente u
            for i in range(self.numParticles):
                x = self.particlePos[2*i]; y = self.particlePos[2*i+1]
                x = clamp(x, self.h, (self.fNumX - 1) * self.h)
                y = clamp(y, self.h, (self.fNumY - 1) * self.h)
                x0 = int((x - 0.0) * self.fInvSpacing)
                if x0 > self.fNumX-2: x0 = self.fNumX-2
                if x0 < 0: x0 = 0
                tx = ((x - 0.0) - x0*self.h) * self.fInvSpacing
                x1 = x0 + 1
                y0 = int((y - h2) * self.fInvSpacing)
                if y0 > self.fNumY-2: y0 = self.fNumY-2
                if y0 < 0: y0 = 0
                ty = ((y - h2) - y0*self.h) * self.fInvSpacing
                y1 = y0 + 1
                sx = 1.0 - tx; sy = 1.0 - ty
                d0 = sx*sy; d1 = tx*sy; d2 = tx*ty; d3 = sx*ty
                nr0 = x0*n + y0; nr1 = x1*n + y0; nr2 = x1*n + y1; nr3 = x0*n + y1
                pv = self.particleVel[2*i]
                self.u[nr0] += pv * d0; self.du[nr0] += d0
                self.u[nr1] += pv * d1; self.du[nr1] += d1
                self.u[nr2] += pv * d2; self.du[nr2] += d2
                self.u[nr3] += pv * d3; self.du[nr3] += d3
            # acumular componente v
            for i in range(self.numParticles):
                x = self.particlePos[2*i]; y = self.particlePos[2*i+1]
                x = clamp(x, self.h, (self.fNumX - 1) * self.h)
                y = clamp(y, self.h, (self.fNumY - 1) * self.h)
                x0 = int((x - h2) * self.fInvSpacing)
                if x0 > self.fNumX-2: x0 = self.fNumX-2
                if x0 < 0: x0 = 0
                tx = ((x - h2) - x0*self.h) * self.fInvSpacing
                x1 = x0 + 1
                y0 = int((y - 0.0) * self.fInvSpacing)
                if y0 > self.fNumY-2: y0 = self.fNumY-2
                if y0 < 0: y0 = 0
                ty = ((y - 0.0) - y0*self.h) * self.fInvSpacing
                y1 = y0 + 1
                sx = 1.0 - tx; sy = 1.0 - ty
                d0 = sx*sy; d1 = tx*sy; d2 = tx*ty; d3 = sx*ty
                nr0 = x0*n + y0; nr1 = x1*n + y0; nr2 = x1*n + y1; nr3 = x0*n + y1
                pv = self.particleVel[2*i+1]
                self.v[nr0] += pv * d0; self.dv[nr0] += d0
                self.v[nr1] += pv * d1; self.dv[nr1] += d1
                self.v[nr2] += pv * d2; self.dv[nr2] += d2
                self.v[nr3] += pv * d3; self.dv[nr3] += d3

        # normaliza (divide por soma dos pesos)
        for k in range(len(self.u)):
            if self.du[k] > 0.0:
                self.u[k] /= self.du[k]
        for k in range(len(self.v)):
            if self.dv[k] > 0.0:
                self.v[k] /= self.dv[k]

        # restore solid cells (mesma lógica do original)
        for i in range(self.fNumX):
            for j in range(self.fNumY):
                idx = i * n + j
                solid = (self.cellType[idx] == 2)
                if solid or (i > 0 and self.cellType[(i-1)*n + j] == 2):
                    self.u[idx] = self.prevU[idx]
                if solid or (j > 0 and self.cellType[i*n + j - 1] == 2):
                    self.v[idx] = self.prevV[idx]


    def transfer_velocities_from_grid(self, flipRatio):
            # blended PIC/FLIP mapping com proteção contra valores extremos
            n = self.fNumY
            h2 = 0.5 * self.h
            for i in range(self.numParticles):
                x = self.particlePos[2*i]; y = self.particlePos[2*i+1]
                x = clamp(x, self.h, (self.fNumX - 1) * self.h)
                y = clamp(y, self.h, (self.fNumY - 1) * self.h)
                # ambos componentes
                for component in (0, 1):
                    dx = 0.0 if component == 0 else h2
                    dy = h2 if component == 0 else 0.0
                    x0 = int((x - dx) * self.fInvSpacing); x0 = max(0, min(self.fNumX-2, x0))
                    tx = ((x - dx) - x0*self.h) * self.fInvSpacing; x1 = x0 + 1
                    y0 = int((y - dy) * self.fInvSpacing); y0 = max(0, min(self.fNumY-2, y0))
                    ty = ((y - dy) - y0*self.h) * self.fInvSpacing; y1 = y0 + 1
                    sx = 1.0 - tx; sy = 1.0 - ty
                    d0 = sx*sy; d1 = tx*sy; d2 = tx*ty; d3 = sx*ty
                    nr0 = x0*n + y0; nr1 = x1*n + y0; nr2 = x1*n + y1; nr3 = x0*n + y1

                    # validade simples: trate AIR como inválido (1 == AIR)
                    valid0 = 1.0 if (self.cellType[nr0] != 1) else 0.0
                    valid1 = 1.0 if (self.cellType[nr1] != 1) else 0.0
                    valid2 = 1.0 if (self.cellType[nr2] != 1) else 0.0
                    valid3 = 1.0 if (self.cellType[nr3] != 1) else 0.0

                    d = valid0*d0 + valid1*d1 + valid2*d2 + valid3*d3
                    if d > 0.0:
                        fvals = 0.0
                        prevvals = 0.0
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

            # proteção final: clamp magnitude (evita picos numéricos)
            maxVel = 6.0
            vx = self.particleVel[0::2]; vy = self.particleVel[1::2]
            speed = np.sqrt(vx*vx + vy*vy)
            if np.any(speed > maxVel):
                mask = speed > maxVel
                scale = maxVel / (speed[mask] + 1e-12)
                vx[mask] *= scale; vy[mask] *= scale


    def solve_incompressibility(self, numIters, dt, overRelaxation, compensateDrift=True):
        if NUMBA_AVAILABLE:
            solve_incompressibility_nb(self.u, self.v, self.p, self.cellType, self.particleDensity,
                                       self.fNumX, self.fNumY, self.density, self.h, dt,
                                       numIters, overRelaxation, compensateDrift, self.particleRestDensity)
        else:
            # fallback python solver (calls similar algorithm)
            solve_incompressibility_nb(self.u, self.v, self.p, self.cellType, self.particleDensity,
                                       self.fNumX, self.fNumY, self.density, self.h, dt,
                                       numIters, overRelaxation, compensateDrift, self.particleRestDensity)

    def update_particle_colors(self):
        # simple color aging + density-based highlight
        for i in range(self.numParticles):
            s = 0.01
            self.particleColor[3*i] = clamp(self.particleColor[3*i] - s, 0.0, 1.0)
            self.particleColor[3*i+1] = clamp(self.particleColor[3*i+1] - s, 0.0, 1.0)
            self.particleColor[3*i+2] = clamp(self.particleColor[3*i+2] + s, 0.0, 1.0)
            # density-based
            h1 = self.fInvSpacing
            x = self.particlePos[2*i]; y = self.particlePos[2*i+1]
            xi = int(x * h1); yi = int(y * h1)
            xi = max(1, min(self.fNumX-1, xi)); yi = max(1, min(self.fNumY-1, yi))
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
        # simple coloring for grid (not used heavily in Python mode)
        self.cellColor.fill(0.0)
        for i in range(self.fNumCells):
            if self.cellType[i] == 2:
                self.cellColor[3*i] = 0.5; self.cellColor[3*i+1] = 0.5; self.cellColor[3*i+2] = 0.5
            elif self.cellType[i] == 0:
                d = self.particleDensity[i]
                if self.particleRestDensity > 0.0:
                    d /= self.particleRestDensity
                # map d in [0,2] to color (simple blue->yellow)
                val = min(max((d - 0.0)/(2.0-0.0), 0.0), 0.9999)
                r = min(max(2.0*(val-0.5), 0.0),1.0)
                g = min(max(1.0 - abs(val-0.5)*2.0, 0.0),1.0)
                b = min(max(1.0 - val*2.0, 0.0),1.0)
                self.cellColor[3*i]=r; self.cellColor[3*i+1]=g; self.cellColor[3*i+2]=b

    def simulate(self, dt, gravity, flipRatio, numPressureIters, numParticleIters,
                 overRelaxation, compensateDrift, separateParticles,
                 obstacleX, obstacleY, obstacleRadius):
        # single substep (no substepping for speed; you can add if needed)
        self.integrate_particles(dt, gravity)
        if separateParticles and self.numParticles > 0:
            self.push_particles_apart(numParticleIters)
        self.handle_particle_collisions(obstacleX, obstacleY, obstacleRadius, 0.0, 0.0)
        self.transfer_velocities_to_grid()
        self.update_particle_density()
        self.solve_incompressibility(numPressureIters, dt, overRelaxation, compensateDrift)
        self.transfer_velocities_from_grid(flipRatio)
        self.update_particle_colors()
        self.update_cell_colors()


# ----- simple pygame renderer & scene setup -----

def world_to_screen(x, y, sim_w, sim_h, screen_w, screen_h):
    sx = int(x / sim_w * screen_w)
    sy = int(screen_h - (y / sim_h * screen_h))
    return sx, sy

def run_simulation(resolution=100, particle_scale=0.3, show_particles=True, max_fps=30):
    pygame.init()
    screen_w, screen_h = 800, 600  # resolução menor
    screen = pygame.display.set_mode((screen_w, screen_h))
    clock = pygame.time.Clock()

    simHeight = 3.0
    cScale = screen_h / simHeight
    simWidth = screen_w / cScale

    tankHeight = 1.0 * simHeight
    tankWidth = 1.0 * simWidth
    h = tankHeight / resolution
    density = 1000.0

    relWaterHeight = 0.5  # reduzido
    relWaterWidth = 0.4   # reduzido
    r = particle_scale * h
    dx = 3.0 * r
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

    # setup grid cells for tank (s: 1 fluid, 0 solid)
    n = fluid.fNumY
    for i in range(fluid.fNumX):
        for j in range(fluid.fNumY):
            s = 1.0
            if i == 0 or i == fluid.fNumX-1 or j == 0:
                s = 0.0
            fluid.s[i*n + j] = s

    # scene params
    gravity = -9.81
    dt = 1.0 / 60.0
    flipRatio = 0.2
    numPressureIters = 30
    numParticleIters = 1
    overRelaxation = 1.5
    compensateDrift = True
    separateParticles = True
    obstacleX = 3.0
    obstacleY = 2.0
    obstacleRadius = 0.15

    running = True
    paused = False
    dragging = True  # NOVO: controle de arrasto
    last = time.time()
    frame = 0

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
            
            # NOVO: Sistema de arrasto
            elif ev.type == pygame.MOUSEBUTTONDOWN:
                if ev.button == 1:  # botão esquerdo
                    mx, my = pygame.mouse.get_pos()
                    ox = mx / screen_w * simWidth
                    oy = (screen_h - my) / screen_h * simHeight
                    
                    # verifica se clicou próximo ao obstáculo
                    dx = ox - obstacleX
                    dy = oy - obstacleY
                    dist = math.sqrt(dx*dx + dy*dy)
                    
                    if dist <= obstacleRadius * 1.5:  # área de clique um pouco maior
                        dragging = True
                    else:
                        # clique fora = teleporta obstáculo (comportamento anterior)
                        obstacleX, obstacleY = ox, oy
            
            elif ev.type == pygame.MOUSEBUTTONUP:
                if ev.button == 1:
                    dragging = False
            
            elif ev.type == pygame.MOUSEMOTION:
                if dragging:
                    mx, my = pygame.mouse.get_pos()
                    obstacleX = mx / screen_w * simWidth
                    obstacleY = (screen_h - my) / screen_h * simHeight
                    
                    # limita obstáculo dentro do tanque
                    obstacleX = clamp(obstacleX, obstacleRadius, tankWidth - obstacleRadius)
                    obstacleY = clamp(obstacleY, obstacleRadius, tankHeight - obstacleRadius)

        if not paused:
            fluid.simulate(dt, gravity, flipRatio, numPressureIters, numParticleIters,
                           overRelaxation, compensateDrift, separateParticles,
                           obstacleX, obstacleY, obstacleRadius)
        
        # draw
        screen.fill((0,0,0))
        if show_particles:
            for i in range(fluid.numParticles):
                x = fluid.particlePos[2*i]; y = fluid.particlePos[2*i+1]
                sx, sy = world_to_screen(x, y, simWidth, simHeight, screen_w, screen_h)
                c = (int(255*fluid.particleColor[3*i]), int(255*fluid.particleColor[3*i+1]), int(255*fluid.particleColor[3*i+2]))
                pygame.draw.circle(screen, c, (sx, sy), max(1, int(1 + fluid.particleRadius / simWidth * screen_w)))
        
        # obstacle - muda cor se está sendo arrastado
        ox_s, oy_s = world_to_screen(obstacleX, obstacleY, simWidth, simHeight, screen_w, screen_h)
        color = (255, 255, 0) if dragging else (255, 0, 0)  # amarelo quando arrastando
        pygame.draw.circle(screen, color, (ox_s, oy_s), max(2, int((obstacleRadius)/simWidth*screen_w)), 2)
        
        # info
        drag_text = " [ARRASTANDO]" if dragging else ""
        info = f"particles: {fluid.numParticles}  numba: {NUMBA_AVAILABLE}  fps: {int(clock.get_fps())}{drag_text}"
        txt = font.render(info, True, (255,255,255))
        screen.blit(txt, (10,10))
        
        pygame.display.flip()
        clock.tick(max_fps)
        frame += 1

    pygame.quit()


if __name__ == "__main__":
    # prints iniciais para debug (úteis ao rodar pelo terminal)
    print("Iniciando flip_fast.py")
    print(f"NumPy versão: {np.__version__}")
    print(f"Numba disponível: {NUMBA_AVAILABLE}")

    # Rode a simulação dentro de try/except para capturar erros
    try:
        # parametros que você já usa; ajuste se quiser
        run_simulation(resolution=90, particle_scale=0.25, show_particles=True, max_fps=60)

    except Exception as e:
        # registra o traceback em arquivo e no console
        tb = traceback.format_exc()
        print("Erro durante a execução. Traceback abaixo:\n")
        print(tb)
        with open("error.log", "w", encoding="utf-8") as f:
            f.write("Erro durante a execução do flip_fast.py\n")
            f.write(tb)

        # tenta abrir uma janela Pygame simples pra mostrar o erro e aguardar fechar manualmente
        try:
            pygame.init()
            screen = pygame.display.set_mode((800, 200))
            pygame.display.set_caption("Erro - veja error.log")
            font = pygame.font.SysFont("Arial", 14)
            lines = tb.splitlines()
            # limita linhas pra caber na tela
            lines = lines[-8:]
            running = True
            while running:
                for ev in pygame.event.get():
                    if ev.type == pygame.QUIT:
                        running = False
                    elif ev.type == pygame.KEYDOWN and ev.key == pygame.K_ESCAPE:
                        running = False
                screen.fill((20, 20, 20))
                y = 10
                screen.blit(font.render("Erro durante a execução. Veja error.log (ou console).", True, (255, 100, 100)), (10, y)); y += 22
                for ln in lines:
                    screen.blit(font.render(ln[:120], True, (220, 220, 220)), (10, y))
                    y += 18
                pygame.display.flip()
                pygame.time.wait(100)
            pygame.quit()
        except Exception:
            # se nem isso funcionar, apenas finalize
            print("Falha ao exibir janela de erro.")

    else:
        # execução terminou sem lançar exceção (programa fechou normalmente)
        print("Execução finalizada sem exceções. Caso você não esperasse o encerramento, verifique os parâmetros (dt, resolução, numPressureIters).")