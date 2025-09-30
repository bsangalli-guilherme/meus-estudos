import numpy as np  # Biblioteca para computação numérica, essencial para trabalhar com arrays (matrizes)
import matplotlib  # Biblioteca principal para criação de gráficos e visualizações

matplotlib.use("TkAgg")

# Importa o submódulo 'pyplot' para criar os gráficos de forma mais simples.
from matplotlib import pyplot as plt


# Define a frequência com que o gráfico será atualizado.
# Um novo quadro da animação será gerado a cada 20 iterações do loop principal.
plot_every = 20


# --- Função Principal ---
def main():
    
    # --- Parâmetros da Simulação ---
    Nx = 400  # Número de pontos na grade (lattice) na direção x
    Ny = 100  # Número de pontos na grade (lattice) na direção y
    tau = 0.53  # Tempo de relaxação. Relacionado à viscosidade do fluido.
    Nt = 3000  # Número total de passos de tempo (iterações) da simulação.

    # --- Configuração do Modelo Lattice Boltzmann (D2Q9) ---
    # D2Q9 significa "2 Dimensões, 9 Velocidades". É um modelo padrão para LBM.
    NL = 9  # Número de direções de velocidade discretas.
    
    # Índices das direções opostas. Usado para a condição de contorno "bounce-back" no obstáculo.
    # Ex: a direção 1 (Norte) é oposta à 5 (Sul). Isso torna o código mais legível.
    OPPOSITE_INDICES = [0, 5, 6, 7, 8, 1, 2, 3, 4]
    
    # Componentes x das 9 velocidades.
    cxs = np.array([0, 0, 1, 1, 1, 0, -1, -1, -1])
    # Componentes y das 9 velocidades.
    cys = np.array([0, 1, 1, 0, -1, -1, -1, 0, 1])
    # Pesos associados a cada direção.
    weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])

    # --- Condições Iniciais ---
    # F é a função de distribuição de partículas.
    # Começa com densidade 1, mais um pequeno ruído aleatório para quebrar a simetria.
    F = np.ones((Ny, Nx, NL)) + 0.1 * np.random.randn(Ny, Nx, NL)
    
    # Aumenta a densidade de partículas se movendo para a direita (direção 3)
    # para criar um fluxo inicial da esquerda para a direita.
    F[:, :, 3] = 2.9
    
    # --- Definição do Obstáculo (Forma Vetorizada) ---
    # Cria uma grade de coordenadas X e Y. Isso evita o uso de um loop for duplo.
    Y, X = np.ogrid[0:Ny, 0:Nx]
    # Calcula a distância de todos os pontos ao centro do cilindro de uma só vez.
    dist_sq = (X - Nx // 4)**2 + (Y - Ny // 2)**2
    # O cilindro é True onde a distância ao quadrado é menor que o raio ao quadrado (13*13).
    cylinder = dist_sq < 13**2

    # --- Configurações para a Animação do Gráfico ---
    plt.ion()
    fig, ax = plt.subplots()
    # Cria a imagem inicial. Usaremos 'im.set_data' para atualizá-la, o que é mais rápido.
    im = ax.imshow(np.zeros((Ny, Nx)), origin='lower', cmap='plasma')
    fig.colorbar(im, label='Velocidade') # Adiciona uma barra de cores.
    
    # --- Loop Principal da Simulação ---
    for it in range(Nt):
        if it % 10 == 0: # Imprime o progresso a cada 10 iterações
             print(it)

        # --- Passo de Streaming (Propagação) ---
        F[:, -1, [6, 7, 8]] = F[:, -2, [6, 7, 8]]
        F[:, 0, [2, 3, 4]] = F[:, 1, [2, 3, 4]]

        for i, cx, cy in zip(range(NL), cxs, cys):
            F[:, :, i] = np.roll(np.roll(F[:, :, i], cx, axis=1), cy, axis=0)

        # --- Condição de Contorno do Obstáculo (Bounce-back) ---
        bndryF = F[cylinder, :]
        bndryF = bndryF[:, OPPOSITE_INDICES] # Usa a constante para inverter as direções.
        F[cylinder, :] = bndryF

        # --- Cálculo das Variáveis Macroscópicas ---
        rho = np.sum(F, 2)
        ux = np.sum(F * cxs, 2) / rho
        uy = np.sum(F * cys, 2) / rho
        
        # Força a velocidade a ser zero dentro do obstáculo.
        ux[cylinder] = 0
        uy[cylinder] = 0

        # --- Passo de Colisão (Forma Vetorizada) ---
        # Esta é a principal otimização de velocidade. O cálculo é feito para todos os
        # pontos e direções de uma vez, sem um loop Python.
        
        # Termos da equação de equilíbrio
        usq = ux**2 + uy**2
        cu = 3 * (cxs[np.newaxis, np.newaxis, :] * ux[:, :, np.newaxis] + cys[np.newaxis, np.newaxis, :] * uy[:, :, np.newaxis])
        
        # Calcula Feq para todas as células de uma só vez usando broadcasting do NumPy.
        Feq = weights * rho[:, :, np.newaxis] * (1 + cu + 0.5 * cu**2 - 1.5 * usq[:, :, np.newaxis])
        
        # Relaxa a distribuição F em direção ao equilíbrio Feq.
        F += -(1 / tau) * (F - Feq)

        # --- Visualização ---
        if (it % plot_every == 0):
            speed = np.sqrt(usq)
            # Mascara o cilindro para que ele seja desenhado com uma cor distinta (preto, por padrão).
            speed_masked = np.ma.masked_array(speed, mask=cylinder)
            
            im.set_data(speed_masked) # Atualiza os dados da imagem (mais rápido que recriar).
            im.set_clim(vmin=0, vmax=np.max(speed)) # Atualiza os limites de cor.
            ax.set_title(f"Iteração = {it}")
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(0.001) # Pausa pequena para a GUI processar

    # --- Finalização ---
    plt.ioff()
    # Adiciona um título final e mostra o último quadro.
    ax.set_title(f"Simulação Finalizada em {Nt} iterações")
    plt.show()


# Ponto de entrada padrão para um script Python.
if __name__ == "__main__":
    main()

