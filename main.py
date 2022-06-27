import numpy as np
from ann_criterion import optimality_criterion

def func(x):
    return 0.01*((x[0]-1)**2+2 *(x[1]-1)**2)* ((x[0]+1)**2+2 *(x[1]+1)**2+0.5)* ((x[0]+2)**2+2* (x[1]-2)**2+0.7)
def PSO(func, dim, N):  # parametri: (funkcija, dimenzionalnost, broj iteracija)
    p_best_val_list = []
    n_particles = 600  # broj cestica
    X = np.random.rand(dim, n_particles) * 5  # koordinate cestica (nasumicno izabrane)
    V = np.random.randn(dim, n_particles) * 0.1  # brzine cestica (nasumicno izabrane)
    p_best = X  # najbolje koordinate svake cestice
    for i in range(n_particles):
        p_best_val_list.append(
            func(np.array(X[:, i])))  # izracunavanje vrednosti funkcije za date najbolje koordinate cestica

    p_best_val = np.array(p_best_val_list)
    g_best = p_best[:, p_best_val.argmin()]  # globalno najbolja koordinata (sa najmanjom vrednoscu funkcije)
    g_best_val = p_best_val.min()  # globalno najbolja (najmanja) vrednost funkcije
    for i in range(N):
        w = 0.9 - i / N * (0.9 - 0.4)  # konstanta inertnosti
        cp = 2.5 - i / N * (2.5 - 0.5)  # kognitivna konstanta
        cg = 0.5 + i / N * (2.5 - 0.5)  # socijalna konstanta
        rp, rg = np.random.rand(2)  # random vrednosti izmedju 0 i 1
        V = w * V + cp * rp * (p_best - X) + cg * rg * (g_best.reshape(-1, 1) - X)  # izracunavanje novih brzina
        X = X + V  # izracunavanje novih koordinata
        new_val_list = []
        for i in range(n_particles):
            new_val_list.append(func(np.array(X[:, i])))  # izracunavanje vrednosti funkcija za nove koordinate
        new_val = np.array(new_val_list)
        p_best[:, (p_best_val >= new_val)] = X[:, (p_best_val >= new_val)]  # maska koja postavlja nove licne najbolje koordinate ukoliko je nova
                                                                            # vrednost funkcije bolja(manja) od prethodne
        p_best_val = np.array([p_best_val, new_val]).min(axis=0)  # odabir bolje(manje) vrednosti funkcije svake cestice
        g_best = p_best[:, p_best_val.argmin()]  # azuriranje globalno najboljih koordinata na osnovu vrednosti funkcije
        g_best_val = p_best_val.min()  # azuriranje globalno najbolje vrednosti funkcije
    return g_best, g_best_val


if __name__ == '__main__':
    g_best,g_best_val=PSO(lambda x:func(x),2,50)
    print(g_best)
    print(g_best_val)
