import numpy as np
import time
import os
import shutil

def clear_terminal():
    print('\033[2J\033[H', end='')

#powinno byc rozwiaz macierz, ale w Vimie nie podkresla mi gdzie sie pojawia ta funkcja wiec boje sie zmieniac
def stworz_macierz(t, A_inv, var_to_coords, coords_to_var, delta_x, delta_z, N, L, h_val, T):
    """
    Parametry:
      - t: czas,
      - A_inv:macierz odwrotna do macierzy uzywanej w stworz macierz,
      - var_to_coords, coords_to_var: mapowanie danych,
      - delta_x, delta_z: kroki siatki w kierunkach x i z,
      - N: liczba punktów podzialu każdym kierunku,
      - L: długość fali,
      - h_val: głębokość cieczy,
      - T: okres,
    """
    g = 9.81
    A_coeff = (g * h_val * T) / (4 * np.pi)

    def licz_B(z, h, L):
        arg = 2 * np.pi / L
        return (np.exp(arg * (z + h)) + np.exp(-arg * (z + h))) / (np.exp(arg * h) + np.exp(-arg * h))

    def licz_C(x, t, L, T):
        return np.sin(2 * np.pi * (x / L + t / T))

    def phi_analityczne(x, z, t):
        return A_coeff * licz_B(z, h_val, L) * licz_C(x, t, L, T)

    num_vars = A_inv.shape[0]

    #wektor rozwiazan (metoda reszt skonczonych)
    b = np.zeros(num_vars)
    dx2 = delta_x ** 2
    dz2 = delta_z ** 2

    for k in range(num_vars):
        i, j = var_to_coords[k]
        neighbors = [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]
        for ni, nj in neighbors:
            x = ni * delta_x
            z = h_val + nj * delta_z
            if not (1 <= ni <= N - 2 and 1 <= nj <= N - 2):
                phi_brzeg = phi_analityczne(x, z, t)
                # Określ kierunek sąsiada i współczynnik
                if ni != i:  # Sąsiad w kierunku x
                    coeff = 1 / dx2
                else:         # Sąsiad w kierunku z
                    coeff = 1 / dz2
                b[k] -= coeff * phi_brzeg
    
    return np.linalg.solve(A_inv, b) #rozwiazywanie macierzy AKA obliczanie nowego wektora zawierajacego aktualna wartosc phi dla kazdego punktu


def porownaj_rozwiazania(t, A_inv, var_to_coords, coords_to_var, delta_x, delta_z, N, L, h_val, T):
          
        phi_num = stworz_macierz(t, A_inv, var_to_coords, coords_to_var, delta_x, delta_z, N, L, h_val, T)

        g = 9.81
        A_coeff = (g * h_val * T) / (4 * np.pi)

        def licz_B(z, h, L):
            arg = 2 * np.pi / L
            return (np.exp(arg * (z + h)) + np.exp(-arg * (z + h))) / (np.exp(arg * h) + np.exp(-arg * h))

        def licz_C(x, t, L, T):
            return np.sin(2 * np.pi * (x / L + t / T))

        def phi_analityczne(x, z, t):
            return A_coeff * licz_B(z, h_val, L) * licz_C(x, t, L, T)

        num_vars = len(phi_num)
        phi_ana = np.zeros(num_vars)
        for k in range(num_vars):
            i, j = var_to_coords[k]
            x = i * delta_x
            z = h_val + j * delta_z
            phi_ana[k] = phi_analityczne(x, z, t)

        error = np.abs(phi_num - phi_ana)
        
        # Obliczanie błędu procentowego z zabezpieczeniem przed dzieleniem przez zero
        epsilon = 1e-12
        error_percent = 100 * np.abs((phi_num - phi_ana) / (np.abs(phi_ana) + epsilon))
        
        max_error = np.max(error)
        max_percent = np.max(error_percent)
        mean_error = np.mean(error)
        mean_percent = np.mean(error_percent)

        print("Dla N =", N)
        print(f"Maksymalny błąd bezwzględny: {max_error:.4e} | procentowy: {max_percent:.2f}%")
        print(f"Średni błąd bezwzględny: {mean_error:.4e} | procentowy: {mean_percent:.2f}%\n")
        
        return phi_num, phi_ana, error

def main():
    # Parametry symulacji

    print("Porównanie rozwiązań dla t =", 0.0)
    for N in range(15, 40):# Liczba punktów w obu kierunkach (siatka NxN)
        L = 40.0        # Długość domeny w kierunku x
        h_val = 10.0    # Głębokość cieczy
        T = 5.0         # okres fali
    
        #na ile krokow dzielimy siatke
        Nx = Nz = N
        delta_x = L / (Nx - 1)
        delta_z = h_val / (Nz - 1)
    
    # Liczba punktów wnętrza – patrz Metoda Skonczonych roznic
        num_vars = (Nx - 2) * (Nz - 2)

    # Budujemy macierz A dla roznic skonczonych
        A = np.zeros((num_vars, num_vars)) #AKA dyskretyzacja operatora Laplace'a
        var_to_coords = {}
        coords_to_var = {}
        idx = 0
        for j in range(1, Nz - 1):
            for i in range(1, Nx - 1):
                var_to_coords[idx] = (i, j)
                coords_to_var[(i, j)] = idx
                idx += 1


    # Wypelnianie tablicy zaleznosciami punktow miedzy soba
        dx2 = delta_x ** 2
        dz2 = delta_z ** 2

# Wypełnianie macierzy A – dyskretyzacja operatora Laplace’a z różnym dx i dz
        for k in range(num_vars):
            i, j = var_to_coords[k]
            A[k, k] = -2 * (1/dx2 + 1/dz2)  # element centralny

            neighbors = [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]
            for ni, nj in neighbors:
                if 1 <= ni <= Nx - 2 and 1 <= nj <= Nz - 2:
                    neighbor_idx = coords_to_var[(ni, nj)]
                    if ni != i:
                        A[k, neighbor_idx] = 1/dx2  # Poprawiony współczynnik dla X
                    elif nj != j:
                        A[k, neighbor_idx] = 1/dz2  # Poprawiony współczynnik dla Z# sąsiad w pionie   
    #macierz odwrotna dodana inicjowana tutaj w celu optymalizacji, gdyz stosunki miedzy punktami sie nie zmieniaja w zaleznosci od czasu
        A_inv = np.linalg.inv(A) #odwrotna dlatego, ze A^(-1) * b = phi


    #Inicjalizacja fali t - czas startu, phi_interior > positions - pozycja startowa czastek fali,
        t = 0.0
        phi_interior = stworz_macierz(t, A_inv, var_to_coords, coords_to_var, delta_x, delta_z, N, L, h_val, T)
        phi_num, phi_ana, error = porownaj_rozwiazania(t, A, var_to_coords, coords_to_var, delta_x, delta_z, N, L, h_val, T)
 
if __name__ == "__main__":
    main()
