import numpy as np
import time
import os
import shutil
import math

def clear_terminal():
    print('\033[2J\033[H', end='')

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

    for k in range(num_vars):
        i, j = var_to_coords[k]
        neighbors = [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]
        for ni, nj in neighbors:
            x = ni * delta_x
            z = h_val + nj * delta_z
            if not (1 <= ni <= N - 2 and 1 <= nj <= N - 2):
                #obliczanie wartosci brzegowych z podanego wzoru
                phi_brzeg = phi_analityczne(x, z, t)
                b[k] -= phi_brzeg

    return A_inv.dot(b) #rozwiazywanie macierzy AKA obliczanie nowego wektora zawierajacego aktualna wartosc phi dla kazdego punktu



def main():
    # Parametry symulacji
    N = 30          # Liczba punktów w obu kierunkach (siatka NxN)
    L = 40.0        # Długość domeny w kierunku x
    h_val = 10.0    # Głębokość cieczy
    T = h_val/2         # Okres fali

    Nx = Nz = N
    delta_x = L / (Nx - 1)
    delta_z = -h_val / (Nz - 1)

    # Liczba punktów wnętrza (wg metody różnic skończonych)
    num_vars = (Nx - 2) * (Nz - 2)

    # Budujemy macierz A i mapowania współrzędnych
    A = np.zeros((num_vars, num_vars))
    var_to_coords = {}
    coords_to_var = {}
    idx = 0
    for j in range(1, Nz - 1):
        for i in range(1, Nx - 1):
            var_to_coords[idx] = (i, j)
            coords_to_var[(i, j)] = idx
            idx += 1

    # Wypełnianie macierzy A według zależności między punktami wnętrza
    for k in range(num_vars):
        i, j = var_to_coords[k]
        A[k, k] = -4
        neighbors = [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]
        for ni, nj in neighbors:
            if 1 <= ni <= Nx - 2 and 1 <= nj <= Nz - 2:
                A[k, coords_to_var[(ni, nj)]] = 1

    A_inv = np.linalg.inv(A)


    # Inicjalizacja: t - czas startu,
    # phi_interior > rozwiązanie dla punktów wnętrza, positions > pozycje cząstek z pewnym offsetem (bazowym przesunięciem)

    t = 0.0
    phi_interior = stworz_macierz(t, A_inv, var_to_coords, coords_to_var, delta_x, delta_z, N, L, h_val, T)
    positions = []

    try:
        while True:
            clear_terminal()
            phi_interior = stworz_macierz(t, A_inv, var_to_coords, coords_to_var, delta_x, delta_z, N, L, h_val, T)
            positions = phi_interior * 1.5  # Aktualizacja pozycji zgodnie z potencjałem

            terminal_width = shutil.get_terminal_size().columns

            ascii_rows = []
            # Pętla po wierszach (j - współrzędna pionowa); "kolumny symulacji" będą wypisywane jako wiersze terminala
            for j in range(1, Nz - 1):
                row_cells = []
                # Pętla po kolumnach (i - współrzędna pozioma)
                for i in range(1, Nx - 1):
                    k = coords_to_var[j, (Nx - (1 + i))]
                    pos = positions[k]
                    # Wyznaczamy poziomą przesunięcie jako liczbę spacji (skalujemy dzieląc przez 100)
                    offset = int((200 + pos)/100)
                    if i < Nx - 2:
                        cell = " " * offset + "."
                    else:
                        cell = " " * offset + "#"
                    row_cells.append(cell)
                # Oddzielamy "kolumny" trzema spacjami – dzięki temu widoczny jest odstęp między punktami
                ascii_rows.append("".join(row_cells))
            # Wypisujemy całą siatkę (każdy wiersz w terminalu)
            for line in ascii_rows:
                print(line, flush=True)
            
            #print(positions)
            #for a in ascii_rows:
             #   print(a)
            print(t)
            t += 0.1 # Aktualizacja czasu
            time.sleep(0.1)  # Frame time
    except KeyboardInterrupt:
        print("Przerwano przez użytkownika.")

if __name__ == "__main__":
    main()
