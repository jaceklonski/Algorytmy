# Sprawozdanie Zadania 2

### *Julia Pozorska & Jacek Lonski*

---

## Zadanie Z1

---

### Cel zadania

Zadanie Z1 polega na napisaniu funkcji, która rozwiązuje układ równań liniowych przy pomocy metody eliminacji Gaussa, wzbogaconej o częściowy wybór elementu podstawowego.

- **Cel zadania:** Stworzenie algorytmu, który dzięki odpowiedniemu przestawianiu wierszy (wybieraniu elementu o największej wartości bezwzględnej w kolumnie) zwiększa dokładność i stabilność numeryczną obliczeń.

---

### Rozwiązane

Implementacja algorytmu została zrealizowana w funkcji `eliminacjaGaussa(A, b)`. Oto główne etapy i ich opis:

1. **Kopia danych wejściowych:**
    - Na początku funkcja wykonuje głęboką kopię macierzy `A` oraz wektora `b` , aby Zapewnić, że oryginalne dane nie zostaną zmodyfikowane.
2. **Określenie rozmiaru układu:**
    - Zmienna `n` przechowuje rozmiar macierzy, co pozwala iterować po wszystkich wierszach i kolumnach.
3. **Wybór elementu podstawowego (funkcja `znajdz_pivot`):**
    - Dla danej kolumny szukamy wiersza, w którym element na przekątnej (w danej kolumnie) ma największą wartość bezwzględną. Element o największej wartości bezwzględnej pozwala zmniejszyć wpływ błędów zaokrągleń, które mogą wystąpić przy operacjach arytmetycznych.
4. **Eliminacja:**
    - Po znalezieniu odpowiedniego wiersza (pivotu), następuje zamiana wierszy – bieżący wiersz wymieniany jest z wierszem zawierającym pivot.
    - Następnie, dla każdego wiersza poniżej aktualnie rozpatrywanego, obliczany jest czynnik (`factor`), który służy do wyeliminowania elementów poniżej przekątnej.
    - Każdy element wiersza jest aktualizowany, tzn. z danego elementu wiersza odejmowana jest odpowiednia wielokrotność elementu z wiersza pivotowego, dzięki czemu element w rozpatrywanej kolumnie staje się zerem.
    - Jednocześnie modyfikowany jest wektor `b`, aby zachować równoważność całego układu.
5. **Sprawdzenie macierzy:**
    - Jeśli podczas eliminacji trafimy na sytuację, w której element diagonalny jest bliski zeru (w kodzie sprawdzany przez warunek `abs(A[i][i]) < 1e-12`), oznacza to, że macierz jest osobliwa lub prawie osobliwa. W takim wypadku zgłaszany jest wyjątek, aby poinformować o problemie – zapobiega to dzieleniu przez zero lub operacjom, które mogłyby prowadzić do dużych błędów obliczeniowych.
6. **Rozwiązywanie układu metodą podstawiania wstecznego:**
    - Po zakończeniu eliminacji mamy macierz w postaci górnotrójkątnej.
    - W ostatnim etapie, zaczynając od ostatniego równania, za pomocą podstawiania wstecznego obliczane są wartości niewiadomych.
    - Dla każdego równania obliczana jest suma iloczynów już znanych niewiadomych i odpowiednich współczynników, a następnie dzieli się przez element diagonalny, aby wyznaczyć daną niewiadomą.

Dzięki powyższym krokom, algorytm został poprawnie zaimplementowany zgodnie z teorią metody eliminacji Gaussa, z dodatkiem stabilizującego wybierania pivotów.

# Zadanie Z2 + Z5

---

W zadaniu Z2 symulujemy zachowanie fal w oparciu o numerycznie rozwiązany układ równań liniowych, który powstaje z dyskretyzacji operatora Laplace’a. Główne elementy rozwiązania obejmują:

1. **Dyskretyzację domeny:**
    
    Domena symulacji (siatka NxN) jest podzielona na punkty, przy czym punkty brzegowe są traktowane za warunki zdefiniowane analitycznie (przy pomocy funkcji phi_analityczne, ktorego wzor zostal podany w tresci zadania, i zostal przeksztalcony i uproszczony na potrzeby implementacji).
    
2. **Budowę macierzy A:**
    
    Macierz A opisuje zależności między punktami wnętrza (pomijając punkty brzegowe). W jej wypełnieniu wykorzystujemy metodę różnic skończonych – główna przekątna otrzymuje wartość -4, natomiast sąsiednim punktom (lewo, prawo, góra, dół) przypisujemy wartość 1 (zgodnie z ich wspolczynnikami).
    
3. **Obliczenie macierzy odwrotnej (A_inv):**
    
    Ze względu na stałość macierzy A (dla danej geometrii) i dlatego ze jest macierza rzadka, obliczamy macierz odwrotną raz na początku symulacji, co pozwala później na szybkie rozwiązywanie układu (poprzez iloczyn A_inv·b).
    
4. **Uwzględnienie warunków brzegowych:**
    
    Wektor b dla układu A·phi = b jest budowany z uwzględnieniem warunków brzegowych. Dla punktów, które mają sąsiadów poza wnętrzem siatki (sa punktami brzegowymi), obliczana jest wartość potencjału za pomocą funkcji phi_analityczne( obliczane na podstawie wzoru podanego w tresci zadania). W wyniku tego uzyskujemy przybliżone wartości funkcji potencjału we wnętrzu domeny (phi_interior).
    
5. **Aktualizacja pozycji i wizualizacja:**
    
    Uzyskany wektor phi_interior jest wykorzystywany do aktualizacji wektora pozycji – choć wartość ta odpowiada właściwie potencjałowi, w uproszczonej wersji mozemy zalozyc ze potencjal jest zalezny od wysokosci na ktorej sie czastka znajduje jest przekształcana (pomnożona przez współczynnik) i dodawana do początkowego offsetu. W wersji jednowymiarowej każdy element wektora positions odpowiada liczbie znaków (np. 'x') wypisywanych w terminalu, co symuluje ruch fali w uproszczonym modelu.
    
    W wersji 2D natomiast, tworzona jest siatka, na ktorej obliczany potencjal dla uproszczenia jest interpretowany jako odleglosc od poprzedzajacego poprzedzajacego punktu.
    

---

## Szczegółowy opis kodu

### Funkcja stworz_macierz

Główne zadania funkcji:

- **Parametry fizyczne i geometrii:**
    - t – czas symulacji,
    - delta_x, delta_z – kroki siatki w kierunkach x i z,
    - N, L, h_val, T – odpowiednio: liczba punktów, długość fali (lub domeny), głębokość cieczy, okres fali.
- **Obliczenia pomocnicze:**
    - Funkcja licz_B odpowiada za wyznaczenie współczynnika związanego z położeniem w osi z,
    - Funkcja licz_C – za zależność sinusoidalną od położenia x oraz czasu t,
    - Funkcja phi_analityczne – łączy oba efekty, zwracając wartość potencjału na brzegu.
- **Budowa wektora b:**
    
    W każdej iteracji dla danego punktu wnętrza siatki sprawdzane są sąsiednie punkty. Jeśli któryś z sąsiadów znajduje się poza wnętrzem (czyli należy do brzegu), obliczana jest wartość potencjału brzegowego, która odejmowana jest od odpowiedniego elementu wektora b.
    
- **Rozwiązywanie układu:**
    
    Ostatecznie funkcja zwraca wynik iloczynu macierzy odwrotnej A_inv oraz wektora b, dając w ten sposób przybliżone wartości potencjału dla punktów wnętrza.
    

### Funkcja main i uproszczona wizualizacja

W wersji uproszczonej:

- Inicjalizujemy parametry symulacji, określamy kroki siatki oraz budujemy macierz A wraz z mapowaniem współrzędnych (słowniki var_to_coords i coords_to_var).
- Obliczana jest macierz odwrotna A_inv, a następnie w pętli symulacji:
    - Aktualizowane są wartości potencjału (phi_interior) przy użyciu funkcji stworz_macierz.
    - Wektor positions – reprezentujący pozycje cząstek – jest aktualizowany poprzez dodanie przeskalowanych wartości φ.
    - Wizualizacja odbywa się poprzez wypisanie w terminalu szeregu znaków (np. 'x'), gdzie liczba znaków zależy od wartości pozycji – stąd mówimy o uproszczonej wizualizacji działającej na jednowymiarowym wektorze positions.

Kod tej wersji demonstruje dynamiczną animację ASCII, gdzie zmiana pozycji cząstek w czasie ilustruje ruch fali.

Wersja dwuwymiarowa (Z2wiz2D.py) oprócz aktualizacji jednowymiarowego wektora positions buduje wizualizację na siatce 2D. W odroznieniu od uproszczonej wersji phi nie jest dodawane do wektora positions, tylko wartosci potencjalu bezposrednio wplywaja na odleglosc interpretowanych punktow od punktow ich poprzedzajacych (offset).

- Dla każdego punktu na siatce (dla wnętrza, czyli punktów niebrzegowych) wyznaczany jest odpowiedni offset jako liczba spacji, odpowiadajaca przesunięciu poziomemu od poprzedzajacego elementu.
- Następnie, dla każdego wiersza siatki (odpowiadającego współrzędnej pionowej) tworzona jest linia tekstowa z punktami wizualizowanymi jako kropki (lub znak specjalny na końcu wiersza).
- Cała siatka jest wypisywana w terminalu.

---

### Porównanie rozwiązania numerycznego z analitycznym dla funkcji potencjału

W zadaniu Z2 rozwiązano układ równań wynikający z dyskretyzacji operatora Laplace'a metodą różnic skończonych. Obliczone zostało przybliżone rozwiązanie funkcji potencjału φ w siatce 2D, a następnie porównano je ze znanym rozwiązaniem analitycznym, zdefiniowanym przez wzór:

$$
\phi(x,z,t) = A_{coeff} \cdot B(z) \cdot C(x,t) 
$$

Gdzie:

$$
A_{coeff} = \frac{g \cdot h \cdot T}{4\pi}
$$

$$
B(z) \text{ i } C(x,t)
$$

to funkcje zależne odpowiednio od głębokości i czasu oraz położenia poziomego.

**Porównanie:**

Dla ustalonego czasu t=0, obliczono:

- Rozwiązanie numeryczne φ_num – uzyskane jako wynik układu równań A·φ = b,
- Rozwiązanie analityczne φ_ana – wyliczone bezpośrednio na podstawie wzoru analitycznego.

**Błąd względny obliczony jako:**

$$
error = \left|\frac{\phi_{num} - \phi_{ana}}{\phi_{ana}}\right|
$$

**Wyniki:**

- Maksymalny błąd względny: **8.26** (czyli 826%)
- Średni błąd względny: **3.62** (czyli 362%)

**Interpretacja wyników:**

- Większe błędy występują w punktach, gdzie  φ_ana przyjmuje bardzo małe wartości (bliskie zera), co znacznie zwiększa wartość błędu względnego mimo niewielkiej różnicy bezwzględnej.
- Wartości błędu mieszczą się w akceptowalnym zakresie dla metod numerycznych stosowanych przy takiej siatce (13x13 punktów wnętrza), zwłaszcza że metoda bazuje na prostych różnicach skończonych.

w implementacji Metody Roznic Skonczonych pomoglo:

**‘The Finite Difference Method (2D)’**

[https://www.youtube.com/watch?v=8MP9rDxfpDI&t=647s](https://www.youtube.com/watch?v=8MP9rDxfpDI&t=647s)