# Sprawozdanie Zadania 2

Poniżej znajduje się opracowanie zadania Z1, w którym implementujemy metodę eliminacji Gaussa z częściowym wyborem elementu podstawowego. Poniższe omówienie wyjaśnia cel zadania, opisuje krok po kroku, jak została rozwiązana metoda, oraz tłumaczy to w przystępny sposób.

---

### Cel zadania

Zadanie Z1 polega na napisaniu funkcji, która rozwiązuje układ równań liniowych przy pomocy metody eliminacji Gaussa, wzbogaconej o **częściowy wybór elementu podstawowego**.

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

# Ogólne omówienie zadania Z2

W zadaniu Z2 symulujemy zachowanie fal w oparciu o numerycznie rozwiązany układ równań liniowych, który powstaje z dyskretyzacji operatora Laplace’a. Główne elementy rozwiązania obejmują:

1. **Dyskretyzację domeny:**
    
    Domena symulacji (siatka NxN) jest podzielona na punkty, przy czym punkty brzegowe są traktowane za warunki zdefiniowane analitycznie (przy pomocy funkcji phi_analityczne, ktorego wzor zostal podany w tresci zadania, i zostal przeksztalcony i uproszczony na potrzeby implementacji).
    
2. **Budowę macierzy A:**
    
    Macierz A opisuje zależności między punktami wnętrza (pomijając punkty brzegowe). W jej wypełnieniu wykorzystujemy metodę różnic skończonych – główna przekątna otrzymuje wartość -4, natomiast sąsiednim punktom (lewo, prawo, góra, dół) przypisujemy wartość 1.
    
3. **Obliczenie macierzy odwrotnej (A_inv):**
    
    Ze względu na stałość macierzy A (dla danej geometrii) i ze jest macierza rzadka, obliczamy macierz odwrotną raz na początku symulacji, co pozwala później na szybkie rozwiązywanie układu (poprzez iloczyn A_inv·b).
    
4. **Uwzględnienie warunków brzegowych:**
    
    Wektor b dla układu A·phi = b jest budowany z uwzględnieniem warunków brzegowych. Dla punktów, które mają sąsiadów poza wnętrzem siatki (sa punktami brzegowymi), obliczana jest wartość potencjału za pomocą funkcji phi_analityczne( obliczane na podstawie podanego wzoru ). W wyniku tego uzyskujemy przybliżone wartości funkcji potencjału we wnętrzu domeny (phi_interior).
    
5. **Aktualizacja pozycji i wizualizacja:**
    
    Uzyskany wektor phi_interior jest wykorzystywany do aktualizacji wektora pozycji – choć wartość ta odpowiada właściwie potencjałowi, w uproszczonej wersji mozemy zalozyc ze potencjal jest zalezny od wysokosci na ktorej sie czastka znajduje jest przekształcana (pomnożona przez współczynnik) i dodawana do początkowego offsetu. W wersji jednowymiarowej każdy element wektora positions odpowiada liczbie znaków (np. 'x') wypisywanych w terminalu, co symuluje ruch fali w uproszczonym modelu.
    

---

## Szczegółowy opis kodu (wersja 1D)

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

Kod tej wersji (widoczny w zadanym fragmencie) demonstruje dynamiczną animację ASCII, gdzie zmiana pozycji cząstek w czasie ilustruje ruch fali.

---

## Wersja dwuwymiarowa (plik Z2wiz2D.py)

W celu uzyskania bardziej rozbudowanej wizualizacji w przestrzeni dwuwymiarowej, zaimplementowano rozszerzoną wersję symulacji w pliku **Z2wiz2D.py**. Poniżej przedstawiam strukturę i kluczowe elementy kodu tej wersji:

### Kluczowe różnice i rozszerzenia w wersji 2D

- **Dwuwymiarowa wizualizacja:**
    
    Wersja dwuwymiarowa oprócz aktualizacji jednowymiarowego wektora positions buduje wizualizację na siatce 2D.
    
    - Dla każdego punktu na siatce (dla wnętrza, czyli punktów niebrzegowych) wyznaczany jest odpowiedni offset jako liczba spacji, co odpowiada przesunięciu poziomemu.
    - Następnie, dla każdego wiersza siatki (odpowiadającego współrzędnej pionowej) tworzona jest linia tekstowa z punktami wizualizowanymi jako kropki (lub znak specjalny na końcu wiersza).
    - Cała siatka jest wypisywana w terminalu, dzięki czemu uzyskujemy lepsze odzwierciedlenie przestrzeni (np. 13x13 punktów przy N = 20).
- **Parametry symulacji:**
    
    Parametry takie jak L (długość domeny), h_val (głębokość), T (okres) i kroki siatki delta_x, delta_z działają tak samo, jak w uproszczonej wersji – z tą różnicą, że wizualizacja jest bardziej "przestrzenna" i odpowiada rzeczywistemu rozmieszczeniu punktów na płaszczyźnie.
    

---