Kvantummechanika-A projektfeladat
Clebsch-Gordan együtthatók generálása
Bódy Lőrinc és Hajnal Dániel

## Telepítés:

1. juli_q letöltése: `git clone https://github.com/B-Lorentz/juli_q.git`

2. AlgebraicNumbers installálása: terminálból lépjünk be a juli_q mappába, ahol `git clone https://github.com/anj1/AlgebraicNumbers.jl`, ezután juliából pkg módban (ez a "]" gomb megnyomásával elérhető) kell a következő parancsot futtatni: `add [path]/AlgebraicNumbers.jl`
ahol [path] az abszolút elérési útvonala az AlegbraicNumbers.jl-nek.

3. Egyéb csomagok, szintén pkg módban kell telepíteni:

`add PolynomialRoots`
`add Nemo`
`add Memoize`

## Használat:

A `scripts` mappában vannak parancssorból futtatható verziók. Az törteket racionáli alakban (pl.: `-3/2`) lehet beírni.

### Közvetlen lekérdezése egy CG-együtthatónak

`$ julia CG_calc.jl J M j1 m1 j2 m2`

Kiírja a megfelelő CG együthatót a rekurzió szerint.

ellenőrzés:
` julia CG_test.jl N pythonfile.py`

### Rekurzív algoritmus egy J,M állapot j1 és j2-re bontására

szignatúra: julia decompose.jl J M j1 j2 filename.md
kiírja a képernyőre a dekompozíciót a megfelelő együtthatókkal, a filename.md-be pedig egy olvashatóbb Markdown kód formájában kimenti

ellenőrzés: julia decompose_test.jl J M j1 j2 pythonfile.py
a python sympy csomagjában van egy CG táblázat, azzal összehasonlítja az általunk számoltakat

## J, M állapot 3 alrendszerre bontása:

szignatúra: julia three_decompose.jl J M j1 j2 j3 filename.py
kiírja a képernyőre a dekompozíciót a megfelelő együtthatókkal, a filename.py-be pedig egy olvashatóbb Jupyter Notebook formájában kimenti

ellenőrzés: julia decompose_test.jl J M j1 j2 pythonfile.py
a python sympy csomagjában van egy CG táblázat, azzal összehasonlítja az általunk számoltakat


