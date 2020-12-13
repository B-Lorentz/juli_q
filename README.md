Kvantummechanika-A projektfeladat
Clebsch-Gordan együtthatók generálása
Bódy Lőrinc és Hajnal Dániel

# Telepítés:

-juli_q letöltése:

terminálból a következő parancsot kell futtatni, a célmappából:

git clone https://github.com/B-Lorentz/juli_q.git

-AlgebraicNumbers installálása:

terminálból lépjünk be a juli_q mappába, majd futtassuk a következő parancsot:

git clone https://github.com/anj1/AlgebraicNumbers.jl

juliából pkg módban (ez a "]" gomb megnyomásával elérhető) kell a következő parancsot futtatni:

add [path]/AlgebraicNumbers.jl

ahol [path] az abszolút elérési útvonala az AlegbraicNumbers.jl-nek (ezt a pwd futtatásával lehet megkapni)

Ez letölti az AlgebraicNumbers csomagot.
(az eredeti régebbi, és nincs benne project file, ehhez hozzátettük)

-Egyéb csomagok, szintén pkg módban kell telepíteni:

add PolynomialRoots
add Nemo
add Memoize




# Használat:

## Rekurzív algoritmus egy J,M állapot j1 és j2-re bontására: decompose.jl

szignatúra: julia decompose.jl J M j1 j2 filename.md
kiírja a képernyőre a dekompozíciót a megfelelő együtthatókkal, a filename.md-be pedig egy olvashatóbb Markdown kód formájában kimenti

ellenőrzés: julia decompose_test.jl J M j1 j2 pythonfile.py
a python sympy csomagjában van egy CG táblázat, azzal összehasonlítja az általunk számoltakat

## J, M állapot 3 alrendszerre bontása:

szignatúra: julia three_decompose.jl J M j1 j2 j3 filename.py
kiírja a képernyőre a dekompozíciót a megfelelő együtthatókkal, a filename.py-be pedig egy olvashatóbb Jupyter Notebook formájában kimenti

ellenőrzés: julia decompose_test.jl J M j1 j2 pythonfile.py
a python sympy csomagjában van egy CG táblázat, azzal összehasonlítja az általunk számoltakat


