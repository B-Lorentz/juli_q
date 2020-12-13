**Kvantummechanika-A projektfeladat**
# Clebsch-Gordan együtthatók generálása

Bódy Lőrinc és Hajnal Dániel

## Telepítés:

1. juli_q letöltése: `$ git clone https://github.com/B-Lorentz/juli_q.git`

2. AlgebraicNumbers installálása: terminálból lépjünk be a juli_q mappába, ahol `$ git clone https://github.com/Snowgun/AlgebraicNumbers.jl.git` -t futtassuk, ezután juliából pkg módban (ez a "]" gomb megnyomásával elérhető) kell a következő parancsot futtatni: `add [path]/AlgebraicNumbers.jl`
ahol [path] az abszolút elérési útvonala az AlegbraicNumbers.jl-nek. (ezt nem mi írtuk, csak hiányzott belőle egy ún. projektfile, és ezt hozzáadtuk)

3. Egyéb csomagok, szintén pkg módban kell telepíteni: `add Memoize`

## Használat:

A `scripts` mappában vannak parancssorból futtatható verziók. Az törteket racionális alakban (pl.: `-3/2`) lehet beírni.

### Közvetlen lekérdezése egy CG-együtthatónak

`$ julia CG_calc.jl J M j1 m1 j2 m2`

Kiírja a megfelelő CG együthatót a rekurzió szerint.

ellenőrzés:
`$ julia CG_test.jl N pythonfile.py`

N különböző (J M j1 m1 j2 m2) 6-ost generál, kiszámítja a megfelelő CG-együtthatókat, és egy python file-t csinál, amely az eredményt összeveti a `sympy` csomag `sympy.physics.quantum.spin.CG` függvényének kimenetével, amely lényegében egy CG-táblázat.

Ezután a python file-t futtatva elvégződik az összehasonlítás (példakiement `output/cgtest.py`)

### Rekurzív algoritmus egy J,M állapot j1 és j2-re bontására

`$ julia decompose.jl J M j1 j2 filename.md`

kiírja a képernyőre a dekompozíciót a megfelelő együtthatókkal, a filename.md-be pedig egy olvashatóbb Markdown kód formájában kimenti. (példakiement `output/dec.md`)

ellenőrzés: 
`$ julia decompose_test.jl J M j1 j2 pythonfile.py`
(példakiement `output/dec_test.py`)

### J, M állapot 3 alrendszerre bontása:
`$ julia three_decompose.jl J M j1 j2 j3 filename.py`

Markdownként menti külön egyre normálva a különféle $J_{12}$ résszösszegeket feltételező eseteket.

(példakiement `output/3d_test.md`)

## 2.1. Mi lehet a három komponensű rendszer teljes impulzusmomentuma?

Két komponensű rendszer teljes impulzusmomentuma

$J_{12}$
