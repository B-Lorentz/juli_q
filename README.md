Kvantummechanika-A projektfeladat
Clebsch-Gordan együtthatók generálása
Bódy Lőrinc és Hajnal Dániel

# Telepítés:

-AlgebraicNumbers (https://github.com/anj1/AlgebraicNumbers.jl) installálása:

juliából pkg módban (ez a "]" gomb megnyomásával elérhető) kell a következő parancsot futtatni:

"add https://github.com/Snowgun/AlgebraicNumbers.jl.git"

Ez letölti az AlgebraicNumbers csomagot.
(az eredeti régebbi, és nincs benne project file, ehhez hozzátettük)

-Egyéb csomagok, szintén pkg módban kell telepíteni:

add PolynomialRoots
add Nemo
add Memoize

-juli_q letöltése:

terminálból a célmappábban a következő parancsot kell futtatni:

git clone https://github.com/B-Lorentz/juli_q.git


# Használat:

-Rekurzív algoritmus egy J állapot j1 és j2-re bontására: cg.jl
-Tetszőleges számú alrendszerre bontás: cg_unlimited.jl

