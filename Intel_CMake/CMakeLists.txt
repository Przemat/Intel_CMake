# CMakeList.txt: projekt CMake dla elementu Intel_CMake, dołącz źródło i zdefiniuj
# tutaj logikę specyficzną dla projektu.
#
cmake_minimum_required (VERSION 3.8)

# Dodaj źródło do pliku wykonywalnego tego projektu.
add_executable (Intel_CMake "Intel_CMake.cpp"  "Matrix.h")

# TODO: Dodaj testy i zainstaluj elementy docelowe w razie potrzeby.
