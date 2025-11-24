# Programa 3 – Factorización LU, Cholesky y Métodos de Potencias
**Materia:** Métodos Numéricos  
**Alumno:** Max Ontiveros  
**Lenguaje:** C  
**Archivos incluidos:** `programa3.c`

---

##  Contenido del programa

Este programa cumple con los requisitos solicitados para el Programa 3:

### ✔ Solución de ecuaciones
- Método de la Secante

### ✔ Solución de sistemas de ecuaciones
- Lectura de matriz y verificación de EDD
- Cálculo de determinante
- Método de Gauss–Seidel
- Factorización LU (Doolittle)
- Factorización de Cholesky

### ✔ Obtención de valores propios
- Método de potencias
- Método de potencias inversa
- Secuencia de aproximaciones
- Valor propio máximo y mínimo

---

##  Ejecución del programa en Colab

Para ejecutar el programa en Google Colab:

1. Abrir el siguiente enlace al notebook:  
👉 **(Aquí colocarás TU enlace de Colab)**

2. El notebook:
   - Monta Google Drive  
   - Copia `programa3.c`  
   - Compila el código con `gcc`  
   - Ejecuta el programa directamente

3. No requiere configuraciones adicionales.

---

## Cómo compilar manualmente (Linux/Mac/Colab)

```bash
gcc programa3.c -o programa3 -lm
./programa3
```

## Cómo compilar manualmente (Windows)

```bash
gcc programa3.c -o programa3.exe -lm
programa3.exe
```

---

## Estructura del repositorio

```
programa3-potencias/
│── programa3.c
│── README.md
```

---

## ✔ Requisitos cubiertos

- [x] Cholesky  
- [x] LU  
- [x] Gauss-Seidel  
- [x] Secante  
- [x] Potencias  
- [x] Potencia inversa  
- [x] Valores propios máximo y mínimo  
- [x] Todas las aproximaciones  
- [x] Menú integrado  
- [x] Funcionamiento directo  
