/* Programa 3. Método de potencias

   INTEGRANTES:
   - Albaran Galindo Luis Sebastian
   - León Ramirez Damian 
   - Ontiveros Martínez Maximiliano
   - Perez García Éden Yohualli Tonatiuh
   

   
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ---------- UTILIDADES DE MEMORIA Y I/O ---------- */

double **crearMatriz(int n) {
    double **M = (double**) malloc(n * sizeof(double*));
    if (!M) { perror("malloc"); exit(EXIT_FAILURE); }
    for (int i = 0; i < n; ++i) {
        M[i] = (double*) calloc(n, sizeof(double));
        if (!M[i]) { perror("calloc"); exit(EXIT_FAILURE); }
    }
    return M;
}

void liberarMatriz(double **M, int n) {
    if (!M) return;
    for (int i = 0; i < n; ++i) free(M[i]);
    free(M);
}

double *crearVector(int n) {
    double *v = (double*) calloc(n, sizeof(double));
    if (!v) { perror("calloc"); exit(EXIT_FAILURE); }
    return v;
}

void liberarVector(double *v) { if (v) free(v); }

void imprimirMatriz(double **M, int n) {
    for (int i = 0; i < n; ++i) {
        printf("| ");
        for (int j = 0; j < n; ++j) printf("%11.6f ", M[i][j]);
        printf("|\n");
    }
}

void imprimirVector(double *v, int n) {
    printf("[");
    for (int i = 0; i < n; ++i) {
        printf("%.6f", v[i]);
        if (i < n-1) printf(", ");
    }
    printf("]\n");
}

/* ---------- FUNCION PARA EL METODO DE LA SECANTE ----------
   Cambiar la funcion f(x) segun el problema.
*/
double f(double x) {
    // EJEMPLO: f(x) = x^2 - 5x + 6  -> raices 2 y 3
    // Cambia esta funcion segun tu ejercicio.
    return x*x - 5.0*x + 6.0;
}

double metodo_secante(double x0, double x1, double tol, int maxIter) {
    double f0 = f(x0), f1 = f(x1);
    double x2 = x1;

    printf("\n--- Metodo de la Secante ---\n");
    printf("Iter\t   x_k\t\t   f(x_k)\n");
    printf("0\t% .10f\t% .6e\n", x0, f0);
    printf("1\t% .10f\t% .6e\n", x1, f1);

    for (int k = 2; k <= maxIter; ++k) {
        if (fabs(f1 - f0) < 1e-16) {
            printf("Denominador demasiado pequeno (f1 - f0 ~ 0). Abortando.\n");
            return x1;
        }
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        double f2 = f(x2);
        printf("%d\t% .10f\t% .6e\n", k-1, x2, f2);

        if (fabs(x2 - x1) < tol) {
            printf("Convergio en %d iteraciones.\n", k-1);
            return x2;
        }
        x0 = x1; f0 = f1;
        x1 = x2; f1 = f2;
    }

    printf("No convergio en %d iteraciones, ultima aproximacion = % .10f\n", maxIter, x2);
    return x2;
}

/* ---------- VALIDACIONES: EDD y SIMETRIA ---------- */

int es_diagonal_dominante(double **M, int n) {
    for (int i = 0; i < n; ++i) {
        double suma = 0.0;
        for (int j = 0; j < n; ++j) if (j != i) suma += fabs(M[i][j]);
        if (fabs(M[i][i]) <= suma) return 0;
    }
    return 1;
}

int es_simetrica(double **M, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j)
            if (fabs(M[i][j] - M[j][i]) > 1e-9) return 0;
    return 1;
}

/* ---------- DETERMINANTE (Eliminacion Gaussiana con pivoteo parcial) ---------- */

double determinante(double **A, int n) {
    double **T = crearMatriz(n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            T[i][j] = A[i][j];

    double det = 1.0;

    for (int i = 0; i < n; ++i) {
        int piv = i;
        for (int r = i+1; r < n; ++r)
            if (fabs(T[r][i]) > fabs(T[piv][i])) piv = r;

        if (fabs(T[piv][i]) < 1e-16) {
            liberarMatriz(T, n);
            return 0.0;
        }

        if (piv != i) {
            double *tmp = T[i]; T[i] = T[piv]; T[piv] = tmp;
            det = -det;
        }

        det *= T[i][i];

        for (int r = i+1; r < n; ++r) {
            double factor = T[r][i] / T[i][i];
            for (int c = i; c < n; ++c) T[r][c] -= factor * T[i][c];
        }
    }

    liberarMatriz(T, n);
    return det;
}

/* ---------- SUSTITUCION ADELANTE Y ATRAS ---------- */

double *sustitucion_adelante(double **L, double *b, int n) {
    double *y = crearVector(n);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < i; ++j) s += L[i][j] * y[j];
        y[i] = (b[i] - s) / L[i][i];
    }
    return y;
}

double *sustitucion_atras(double **U, double *y, int n) {
    double *x = crearVector(n);
    for (int i = n-1; i >= 0; --i) {
        double s = 0.0;
        for (int j = i+1; j < n; ++j) s += U[i][j] * x[j];
        x[i] = (y[i] - s) / U[i][i];
    }
    return x;
}

/* ---------- CHOLESKY ---------- */

double **cholesky(double **A, int n) {
    double **L = crearMatriz(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double s = 0;
            for (int k = 0; k < j; ++k) s += L[i][k]*L[j][k];
            if (i == j) {
                double val = A[i][i] - s;
                if (val <= 0) { liberarMatriz(L,n); return NULL; } // no definida positiva
                L[i][i] = sqrt(val);
            } else {
                L[i][j] = (A[i][j] - s) / L[j][j];
            }
        }
    }
    return L;
}

/* ---------- DOOLITTLE (LU) ---------- */

int doolittle(double **A, double **L, double **U, int n) {
    // L and U must be allocated
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            L[i][j] = 0.0;
            U[i][j] = 0.0;
        }
    }
    for (int i = 0; i < n; ++i) L[i][i] = 1.0;

    for (int k = 0; k < n; ++k) {
        // U row k
        for (int j = k; j < n; ++j) {
            double s = 0.0;
            for (int p = 0; p < k; ++p) s += L[k][p] * U[p][j];
            U[k][j] = A[k][j] - s;
        }
        if (fabs(U[k][k]) < 1e-16) return 0; // fallo (pivote cero)
        // L column k
        for (int i = k+1; i < n; ++i) {
            double s = 0.0;
            for (int p = 0; p < k; ++p) s += L[i][p] * U[p][k];
            L[i][k] = (A[i][k] - s) / U[k][k];
        }
    }
    return 1;
}

/* ---------- MULTIPLICACION MATRIZ- VECTOR ---------- */

double *mult_mat_vec(double **M, double *v, int n) {
    double *r = crearVector(n);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < n; ++j) s += M[i][j] * v[j];
        r[i] = s;
    }
    return r;
}

/* ---------- METODO DE POTENCIAS (CON SUCESION) ---------- */

void metodo_potencias(double **A, int n, double tol, int maxIter, double *out_lambda, double **out_vec) {
    double *x = crearVector(n);
    for (int i = 0; i < n; ++i) x[i] = 1.0; // init

    double lambda = 0.0, lambdaAnt = 0.0;

    printf("\n--- Metodo de Potencias (sucesion) ---\n");
    for (int iter = 1; iter <= maxIter; ++iter) {
        double *y = mult_mat_vec(A, x, n);

        // aproximacion lambda: elemento de mayor magnitud en y
        double maxAbs = 0.0;
        for (int i = 0; i < n; ++i) if (fabs(y[i]) > maxAbs) { maxAbs = fabs(y[i]); lambda = y[i]; }

        // normalizar x = y / lambda
        if (fabs(lambda) < 1e-16) {
            printf("lambda aproximado ~ 0. Interrumpiendo.\n");
            liberarVector(y); break;
        }
        for (int i = 0; i < n; ++i) x[i] = y[i] / lambda;

        // imprimir la aproximacion actual (lambda y vector)
        printf("Iter %3d: lambda aprox = %.10f  | vector aprox = ", iter, lambda);
        imprimirVector(x, n);

        if (iter > 1 && fabs(lambda - lambdaAnt) < tol) {
            printf("Convergio en iteracion %d\n", iter);
            liberarVector(y);
            break;
        }
        lambdaAnt = lambda;
        liberarVector(y);
    }

    *out_lambda = lambda;
    *out_vec = x; // devolver vector (liberar fuera)
}

/* ---------- METODO DE POTENCIA INVERSA (CON SUCESION) ---------- */

void metodo_potencia_inversa(double **A, int n, double tol, int maxIter, double *out_lambda_min, double **out_vec_min) {
    // descomponer A = L*U (Doolittle)
    double **L = crearMatriz(n);
    double **U = crearMatriz(n);
    int ok = doolittle(A, L, U, n);
    if (!ok) {
        printf("LU no pudo factorizarse (posible pivote nulo). Metodo inverso no disponible.\n");
        liberarMatriz(L,n); liberarMatriz(U,n);
        *out_lambda_min = 0.0; *out_vec_min = NULL;
        return;
    }

    double *x = crearVector(n);
    for (int i = 0; i < n; ++i) x[i] = 1.0;

    double lamInv = 0.0, lamInvAnt = 0.0;

    printf("\n--- Metodo de Potencia Inversa (sucesion) ---\n");
    for (int iter = 1; iter <= maxIter; ++iter) {
        // resolver A * y = x  ->  (L U) y = x
        double *z = sustitucion_adelante(L, x, n);
        double *y = sustitucion_atras(U, z, n);

        // encontrar max en y
        double maxAbs = 0.0;
        for (int i = 0; i < n; ++i) if (fabs(y[i]) > maxAbs) { maxAbs = fabs(y[i]); lamInv = y[i]; }

        if (fabs(lamInv) < 1e-16) {
            printf("lamInv ~ 0. Interrumpiendo.\n");
            liberarVector(z); liberarVector(y); break;
        }

        // normalizar x = y / lamInv
        for (int i = 0; i < n; ++i) x[i] = y[i] / lamInv;

        double eigenA = 1.0 / lamInv;
        printf("Iter %3d: lambda(min) aprox = %.10f  | vector aprox = ", iter, eigenA);
        imprimirVector(x, n);

        if (iter > 1 && fabs(lamInv - lamInvAnt) < tol) {
            printf("Convergio en iteracion %d\n", iter);
            liberarVector(z); liberarVector(y);
            break;
        }
        lamInvAnt = lamInv;
        liberarVector(z); liberarVector(y);
    }

    *out_lambda_min = 1.0 / lamInv;
    *out_vec_min = x;

    liberarMatriz(L, n);
    liberarMatriz(U, n);
}

/* ---------- GAUSS-SEIDEL ---------- */

void gauss_seidel(double **A, double *b, int n, int maxIter, double tol) {
    double *x = crearVector(n);
    double *xOld = crearVector(n); // inicial zeros

    printf("\n--- Metodo Gauss-Seidel ---\n");

    for (int k = 1; k <= maxIter; ++k) {
        for (int i = 0; i < n; ++i) xOld[i] = x[i];

        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int j = 0; j < n; ++j) if (j != i) s += A[i][j] * x[j];
            if (fabs(A[i][i]) < 1e-16) { printf("0 en diagonal, no se puede continuar.\n"); liberarVector(x); liberarVector(xOld); return; }
            x[i] = (b[i] - s) / A[i][i];
        }

        double err = 0.0;
        for (int i = 0; i < n; ++i) if (fabs(x[i] - xOld[i]) > err) err = fabs(x[i] - xOld[i]);

        printf("Iter %3d: error = %.6e  | x = ", k, err);
        imprimirVector(x, n);

        if (err < tol) { printf("Convergio en iteracion %d\n", k); break; }
    }

    printf("Resultado final (Gauss-Seidel):\n");
    imprimirVector(x, n);
    liberarVector(x); liberarVector(xOld);
}

/* ---------- LECTURA DE MATRIZ Y VECTOR ---------- */

void leer_matriz_y_vector(double ***Amat, double **bvec, int *pn) {
    int n;
    printf("Ingrese n (tamanio de la matriz nxn): ");
    if (scanf("%d", &n) != 1) { printf("Entrada invalida.\n"); exit(EXIT_FAILURE); }

    double **A = crearMatriz(n);
    double *b = crearVector(n);

    printf("Ingrese los elementos de la matriz A (fila por fila):\n");
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            printf("A[%d][%d] = ", i+1, j+1);
            if (scanf("%lf", &A[i][j]) != 1) { printf("Entrada invalida.\n"); exit(EXIT_FAILURE); }
        }

    printf("Ingrese los elementos del vector b:\n");
    for (int i = 0; i < n; ++i) {
        printf("b[%d] = ", i+1);
        if (scanf("%lf", &b[i]) != 1) { printf("Entrada invalida.\n"); exit(EXIT_FAILURE); }
    }

    *Amat = A;
    *bvec = b;
    *pn = n;
}

/* ---------- PROGRAMA PRINCIPAL Y MENUS ---------- */

int main(void) {
    double **A = NULL;
    double *b = NULL;
    int n = 0;

    int op;
    do {
        printf("\n===== MENU PRINCIPAL =====\n");
        printf("1. Solucion de ecuaciones (Metodo de la Secante)\n");
        printf("2. Solucion de sistemas de ecuaciones\n");
        printf("3. Obtencion de valores propios (Potencias / Potencia inversa)\n");
        printf("4. Salir\n");
        printf("Seleccione opcion: ");
        if (scanf("%d", &op) != 1) { printf("Entrada invalida.\n"); return 0; }

        if (op == 1) {
            double x0, x1, tol;
            int maxIter;
            printf("\n--- Metodo de la Secante ---\n");
            printf("Nota: Modifique la funcion f(x) en el codigo si lo requiere.\n");
            printf("Ingrese x0: "); scanf("%lf", &x0);
            printf("Ingrese x1: "); scanf("%lf", &x1);
            printf("Ingrese tolerancia (ej. 1e-6): "); scanf("%lf", &tol);
            printf("Ingrese maximo de iteraciones: "); scanf("%d", &maxIter);
            double raiz = metodo_secante(x0, x1, tol, maxIter);
            printf("Raiz aproximada: %.10f\n", raiz);
        }

        else if (op == 2) {
            int sub;
            do {
                printf("\n--- Submenu: Sistemas de Ecuaciones ---\n");
                printf("1. Lectura de la matriz (Validar EDD / Determinante)\n");
                printf("2. Solucion del sistema (Gauss-Seidel)\n");
                printf("3. Factorizacion LU (Cholesky si aplica)\n");
                printf("4. Volver al menu principal\n");
                printf("Seleccione opcion: ");
                if (scanf("%d", &sub) != 1) { printf("Entrada invalida.\n"); return 0; }

                if (sub == 1) {
                    if (A) { liberarMatriz(A, n); liberarVector(b); A = NULL; b = NULL; n = 0; }
                    leer_matriz_y_vector(&A, &b, &n);
                    printf("\nMatriz A:\n"); imprimirMatriz(A, n);
                    if (es_diagonal_dominante(A, n)) {
                        printf("La matriz ES estrictamente diagonal dominante (EDD).\n");
                    } else {
                        printf("La matriz NO es EDD. Calculando determinante...\n");
                        double det = determinante(A, n);
                        printf("Determinante = %.10f\n", det);
                        if (fabs(det) < 1e-9) printf("ADVERTENCIA: El sistema puede no tener solucion unica (Det ~ 0).\n");
                    }
                }

                else if (sub == 2) {
                    if (!A || !b) { printf("Primero debe leer la matriz y el vector (opcion 2.1).\n"); continue; }
                    int maxIter = 100;
                    double tol = 1e-4;
                    printf("Ingrese tolerancia (ej. 1e-4): "); scanf("%lf", &tol);
                    printf("Ingrese maximo de iteraciones: "); scanf("%d", &maxIter);
                    gauss_seidel(A, b, n, maxIter, tol);
                }

                else if (sub == 3) {
                    if (!A || !b) { printf("Primero debe leer la matriz y el vector (opcion 2.1).\n"); continue; }
                    if (!es_simetrica(A, n)) {
                        printf("La matriz NO es simetrica. Cholesky requiere matriz simetrica definida positiva.\n");
                    } else {
                        double **L = cholesky(A, n);
                        if (!L) {
                            printf("Cholesky fallo: la matriz puede no ser definida positiva.\n");
                        } else {
                            printf("Descomposicion A = L L^T (Matriz L):\n");
                            imprimirMatriz(L, n);
                            double **Lt = crearMatriz(n);
                            for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) Lt[i][j] = L[j][i];
                            double *y = sustitucion_adelante(L, b, n);
                            double *x = sustitucion_atras(Lt, y, n);
                            printf("Sustitucion hacia adelante (Ly = b) => y = "); imprimirVector(y, n);
                            printf("Sustitucion hacia atras (L^T x = y) => solucion x = "); imprimirVector(x, n);
                            liberarVector(y); liberarVector(x);
                            liberarMatriz(L, n); liberarMatriz(Lt, n);
                        }
                    }
                    // Adicional: mostrar LU (Doolittle) tambien como extra
                    double **Llu = crearMatriz(n);
                    double **Ulu = crearMatriz(n);
                    if (doolittle(A, Llu, Ulu, n)) {
                        printf("\nDescomposicion LU (Doolittle):\nL =\n"); imprimirMatriz(Llu, n);
                        printf("U =\n"); imprimirMatriz(Ulu, n);
                    } else {
                        printf("\nDoolittle fallo (pivote nulo o degenerate).\n");
                    }
                    liberarMatriz(Llu, n); liberarMatriz(Ulu, n);
                }

            } while (sub != 4);
        }

        else if (op == 3) {
            if (!A) {
                printf("Primero debe ingresar la matriz A (opcion 2.1 en Sistemas) para calcular valores propios.\n");
                continue;
            }
            double tol; int maxIter;
            printf("Ingrese tolerancia (ej. 1e-6): "); scanf("%lf", &tol);
            printf("Ingrese maximo de iteraciones: "); scanf("%d", &maxIter);

            double lambda_max; double *vec_max = NULL;
            metodo_potencias(A, n, tol, maxIter, &lambda_max, &vec_max);

            double lambda_min; double *vec_min = NULL;
            metodo_potencia_inversa(A, n, tol, maxIter, &lambda_min, &vec_min);

            printf("\n=== RESULTADOS FINALES DE VALORES PROPIOS ===\n");
            printf("Valor propio dominante (maximo) aproximado: %.10f\n", lambda_max);
            if (vec_max) { printf("Vector propio asociado (aprox): "); imprimirVector(vec_max, n); liberarVector(vec_max); }

            printf("Valor propio minimo (aprox): %.10f\n", lambda_min);
            if (vec_min) { printf("Vector propio asociado (aprox): "); imprimirVector(vec_min, n); liberarVector(vec_min); }
        }

    } while (op != 4);

    printf("Saliendo...\n");
    if (A) liberarMatriz(A, n);
    if (b) liberarVector(b);
    return 0;
}
