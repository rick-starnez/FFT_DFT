# Simulador: DFT vs FFT

Aplicación de escritorio en Python para comparar visualmente la Transformada Discreta de Fourier (DFT) calculada de forma directa contra la Transformada Rápida de Fourier (FFT) calculada con el algoritmo de Cooley-Tukey.

## Descripción

El simulador muestra dos columnas lado a lado:

- Columna izquierda: FFT (Transformada Rápida de Fourier)
  - Señal de entrada en el dominio del tiempo
  - Espectro de magnitud de la FFT
  - Tiempo de cómputo medido

- Columna derecha: DFT (Transformada Discreta de Fourier)
  - Señal de entrada en el dominio del tiempo
  - Espectro de magnitud de la DFT
  - Tiempo de cómputo medido

El objetivo didáctico es que el alumno:

- Observe que la DFT y la FFT producen exactamente el mismo espectro.
- Compare los tiempos de cómputo y entienda por qué la FFT es más rápida
  (O(N log N) vs O(N^2)).
- Vea cómo crece la diferencia de tiempos al aumentar N.

## Señales disponibles

- Cuadrada
- Seno
- Triangular
- Random (aleatoria con semilla fija para reproducibilidad)

## Fórmulas e implementación

DFT directa (definición):

    X[k] = sum_{n=0}^{N-1} x[n] * exp(-j*2*pi*k*n/N),   k = 0..N-1

Implementada en la aplicación como multiplicación matricial O(N^2):

    W[k,n] = exp(-j*2*pi*k*n/N)
    X = W @ x

FFT (Cooley-Tukey, radix-2, iterativa):

- Requiere N potencia de 2.
- Reordena los índices con bit-reversal.
- Aplica mariposas por etapas con factores twiddle.
- Complejidad O(N log N).

El espectro de magnitud mostrado es |X[k]| / N para k = 0..N/2-1,
es decir, la mitad del espectro (hasta la frecuencia de Nyquist).

## Controles de la interfaz

- Combo "Señal": selecciona el tipo de señal de entrada.
- Spinbox "N": número de muestras (potencia de 2 recomendada; el programa
  ajusta automáticamente a la potencia de 2 más cercana para la FFT).

Al cambiar cualquier control, ambas columnas se recalculan y redibujan
automáticamente.

## Verificación

En la barra superior se muestra el error máximo absoluto entre los
coeficientes producidos por la FFT y por la DFT. Debe estar en el orden
de 1e-12 o menor (limitado por la precisión de punto flotante), lo cual
confirma que ambos algoritmos son matemáticamente equivalentes.

## Requisitos

- Python 3.9 o superior
- numpy
- matplotlib
- tkinter (incluido en la instalación estándar de Python en Windows y macOS;
  en Linux puede requerir el paquete python3-tk del sistema)

Las versiones recomendadas están en el archivo requirements.txt.

## Instalación

1. Clonar el repositorio:

       git clone <URL_DEL_REPOSITORIO>
       cd <NOMBRE_DEL_REPOSITORIO>

2. (Recomendado) Crear y activar un entorno virtual:

       python -m venv venv

   En macOS/Linux:

       source venv/bin/activate

   En Windows:

       venv\Scripts\activate

3. Instalar dependencias:

       pip install -r requirements.txt

## Ejecución

    python DFT_FFT.py

Se abrirá una ventana con las dos columnas del simulador.

## Nota para Linux

Si al ejecutar aparece un error del tipo "No module named tkinter",
instalar el paquete del sistema correspondiente:

- Debian/Ubuntu:

      sudo apt install python3-tk

- Fedora:

      sudo dnf install python3-tkinter

- Arch:

      sudo pacman -S tk

## Estructura del proyecto

    .
    ├── dft_vs_fft.py         # Aplicación principal
    ├── requirements.txt      # Dependencias de Python
    └── README.md             # Este archivo

## Uso didáctico sugerido

1. Seleccionar la señal "Seno" con N pequeño (por ejemplo 64) y observar
   que ambas columnas muestran el mismo espectro (dos picos claros).
2. Aumentar N progresivamente (128, 256, 512, 1024, 2048) y comparar los
   tiempos de cómputo mostrados:
   - El tiempo de la DFT crece aproximadamente como N^2.
   - El tiempo de la FFT crece aproximadamente como N log N.
3. Con N = 2048 o mayor, la diferencia entre ambos tiempos se hace muy
   notoria, evidenciando por qué la FFT es el algoritmo estándar en la
   práctica.
4. Cambiar entre las distintas señales y observar cómo cambia la
   distribución de energía en el espectro.
5. Verificar que el error máximo entre DFT y FFT es prácticamente cero
   (del orden de 1e-12), confirmando la equivalencia matemática.

## Licencia

Uso libre con fines educativos.
