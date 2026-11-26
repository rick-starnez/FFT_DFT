import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# ---------- Generadores de señales ----------
def gen_square(n, N):
    """Onda cuadrada (50% duty)."""
    return np.where((np.arange(N) % n) < (n // 2), 1.0, -1.0)


def gen_sine(N, freq):
    """Seno de la longitud total N con `freq` ciclos."""
    t = np.linspace(0, 2 * np.pi * freq, N, endpoint=False)
    return np.sin(t)


def gen_triangle(N, freq):
    """Onda triangular con `freq` ciclos."""
    t = np.linspace(0, 1, N, endpoint=False)
    return 2 * np.abs(2 * (t * freq - np.floor(t * freq + 0.5))) - 1


def gen_random(N, seed=None):
    """Señal aleatoria."""
    rng = np.random.default_rng(seed)
    return rng.uniform(-1, 1, N)


# ---------- DFT directa (O(N^2)) ----------
def dft(x):
    N = len(x)
    n = np.arange(N)
    k = n.reshape((N, 1))
    W = np.exp(-2j * np.pi * k * n / N)
    return W @ x


# ---------- FFT iterativa (radix-2, Cooley-Tukey) ----------
def fft(x):
    x = np.asarray(x, dtype=complex)
    N = len(x)
    if N & (N - 1) != 0:
        raise ValueError("FFT radix-2 requiere N potencia de 2")

    # Bit-reversal
    j = 0
    for i in range(1, N):
        bit = N >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            x[i], x[j] = x[j], x[i]

    # Mariposas
    length = 2
    while length <= N:
        ang = -2j * np.pi / length
        wlen = np.exp(ang)
        for i in range(0, N, length):
            w = 1 + 0j
            for k in range(i, i + length // 2):
                u = x[k]
                v = x[k + length // 2] * w
                x[k] = u + v
                x[k + length // 2] = u - v
                w *= wlen
        length <<= 1
    return x


# ---------- Aplicación ----------
class DFTSimulator:
    def __init__(self, root):
        self.root = root
        self.root.title("Simulador DFT vs FFT")
        self.root.geometry("1400x800")

        # --- Barra superior con combo ---
        top = ttk.Frame(root, padding=8)
        top.pack(side=tk.TOP, fill=tk.X)

        ttk.Label(top, text="Señal:").pack(side=tk.LEFT, padx=(0, 6))

        self.signal_var = tk.StringVar(value="Cuadrada")
        combo = ttk.Combobox(
            top,
            textvariable=self.signal_var,
            values=["Cuadrada", "Seno", "Triangular", "Random"],
            state="readonly",
            width=15,
        )
        combo.pack(side=tk.LEFT)
        combo.bind("<<ComboboxSelected>>", lambda e: self.update_plots())

        ttk.Label(top, text="   N (potencia de 2):").pack(side=tk.LEFT, padx=(12, 6))
        self.N_var = tk.IntVar(value=64)
        N_spin = ttk.Spinbox(
            top, from_=8, to=4096, increment=8,
            textvariable=self.N_var, width=8,
            command=self.update_plots,
        )
        N_spin.pack(side=tk.LEFT)
        N_spin.bind("<Return>", lambda e: self.update_plots())

        # --- Contenedor de las dos columnas ---
        body = ttk.Frame(root)
        body.pack(fill=tk.BOTH, expand=True)

        # Columna izquierda: FFT
        left = ttk.LabelFrame(body, text="FFT (rápida)")
        left.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=6, pady=6)

        # Columna derecha: DFT
        right = ttk.LabelFrame(body, text="DFT (directa)")
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=6, pady=6)

        # Figuras de matplotlib
        self.fig_fft, self.ax_fft_sig, self.ax_fft_mag = self._make_figure()
        self.fig_dft, self.ax_dft_sig, self.ax_dft_mag = self._make_figure()

        self.canvas_fft = FigureCanvasTkAgg(self.fig_fft, master=left)
        self.canvas_fft.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.canvas_dft = FigureCanvasTkAgg(self.fig_dft, master=right)
        self.canvas_dft.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Tiempos mostrados
        self.time_label = ttk.Label(top, text="")
        self.time_label.pack(side=tk.RIGHT, padx=10)

        self.update_plots()

    def _make_figure(self):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 6))
        fig.tight_layout(pad=2.0)
        return fig, ax1, ax2

    def _make_signal(self, name, N):
        if name == "Cuadrada":
            return gen_square(N // 4, N)
        elif name == "Seno":
            return gen_sine(N, 5)
        elif name == "Triangular":
            return gen_triangle(N, 5)
        else:
            return gen_random(N, seed=42)

    def update_plots(self):
        name = self.signal_var.get()
        N = int(self.N_var.get())
        # Asegurar potencia de 2 (necesario para la FFT propia)
        N = 2 ** int(np.round(np.log2(max(8, N))))
        self.N_var.set(N)

        x = self._make_signal(name, N)

        # --- FFT ---
        import time
        t0 = time.perf_counter()
        X_fft = fft(x)
        t_fft = (time.perf_counter() - t0) * 1000  # ms

        # --- DFT ---
        t0 = time.perf_counter()
        X_dft = dft(x)
        t_dft = (time.perf_counter() - t0) * 1000  # ms

        # Verificación
        err = np.max(np.abs(X_fft - X_dft))

        self.time_label.config(
            text=f"FFT: {t_fft:.2f} ms   |   DFT: {t_dft:.2f} ms   |   err máx: {err:.2e}"
        )

        # Magnitudes normalizadas
        mag_fft = np.abs(X_fft) / N
        mag_dft = np.abs(X_dft) / N
        freqs = np.arange(N)

        # --- Plot FFT ---
        self.ax_fft_sig.clear()
        self.ax_fft_sig.plot(x, color="tab:blue")
        self.ax_fft_sig.set_title(f"Señal de entrada — {name} (N={N})")
        self.ax_fft_sig.grid(True, alpha=0.3)

        self.ax_fft_mag.clear()
        self.ax_fft_mag.stem(freqs[: N // 2], mag_fft[: N // 2], basefmt=" ")
        self.ax_fft_mag.set_title(f"|FFT(x)| normalizado — {t_fft:.2f} ms")
        self.ax_fft_mag.set_xlabel("k")
        self.ax_fft_mag.set_ylabel("Magnitud")
        self.ax_fft_mag.grid(True, alpha=0.3)

        # --- Plot DFT ---
        self.ax_dft_sig.clear()
        self.ax_dft_sig.plot(x, color="tab:orange")
        self.ax_dft_sig.set_title(f"Señal de entrada — {name} (N={N})")
        self.ax_dft_sig.grid(True, alpha=0.3)

        self.ax_dft_mag.clear()
        self.ax_dft_mag.stem(freqs[: N // 2], mag_dft[: N // 2], basefmt=" ",
                            linefmt="tab:orange", markerfmt="o")
        self.ax_dft_mag.set_title(f"|DFT(x)| normalizado — {t_dft:.2f} ms")
        self.ax_dft_mag.set_xlabel("k")
        self.ax_dft_mag.set_ylabel("Magnitud")
        self.ax_dft_mag.grid(True, alpha=0.3)

        self.fig_fft.tight_layout()
        self.fig_dft.tight_layout()
        self.canvas_fft.draw()
        self.canvas_dft.draw()


if __name__ == "__main__":
    root = tk.Tk()
    app = DFTSimulator(root)
    root.mainloop()