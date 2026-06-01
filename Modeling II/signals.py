import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

DT = 1e-4
FS = 1.0 / DT
INPUT_FILE = "signal_7.npy"
DETECTION_FACTOR = 6.0
MIN_PEAK_DISTANCE_HZ = 3.0
MAX_PEAKS = 20
BAND_HALF_WIDTH_BINS = 2
def load_signal(filename: str) -> np.ndarray:
    x = np.load(filename)
    if x.ndim != 1:
        x = np.ravel(x)
    return x.astype(float)

def local_maxima_indices(a: np.ndarray) -> np.ndarray:
    if len(a) < 3:
        return np.array([], dtype=int)
    return np.where((a[1:-1] > a[:-2]) & (a[1:-1] >= a[2:]))[0] + 1


def robust_threshold(amp: np.ndarray) -> float:
    start = max(5, len(amp) // 5000)
    base = amp[start:]

    med = np.median(base)
    mad = np.median(np.abs(base - med))
    sigma_robust = 1.4826 * mad
    return med + DETECTION_FACTOR * sigma_robust


def select_peaks(freqs: np.ndarray, amp: np.ndarray) -> tuple[np.ndarray, float]:
    threshold = robust_threshold(amp)

    cand = local_maxima_indices(amp)
    cand = cand[amp[cand] > threshold]

    if len(cand) == 0:
        cand = np.array([np.argmax(amp[1:]) + 1], dtype=int)
    cand = cand[np.argsort(amp[cand])[::-1]]

    df = freqs[1] - freqs[0]
    min_bins = max(1, int(round(MIN_PEAK_DISTANCE_HZ / df)))

    selected = []
    for idx in cand:
        if all(abs(idx - j) >= min_bins for j in selected):
            selected.append(idx)
        if len(selected) >= MAX_PEAKS:
            break

    return np.array(sorted(selected), dtype=int), threshold

def build_filter_mask(n_spec: int, peak_idx: np.ndarray, half_width_bins: int) -> np.ndarray:
    mask = np.zeros(n_spec, dtype=bool)
    for idx in peak_idx:
        lo = max(0, idx - half_width_bins)
        hi = min(n_spec - 1, idx + half_width_bins)
        mask[lo:hi + 1] = True
    return mask

def main():
    path = Path(INPUT_FILE)
    if not path.exists():
        raise FileNotFoundError(f"Файл {INPUT_FILE} не найден в текущей папке.")
    x = load_signal(INPUT_FILE)
    N = len(x)
    t = np.arange(N) * DT

    print(f"N = {N}")
    print(f"dt = {DT:.6e} s, fs = {FS:.2f} Hz")
    print(f"Длина записи T = {N * DT:.6f} s")

    x_mean = np.mean(x)
    x0 = x - x_mean
    window = np.hanning(N)
    x_det = x0 * window
    X_det = np.fft.rfft(x_det)
    freqs = np.fft.rfftfreq(N, d=DT)

    amp_det = (2.0 / np.sum(window)) * np.abs(X_det)
    amp_det[0] = np.abs(X_det[0]) / np.sum(window)
    if N % 2 == 0:
        amp_det[-1] = np.abs(X_det[-1]) / np.sum(window)
    peak_idx, threshold = select_peaks(freqs, amp_det)
    peak_freqs = freqs[peak_idx]

    print("\nНайденные частоты полезного сигнала (Гц):")
    for f in peak_freqs:
        print(f"  {f:.6f}")

    np.savetxt("detected_frequencies.txt", peak_freqs, fmt="%.10f")
    X = np.fft.rfft(x0)
    mask = build_filter_mask(len(X), peak_idx, BAND_HALF_WIDTH_BINS)
    mask[0] = False

    X_filt = np.zeros_like(X)
    X_filt[mask] = X[mask]
    x_rec0 = np.fft.irfft(X_filt, n=N)
    x_rec = x_rec0 + x_mean
    mse = np.mean((x - x_rec) ** 2)
    print(f"\nMSE восстановления: {mse:.6e}")
  
    n_view = min(N, 5000)
    plt.figure(figsize=(10, 5))
    plt.plot(t[:n_view], x[:n_view], label="Зашумленный сигнал", linewidth=1)
    plt.plot(t[:n_view], x_rec[:n_view], label="Восстановленный сигнал", linewidth=2)
    plt.xlabel("Время, с")
    plt.ylabel("Амплитуда")
    plt.title("Фрагмент сигнала до и после фильтрации")
    plt.legend()
    plt.tight_layout()
    plt.savefig("time_signal_comparison.png", dpi=200)
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.plot(freqs, amp_det, label="Амплитудный спектр", linewidth=1)
    plt.axhline(threshold, color="red", linestyle="--", label="Порог обнаружения")

    if len(peak_idx) > 0:
        plt.scatter(freqs[peak_idx], amp_det[peak_idx], color="black", zorder=3, label="Найденные пики")

    plt.xlim(0, FS / 2)
    plt.xlabel("Частота, Гц")
    plt.ylabel("Амплитуда")
    plt.title("Спектр сигнала и выделенные частоты")
    plt.legend()
    plt.tight_layout()
    plt.savefig("spectrum_peaks.png", dpi=200)
    plt.show()

    amp_full = (2.0 / N) * np.abs(X)
    amp_filt = (2.0 / N) * np.abs(X_filt)
    amp_full[0] = np.abs(X[0]) / N
    amp_filt[0] = np.abs(X_filt[0]) / N
    if N % 2 == 0:
        amp_full[-1] = np.abs(X[-1]) / N
        amp_filt[-1] = np.abs(X_filt[-1]) / N

    plt.figure(figsize=(10, 5))
    plt.plot(freqs, amp_full, label="Исходный спектр", linewidth=1)
    plt.plot(freqs, amp_filt, label="Отфильтрованный спектр", linewidth=2)
    plt.xlim(0, FS / 2)
    plt.xlabel("Частота, Гц")
    plt.ylabel("Амплитуда")
    plt.title("Сравнение спектров до и после фильтрации")
    plt.legend()
    plt.tight_layout()
    plt.savefig("spectrum_comparison.png", dpi=200)
    plt.show()

    plt.figure(figsize=(10, 4))
    plt.plot(t, x, label="Зашумленный", linewidth=0.8)
    plt.plot(t, x_rec, label="Восстановленный", linewidth=1.2)
    plt.xlabel("Время, с")
    plt.ylabel("Амплитуда")
    plt.title("Полный сигнал")
    plt.legend()
    plt.tight_layout()
    plt.savefig("full_signal_comparison.png", dpi=200)
    plt.show()

    np.save("signal_7_filtered.npy", x_rec)

    print("\nГотово.")
    print("Сохранены файлы:")
    print("  detected_frequencies.txt")
    print("  signal_7_filtered.npy")
    print("  time_signal_comparison.png")
    print("  spectrum_peaks.png")  
    print("  spectrum_comparison.png")
    print("  full_signal_comparison.png")


if __name__ == "__main__":
    main()
