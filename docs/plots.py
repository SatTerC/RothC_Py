from datetime import datetime
from pathlib import Path

from rothc_py.modernc import percent_modern_c

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Error: Matplotlib is not installed. Exiting.")
    raise

docs_dir = Path(__file__).parent


def plot_percent_modern_c() -> None:
    start_date = datetime(1900, 1, 1)
    n_months = 150 * 12  # up to 2050

    dates = [
        datetime(start_date.year + (i // 12), start_date.month + (i % 12), 1)
        for i in range(n_months)
    ]
    values = percent_modern_c(start_date, n_months)

    fig, ax = plt.subplots()
    ax.set_xlabel("dates")
    ax.set_ylabel("% modern C")

    ax.plot(dates, values, label="% modern C")
    ax.axvline(
        datetime(1939, 1, 1), color="orange", linestyle=":", label="First observations"
    )
    ax.axvline(
        datetime(2007, 12, 1), color="orange", linestyle=":", label="Last observations"
    )

    ax.legend()

    target_location = docs_dir / "modernc.png"
    fig.savefig(target_location)


if __name__ == "__main__":
    plot_percent_modern_c()
