"""Shared helpers for the Exercise 3 scripts."""
import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name), dpi=140)
    plt.close(fig)
    print("  wrote figures/" + name)


def check(label, ok):
    print(f"  [{'OK' if ok else 'FAIL'}] {label}")
    assert ok, label
