# Advanced System Theory

Study notes, step-by-step solutions and Python checks for the course
**Advanced System Theory** (Paderborn University, Automatic Control Group RAT, Prof. Dr. Erdal Kayacan), winter semester 2026/27.

The notes are written for someone learning system theory from scratch: every concept is explained with a recipe and worked examples, and no solution steps are skipped.

## Repository layout

```
Advanced-System-Theory/
├── README.md                                   ← you are here
├── course-material/                            ← put your local copies of the course PDFs here (not published)
├── Exercise-02_System-Properties-and-Linearization/
    ├── README.md                               ← overview and progress of this exercise
    ├── concepts.md                             ← Lecture 2 theory: everything needed for Exercise 2
    ├── solutions.md                            ← step-by-step solutions
    ├── pdf/                                    ← PDF versions of concepts.md and solutions.md
    └── python/                                 ← scripts that verify every result, and figures
└── tools/
    ├── build_pdfs.sh                           ← regenerate the PDFs from the Markdown
    └── pdf-header.tex                          ← LaTeX header: bold vectors/matrices as in the lecture
```

## Exercises

| # | Lecture | Topic | Concepts | Solutions | Code | Status |
|---|---|---|---|---|---|---|
| 2 | Week 02 | System properties and linearization | [md](Exercise-02_System-Properties-and-Linearization/concepts.md) · [pdf](Exercise-02_System-Properties-and-Linearization/pdf/concepts.pdf) | [md](Exercise-02_System-Properties-and-Linearization/solutions.md) · [pdf](Exercise-02_System-Properties-and-Linearization/pdf/solutions.pdf) | [python](Exercise-02_System-Properties-and-Linearization/python) | Tasks 1–3 ✅ |

## How to use

1. Read the **concepts** file first. It builds up everything needed for the exercise sheet.
2. Try the tasks yourself, then compare with the **solutions** file.
3. Run the **Python** scripts to check the numbers:
   ```bash
   pip install numpy scipy matplotlib sympy
   python Exercise-02_System-Properties-and-Linearization/python/task1_system_properties.py
   python Exercise-02_System-Properties-and-Linearization/python/task2_linearization.py
   python Exercise-02_System-Properties-and-Linearization/python/task3_discrete_time.py
   ```

The Markdown files use LaTeX math, which GitHub renders directly. The PDFs can be rebuilt with `bash tools/build_pdfs.sh` (needs `pandoc` and `xelatex`).
