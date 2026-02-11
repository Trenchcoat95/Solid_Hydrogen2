import os
import sys
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

# List notebooks in execution order
NOTEBOOKS = [
    "SH_mu_nu_reco_performance.ipynb",
    "SH_sample_properties.ipynb"
]

# Directory where notebooks live
BASE_DIR = os.path.abspath(".")

# Execution timeout per notebook (seconds)
TIMEOUT = 6000  # 100 minutes

def run_notebook(path):
    print(f"Running: {path}")
    
    with open(path, "r", encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)

    ep = ExecutePreprocessor(timeout=TIMEOUT, kernel_name="python3")

    try:
        ep.preprocess(nb, {"metadata": {"path": os.path.dirname(path) or "./"}})
    except Exception as e:
        print(f"Error executing {path}")
        raise

    # # Save executed notebook (with outputs)
    # with open(path, "w", encoding="utf-8") as f:
    #     nbformat.write(nb, f)

    print(f"Finished: {path}\n")


def main():
    for nb_name in NOTEBOOKS:
        nb_path = os.path.join(BASE_DIR, nb_name)
        if not os.path.exists(nb_path):
            print(f"Notebook not found: {nb_path}")
            sys.exit(1)

        run_notebook(nb_path)

    print("All notebooks executed successfully.")


if __name__ == "__main__":
    main()
