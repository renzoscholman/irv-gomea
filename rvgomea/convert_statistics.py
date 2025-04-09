import os

import pandas as pd


def convert_statistics(input_filename: str, output_filename: str):
    rows = []
    with open(input_filename) as f:
        f.readline()
        for line in f.readlines():
            tokens = line.strip().split()
            generations = int(tokens[0])
            evaluations = float(tokens[1])
            seconds = float(tokens[2])
            best_objective = float(tokens[3])
            chol_fails = int(tokens[5])

            rows.append({
                "generations": generations,
                "evaluations": evaluations,
                "seconds": seconds,
                "best_objective": best_objective,
                "cholesky_fails": chol_fails,
            })

    os.system(f"rm {input_filename}")

    df = pd.DataFrame(rows)

    if output_filename is not None:
        df.to_csv(output_filename)

    return df
