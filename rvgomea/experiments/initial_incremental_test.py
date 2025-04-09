from IPython.display import display
import pandas as pd
import os

if __name__ == "__main__":
    dir = os.path.dirname(__file__)
    todo_df = pd.read_csv(os.path.join(dir, "todo.csv"))
    test_template_string = "/bin/bash -c \"parallel -j 24 -n0 ./RV-GOMEA -f -100110110 -i -r {} {} -115 -110 0 0.35 {} 1 0.9 1.0 2000000 1e-10 100 0 180 ::: {{1..100}} | awk '{{SUM+=\\$2;}} END {{print SUM/NR, NR}}'\""
    for index, row in todo_df.iterrows():
        command = test_template_string.format(row["problem_idx"], row["dimensionality"], row["population_size"]/2)
        print(command)
        output = os.popen(command).read()
        print(output.split(" "))
        todo_df.at[index, "average_incremental"] = output.split(" ")[0]
        todo_df.at[index, "num_runs_incremental"] = output.split(" ")[1].replace("\n", "")
    todo_df.to_csv(os.path.join(dir, "incremental_runs.csv"), index=False)
    display(todo_df)