import argparse

from rvgomea.defaults import DEFAULT_LINKAGE_MODEL, DEFAULT_PROBLEM, DEFAULT_DIMENSIONALITY
from rvgomea.experiments.pop_sweep_runner import run_pop_sweep
from rvgomea.run_config import RunConfig


def main():
    parser = argparse.ArgumentParser(
        prog='Single pop sweep')

    parser.add_argument('-l', '--linkage-model', type=str, default=DEFAULT_LINKAGE_MODEL)
    parser.add_argument('-p', '--problem', type=str, default=DEFAULT_PROBLEM)
    parser.add_argument('-d', '--dimensionality', type=int, default=DEFAULT_DIMENSIONALITY)
    parser.add_argument('-i', '--incremental', action=argparse.BooleanOptionalAction)
    parser.add_argument('-g', '--generational-elitist', action=argparse.BooleanOptionalAction)

    args = parser.parse_args()

    base_run_config = RunConfig(
        linkage_model=args.linkage_model,
        population_size=-1,
        random_seed=-1,
        problem=args.problem,
        dimensionality=args.dimensionality,
        lower_init_bound=-115,
        upper_init_bound=-110,
    )
    
    output_dir = f"pop_sweep_test/{args.problem}_{args.dimensionality}_{args.linkage_model}_{'gen_elitist' if args.generational_elitist else 'fos_elitist'}{'_incremental' if args.incremental else ''}"

    run_pop_sweep(
        output_dir, base_run_config, 30, num_cpus=24,
        log_progress=True, incremental=args.incremental, gen_elitist=args.generational_elitist
    )


if __name__ == '__main__':
    main()
