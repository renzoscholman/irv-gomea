export PYTHONPATH="$PYTHONPATH:."

# PROBLEMS=("sphere" "rosenbrock" "reb2-chain-weak" "reb2-chain-strong" "reb2-chain-alternating" "reb5-no-overlap" "reb5-small-overlap" "reb5-small-overlap-alternating" "osoreb" "reb5-large-overlap" "reb5-disjoint-pairs" "reb-grid")
#PROBLEMS=("sphere" "rosenbrock" "reb2-chain-weak" "reb2-chain-strong" "reb2-chain-alternating" "reb5-no-overlap" "reb5-small-overlap" "reb5-small-overlap-alternating" "osoreb" "osoreb-big-strong" "osoreb-small-strong" "reb5-large-overlap" "reb5-disjoint-pairs" "reb-grid")
# PROBLEMS=("reb10-no-overlap" "reb10-small-overlap" "reb10-small-overlap-alternating" "reb10-large-overlap")
PROBLEMS=("sphere" "rosenbrock" "reb2-chain-weak" "reb2-chain-strong" "reb2-chain-alternating" "reb5-no-overlap" "reb5-small-overlap" "reb5-small-overlap-alternating" "osoreb" "osoreb-big-strong" "osoreb-small-strong" "reb5-large-overlap" "reb5-disjoint-pairs" "reb-grid" "reb10-no-overlap" "reb10-small-overlap" "reb10-small-overlap-alternating" "reb10-large-overlap")
#FRAGMENTS=("snellius3") #"shark" "snellius1" "snellius2")
#FRAGMENTS=("snellius3" "shark2" "shark3")
FRAGMENTS=("shark3")
BASEDIR="data" # "gecco-data"

# for problem in "${PROBLEMS[@]}"
# do
#     mkdir ${BASEDIR}/scalability-bisection-${problem}

#     for fragment in "${FRAGMENTS[@]}"
#     do
#         rsync -P -a ${BASEDIR}/scalability-all-results/${fragment}/data/scalability-bisection-${problem}/ data/scalability-bisection-${problem}/
#     done
# done

for problem in "${PROBLEMS[@]}"
do
  echo ${problem}
  # python rvgomea/cmd/aggregate_set_of_bisections.py data/scalability-bisection-incremental-${problem}
 #  python rvgomea/cmd/aggregate_set_of_bisections.py /mnt/constellation/export/scratch2/data/rjs/experiments/conditional/final/scalability-bisection-incremental-${problem}
  #python rvgomea/cmd/aggregate_set_of_bisections.py /mnt/constellation/export/scratch2/data/rjs/experiments/conditional/baseline_generational_elitist/scalability-bisection-${problem}
  python rvgomea/cmd/aggregate_set_of_bisections.py /mnt/constellation/export/scratch2/data/rjs/experiments/conditional/final_bisections/data_rvg/scalability-bisection-${problem}
  #python rvgomea/cmd/aggregate_set_of_bisections.py /mnt/constellation/export/scratch2/data/rjs/experiments/conditional/final_bisections/data_dmd_halved_1.0/scalability-bisection-incremental-${problem}

  # python rvgomea/cmd/aggregate_set_of_bisections.py /mnt/constellation/export/scratch2/data/rjs/experiments/conditional/baseline/scalability-bisection-${problem}

  echo ""
done
