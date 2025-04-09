/**
 *
 * Fitness-based Conditional Real-Valued Gene-pool Optimal Mixing Evolutionary Algorithm
 *
 * Copyright (c) 2024 by Georgios Andreadis, Tanja Alderliesten, Peter A.N. Bosman, Anton Bouter, and Chantal Olieman
 * This code is licensed under CC BY-NC-ND 4.0. A copy of the license is included in the LICENSE file.
 *
 * If you use this software for any purpose, please cite the most recent pre-print titled:
 * "Fitness-based Linkage Learning and Maximum-Clique Conditional Linkage Modelling for Gray-box Optimization
 *  with RV-GOMEA", by Georgios Andreadis, Tanja Alderliesten, and Peter A.N. Bosman. 2024.
 *
 * IN NO EVENT WILL THE AUTHOR OF THIS SOFTWARE BE LIABLE TO YOU FOR ANY
 * DAMAGES, INCLUDING BUT NOT LIMITED TO LOST PROFITS, LOST SAVINGS, OR OTHER
 * INCIDENTIAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR THE INABILITY
 * TO USE SUCH PROGRAM, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGES, OR FOR ANY CLAIM BY ANY OTHER PARTY. THE AUTHOR MAKES NO
 * REPRESENTATIONS OR WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE
 * AUTHOR SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY ANYONE AS A RESULT OF
 * USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES.
 *
 */

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "rv-gomea.h"

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

std::vector<population_t *> populations;

/*-=-=-=-=-=-=-=-=-=-=- Section Interpret Command Line -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Parses and checks the command line.
 */
rvg_t::rvg_t(int argc, char **argv) {
    startTimer();

    haveNextNextGaussian = 0;

    parseCommandLine(argc, argv);

    if (use_guidelines) {
        tau = 0.35;
        // if (maximum_number_of_populations == 1)
        //     base_population_size = (int) (36.1 + 7.58 * log2((double) number_of_parameters));
        // else
        //     base_population_size = 10;
        //base_population_size           = (int) (10.0*pow((double) number_of_parameters,0.5));
        base_population_size           = (int) (17.0 + 3.0*pow((double) number_of_parameters,1.5));
        if (use_incremental)
            base_population_size = (int) 10 + 3 * number_of_parameters;
        //base_population_size           = (int) (4.0*pow((double) number_of_parameters,0.5));
        distribution_multiplier_decrease = use_incremental ? 0.95 : 0.9;
        st_dev_ratio_threshold = 1.0;
        maximum_no_improvement_stretch = 25 + number_of_parameters;
    }

    checkOptions();
}

/**
 * Parses the command line.
 * For options, see printUsage.
 */
void rvg_t::parseCommandLine(int argc, char **argv) {
    int index;

    index = 1;

    parseOptions(argc, argv, &index);

    parseParameters(argc, argv, &index);
}

/**
 * Parses only the options from the command line.
 */
void rvg_t::parseOptions(int argc, char **argv, int *index) {
    double dummy;

    write_generational_statistics = 0;
    write_generational_solutions = 0;
    print_verbose_overview = 0;
    use_vtr = 0;
    use_guidelines = 0;
    black_box_evaluations = 0;

    for (; (*index) < argc; (*index)++) {
        if (argv[*index][0] == '-') {
            /* If it is a negative number, the option part is over */
            if (sscanf(argv[*index], "%lf", &dummy) && argv[*index][1] != '\0')
                break;

            if (argv[*index][1] == '\0')
                optionError(argv, *index);
            else if (argv[*index][2] != '\0')
                optionError(argv, *index);
            else {
                switch (argv[*index][1]) {
                    case 'h':
                        printUsage();
                        break;
                    case 'P':
                        printAllInstalledProblems();
                        break;
                    case 's':
                        write_generational_statistics = 1;
                        break;
                    case 'w':
                        write_generational_solutions = 1;
                        break;
                    case 'd':
                        write_fitness_dependencies = 1;
                        break;
                    case 'v':
                        print_verbose_overview = 1;
                        break;
                    case 'r':
                        use_vtr = 1;
                        break;
                    case 'i':
                        use_incremental = 1;
                        break;
                    case 'g':
                        use_guidelines = 1;
                        break;
                    case 'b':
                        black_box_evaluations = 1;
                        break;
                    case 'f':
                        parseFOSElementSize(index, argc, argv);
                        break;
                    case 'S':
                        parseSeed(index, argc, argv);
                        break;
                    default :
                        optionError(argv, *index);
                }
            }
        } else /* Argument is not an option, so option part is over */
            break;
    }
}

void rvg_t::parseFOSElementSize(int *index, int argc, char **argv) {
    short noError = 1;

    (*index)++;
    noError = noError && sscanf(argv[*index], "%d", &FOS_element_size);

    if (!noError) {
        printf("Error parsing parameters.\n\n");

        printUsage();
    }
}

void rvg_t::parseSeed(int *index, int argc, char **argv) {
    short noError = 1;

    (*index)++;
    noError = noError && sscanf(argv[*index], "%ld", &random_seed);

    if (!noError) {
        printf("Error parsing parameters.\n\n");

        printUsage();
    }
}

/**
 * Checks whether the selected options are feasible.
 */
void rvg_t::checkOptions(void) {
    if (number_of_parameters < 1) {
        printf("\n");
        printf("Error: number of parameters < 1 (read: %d). Require number of parameters >= 1.", number_of_parameters);
        printf("\n\n");

        exit(0);
    }

    if (((int) (tau * base_population_size)) <= 0 || tau >= 1) {
        printf("\n");
        printf("Error: tau not in range (read: %e). Require tau in [1/pop,1] (read: [%e,%e]).", tau,
               1.0 / ((double) base_population_size), 1.0);
        printf("\n\n");

        exit(0);
    }

    if (base_population_size < 1) {
        printf("\n");
        printf("Error: population size < 1 (read: %d). Require population size >= 1.", base_population_size);
        printf("\n\n");

        exit(0);
    }

    if (maximum_number_of_populations < 1) {
        printf("\n");
        printf("Error: number of populations < 1 (read: %d). Require number of populations >= 1.",
               maximum_number_of_populations);
        printf("\n\n");

        exit(0);
    }

    if (installedProblemName(problem_index) == NULL) {
        printf("\n");
        printf("Error: unknown index for problem (read index %d).", problem_index);
        printf("\n\n");

        exit(0);
    }

    /*if( rotation_angle > 0 && ( !learn_linkage_tree && FOS_element_size > 1 && FOS_element_size != block_size && FOS_element_size != number_of_parameters) )
    {
        printf("\n");
        printf("Error: invalid FOS element size (read %d). Must be %d, %d or %d.", FOS_element_size, 1, block_size, number_of_parameters );
        printf("\n\n");

        exit( 0 );
    }*/
}


/**
 * Informs the user of an illegal option and exits the program.
 */
void rvg_t::optionError(char **argv, int index) {
    printf("Illegal option: %s\n\n", argv[index]);

    printUsage();
}

/**
 * Parses only the EA parameters from the command line.
 */
void rvg_t::parseParameters(int argc, char **argv, int *index) {
    if ((argc - *index) < 15) {
        printf("Number of parameters is incorrect, require 15 parameters (you provided %d).\n\n", (argc - *index));

        printUsage();
    }

    // use_conditional_sampling = false;
    selection_during_gom = 0;
    alter_full_fos_gom_acceptance = 0;
    update_elitist_during_gom = 1;
    perform_factorized_gom = 1;
    perform_eda_gom = 0;

    int a, b, c;

    int noError = 1;
    noError = noError && sscanf(argv[*index + 0], "%d", &problem_index);
    noError = noError && sscanf(argv[*index + 1], "%d", &number_of_parameters);
    noError = noError && sscanf(argv[*index + 2], "%lf", &lower_user_range);
    noError = noError && sscanf(argv[*index + 3], "%lf", &upper_user_range);
    noError = noError && sscanf(argv[*index + 4], "%lf", &rotation_angle);
    noError = noError && sscanf(argv[*index + 5], "%lf", &tau);
    noError = noError && sscanf(argv[*index + 6], "%d", &base_population_size);
    noError = noError && sscanf(argv[*index + 7], "%d", &maximum_number_of_populations);
    noError = noError && sscanf(argv[*index + 8], "%lf", &distribution_multiplier_decrease);
    noError = noError && sscanf(argv[*index + 9], "%lf", &st_dev_ratio_threshold);
    noError = noError && sscanf(argv[*index + 10], "%lf", &maximum_number_of_evaluations);
    noError = noError && sscanf(argv[*index + 11], "%lf", &vtr);
    noError = noError && sscanf(argv[*index + 12], "%d", &maximum_no_improvement_stretch);
    noError = noError && sscanf(argv[*index + 13], "%lf", &fitness_variance_tolerance);
    noError = noError && sscanf(argv[*index + 14], "%lf", &maximum_number_of_seconds);
    if (argc - *index > 15) {
        noError = noError && sscanf(argv[*index + 15], "%d", &a);
        noError = noError && sscanf(argv[*index + 16], "%d", &b);
        noError = noError && sscanf(argv[*index + 17], "%d", &c);
        selection_during_gom = (short) a;
        update_elitist_during_gom = (short) b;
        alter_full_fos_gom_acceptance = (short) c;
    }

    if (!noError) {
        printf("Error parsing parameters.\n\n");

        printUsage();
    }
}

/**
 * Prints the settings as read from the command line.
 */
void rvg_t::printVerboseOverview(void) {
    int i;

    printf("### Settings ######################################\n");
    printf("#\n");
    printf("# Statistics writing every generation: %s\n", write_generational_statistics ? "enabled" : "disabled");
    printf("# Population file writing            : %s\n", write_generational_solutions ? "enabled" : "disabled");
    printf("# Use of value-to-reach (vtr)        : %s\n", use_vtr ? "enabled" : "disabled");
    printf("#\n");
    printf("###################################################\n");
    printf("#\n");
    printf("# Problem                 = %s\n", installedProblemName(problem_index));
    printf("# Number of parameters    = %d\n", number_of_parameters);
    printf("# Initialization ranges   = [%e;%e]\n", lower_user_range, upper_user_range);
    printf("# Boundary ranges         = ");
    for (i = 0; i < number_of_parameters; i++) {
        if (i > 0){
            if(fitness->getLowerRangeBound(i) == fitness->getLowerRangeBound(i-1) && 
               fitness->getUpperRangeBound(i) == fitness->getUpperRangeBound(i-1))
               continue;
        }
        printf("x_%d: [%e;%e]", i, fitness->getLowerRangeBound(i), fitness->getUpperRangeBound(i));
        if (i < number_of_parameters - 1)
            printf("\n#                           ");
    }
    printf("\n");
    printf("# Rotation angle            = %e\n", rotation_angle);
    printf("# Tau                       = %e\n", tau);
    printf("# Population size/normal    = %d\n", base_population_size);
    printf("# Selection during gom      = %d\n",selection_during_gom);
    printf("# Update elitist during gom = %d\n",update_elitist_during_gom);
    printf("# FOS element size          = %d\n", FOS_element_size);
    printf("# FOS element UB            = %d\n", FOS_element_ub);
    printf("# FOS max clique size       = %d\n", max_clique_size);
    printf("# FOS seed cliques per var  = %d\n", seed_cliques_per_variable);
    printf("# FOS use set cover         = %d\n", use_set_cover);
    printf("# FOS use generational GOM  = %d\n", perform_eda_gom);
    printf("# FOS use factorized GOM    = %d\n", perform_factorized_gom);
    printf("# FOS use condit sampling   = %d\n", use_conditional_sampling);
    printf("# FOS include full element  = %d\n", include_full_fos_element);
    printf("# FOS learn conditional LT  = %d\n", learn_conditional_linkage_tree);
    printf("# FOS static LT             = %d\n", static_linkage_tree);
    printf("# FOS learn LT              = %d\n", learn_linkage_tree);
    printf("# FOS prune LT              = %d\n", prune_linkage_tree);
    printf("# FOS similarity measure    = %c\n", similarity_measure);
    printf("# Incremental covariance    = %d\n", use_incremental);
    printf("# Max num of populations    = %d\n", maximum_number_of_populations);
    printf("# Dis. mult. decreaser      = %e\n", distribution_multiplier_decrease);
    printf("# St. dev. rat. threshold   = %e\n", st_dev_ratio_threshold);
    printf("# Maximum numb. of eval.    = %lf\n", maximum_number_of_evaluations);
    printf("# Value to reach (vtr)      = %e\n", vtr);
    printf("# Max. no improv. stretch   = %d\n", maximum_no_improvement_stretch);
    printf("# Fitness var. tolerance    = %e\n", fitness_variance_tolerance);
    printf("# Random seed               = %ld\n", random_seed);
    printf("#\n");
    printf("###################################################\n");
}

/**
 * Prints usage information and exits the program.
 */
void rvg_t::printUsage(void) {
    printf("Usage: RV-GOMEA [-?] pro dim low upp rot tau pop nop dmd srt eva vtr imp tol\n");
    printf(" -h: Prints out this usage information.\n");
    printf(" -P: Prints out a list of all installed optimization problems.\n");
    printf(" -s: Enables computing and writing of statistics every generation.\n");
    printf(" -w: Enables writing of solutions and their fitnesses every generation.\n");
    printf(" -v: Enables verbose mode. Prints the settings before starting the run.\n");
    printf(" -r: Enables use of vtr in termination condition (value-to-reach).\n");
    printf(" -b: Enables counting every partial evaluation as a full evaluation.\n");
    printf(" -f %%d: Sets linkage model that is used. Positive: Use a FOS with elements of %%d consecutive variables.\n");
    printf("     Use -1 for full linkage model, -2 for dynamic linkage tree learned from the population, -3 for fixed linkage tree learned from distance measure,\n");
    printf("     -4 for bounded fixed linkage tree learned from distance measure, -5 for fixed bounded linkage tree learned from random distance measure.\n");
    printf(" -g: Uses guidelines to override parameter settings for those parameters\n");
    printf("     for which a guideline is known in literature. These parameters are:\n");
    printf("     tau pop dmd srt imp\n");;
    printf(" -S: A fixed random seed is used.\n");

    printf("\n");
    printf("  pro: Index of optimization problem to be solved (minimization).\n");
    printf("  dim: Number of parameters.\n");
    printf("  low: Overall initialization lower bound.\n");
    printf("  upp: Overall initialization upper bound.\n");
    printf("  rot: The angle by which to rotate the problem.\n");
    printf("  tau: Selection percentile (tau in [1/pop,1], truncation selection).\n");
    printf("  pop: Population size per normal.\n");
    printf("  nop: The number of populations (parallel runs that initially partition the search space).\n");
    printf("  dmd: The distribution multiplier decreaser (in (0,1), increaser is always 1/dmd).\n");
    printf("  srt: The standard-devation ratio threshold for triggering variance-scaling.\n");
    printf("  eva: Maximum number of evaluations allowed.\n");
    printf("  vtr: The value to reach. If the objective value of the best feasible solution reaches\n");
    printf("       this value, termination is enforced (if -r is specified).\n");
    printf("  imp: Maximum number of subsequent generations without an improvement while the\n");
    printf("       the distribution multiplier is <= 1.0.\n");
    printf("  tol: The tolerance level for fitness variance (i.e. minimum fitness variance)\n");
    printf("  sec: The time limit in seconds.\n");
    exit(0);
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Initialize -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Performs initializations that are required before starting a run.
 */
void rvg_t::initialize(void) {
    total_number_of_writes = 0;
    number_of_subgenerations_per_population_factor = 8;//use_incremental ? 12 : 8;

    initializeRandomNumberGenerator();

    initializeProblem(problem_index, vtr);
}

void rvg_t::restartLargestPopulation() {
    int pop_size = populations[populations.size() - 1]->population_size;
    population_t *new_pop = new population_t(fitness, pop_size, lower_user_range, upper_user_range);
    new_pop->maximum_no_improvement_stretch = maximum_no_improvement_stretch;
    new_pop->st_dev_ratio_threshold = st_dev_ratio_threshold;
    new_pop->distribution_multiplier_decrease = distribution_multiplier_decrease;
    new_pop->maximum_no_improvement_stretch = maximum_no_improvement_stretch;
    new_pop->tau = tau;
    delete (populations[populations.size() - 1]);
    populations[populations.size() - 1] = new_pop;
}


void rvg_t::initializeNewPopulation() {
    int new_pop_size = base_population_size;
    if (populations.size() > 0) new_pop_size = 2 * populations[populations.size() - 1]->population_size;
    population_t *new_pop = new population_t(fitness, new_pop_size, lower_user_range, upper_user_range);
    new_pop->maximum_no_improvement_stretch = maximum_no_improvement_stretch;
    new_pop->st_dev_ratio_threshold = st_dev_ratio_threshold;
    new_pop->distribution_multiplier_decrease = distribution_multiplier_decrease;
    new_pop->maximum_no_improvement_stretch = maximum_no_improvement_stretch;
    new_pop->tau = tau;
    populations.push_back(new_pop);

    // Delay to here so all FOS info has also been set for correct verbose output
    if(populations.size() == 1) {
        if (print_verbose_overview)
            printVerboseOverview();
    }
}

void rvg_t::initializeProblem(int problem_index, double vtr) {
    fitness = fitness_t::getFitnessClass(problem_index, number_of_parameters, vtr);
    if (fitness == NULL) {
        printf("Unknown problem index.\n");
        fflush(stdout);
        exit(1);
    }

    // Check if variable interaction graph does not include self-loops
    for (int i = 0; i < number_of_parameters; i++) {
        std::set<int> dependent_set = fitness->variable_interaction_graph[i];
        if (std::find(dependent_set.begin(), dependent_set.end(), i) != dependent_set.end()) {
            printf("Variable interaction graph includes itself.\n");
            fflush(stdout);
            exit(2);
        }
    }

    fitness->use_vtr = use_vtr;
    fitness->black_box_optimization = black_box_evaluations;
    fitness->vtr_hit_status = 0;
    fitness->elitist_objective_value = 1e308;
    fitness->elitist_constraint_value = 1e308;
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/




/*-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Output =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Writes (appends) statistics about the current generation to a
 * file named "statistics.dat".
 */
void rvg_t::writeGenerationalStatisticsForOnePopulation(int population_index) {
    char string[1000];
    FILE *file;

    /* Average, best and worst */
    double population_objective_avg = 0.0;
    double population_constraint_avg = 0.0;
    double population_objective_best = populations[population_index]->individuals[0]->objective_value;
    double population_constraint_best = populations[population_index]->individuals[0]->constraint_value;
    double population_objective_worst = populations[population_index]->individuals[0]->objective_value;
    double population_constraint_worst = populations[population_index]->individuals[0]->constraint_value;
    solution_t * best_ind = populations[population_index]->individuals[0];
    for (int j = 0; j < populations[population_index]->population_size; j++) {
        best_ind = fitness->betterFitness(populations[population_index]->individuals[j], best_ind) ? populations[population_index]->individuals[j] : best_ind;
        population_objective_avg += populations[population_index]->individuals[j]->objective_value;
        population_constraint_avg += populations[population_index]->individuals[j]->constraint_value;
        if (fitness_t::betterFitness(population_objective_worst, population_constraint_worst,
                                     populations[population_index]->individuals[j]->objective_value,
                                     populations[population_index]->individuals[j]->constraint_value)) {
            population_objective_worst = populations[population_index]->individuals[j]->objective_value;
            population_constraint_worst = populations[population_index]->individuals[j]->constraint_value;
        }
        if (fitness_t::betterFitness(populations[population_index]->individuals[j]->objective_value,
                                     populations[population_index]->individuals[j]->constraint_value,
                                     population_objective_best, population_constraint_best)) {
            population_objective_best = populations[population_index]->individuals[j]->objective_value;
            population_constraint_best = populations[population_index]->individuals[j]->constraint_value;
        }
    }
    population_objective_avg = population_objective_avg / ((double) populations[population_index]->population_size);
    population_constraint_avg = population_constraint_avg / ((double) populations[population_index]->population_size);

    /* Variance */
    double population_objective_var = 0.0;
    double population_constraint_var = 0.0;
    for (int j = 0; j < populations[population_index]->population_size; j++) {
        population_objective_var +=
                (populations[population_index]->individuals[j]->objective_value - population_objective_avg) *
                (populations[population_index]->individuals[j]->objective_value - population_objective_avg);
        population_constraint_var +=
                (populations[population_index]->individuals[j]->constraint_value - population_constraint_avg) *
                (populations[population_index]->individuals[j]->constraint_value - population_constraint_avg);
    }
    population_objective_var = population_objective_var / ((double) populations[population_index]->population_size);
    population_constraint_var = population_constraint_var / ((double) populations[population_index]->population_size);

    if (population_objective_var <= 0.0)
        population_objective_var = 0.0;
    if (population_constraint_var <= 0.0)
        population_constraint_var = 0.0;

    /* Then write them */
    file = NULL;
    if (total_number_of_writes == 0) {
        file = fopen("statistics.dat", "w");

        sprintf(string,
                "# Generation  Evaluations  Time(s)  Best-obj. Best-cons. Chol-Fails [Pop.index  Subgen.  Pop.size  Dis.mult.[0]  Pop.best.obj. Pop.avg.obj.  Pop.var.obj. Pop.worst.obj.  Pop.best.con. Pop.avg.con.  Pop.var.con. Pop.worst.con.]\n");
        fputs(string, file);
    } else
        file = fopen("statistics.dat", "a");

    sprintf(string, "%10d %11lf %11.3lf %20.15e %13e %10d ", total_number_of_generations, fitness->number_of_evaluations, getTimer(), fitness->elitist_objective_value, fitness->elitist_constraint_value, cholesky_fails);
    fputs(string, file);

    double min_dist_mult, mean_dist_mult, max_dist_mult;
    double tot_dist_mult = 0;
    for(int i = 0; i < populations.size(); i++){
        min_dist_mult = populations[i]->linkage_model->distributions[0]->distribution_multiplier;
        tot_dist_mult += populations[i]->linkage_model->distributions[0]->distribution_multiplier;
        max_dist_mult = populations[i]->linkage_model->distributions[0]->distribution_multiplier;
        for(int j = 1; j < populations[i]->linkage_model->distributions.size(); j++){
            min_dist_mult = min(min_dist_mult, populations[i]->linkage_model->distributions[j]->distribution_multiplier);
            tot_dist_mult += populations[i]->linkage_model->distributions[j]->distribution_multiplier;
            max_dist_mult = max(max_dist_mult, populations[i]->linkage_model->distributions[j]->distribution_multiplier);
        }
        mean_dist_mult = tot_dist_mult / populations[i]->linkage_model->distributions.size();
    }

    sprintf(string, "%4.8e %4.8e %4.8e \t\t\t", min_dist_mult, mean_dist_mult, max_dist_mult);
    fputs(string, file);

    double min_sdr, mean_sdr, max_sdr;
    double tot_sdr = 0;
    for(int i = 0; i < populations.size(); i++){
        min_sdr = populations[i]->linkage_model->distributions[0]->st_dev_ratio;
        tot_sdr += populations[i]->linkage_model->distributions[0]->st_dev_ratio;
        max_sdr = populations[i]->linkage_model->distributions[0]->st_dev_ratio;
        for(int j = 1; j < populations[i]->linkage_model->distributions.size(); j++){
            min_sdr = min(min_sdr, populations[i]->linkage_model->distributions[j]->st_dev_ratio);
            tot_sdr += populations[i]->linkage_model->distributions[j]->st_dev_ratio;
            max_sdr = max(max_sdr, populations[i]->linkage_model->distributions[j]->st_dev_ratio);
        }
        mean_sdr = tot_sdr / populations[i]->linkage_model->distributions.size();
    }

    sprintf(string, "%4.8e %4.8e %4.8e \t\t\t", min_sdr, mean_sdr, max_sdr);
    fputs(string, file);

    total_number_of_generations++;

    // for(int i = 0; i < fitness->number_of_parameters; i++){
    //     sprintf(string, "% 3.8e ", best_ind->variables[i]);
    //     fputs(string, file);
    // }

    //sprintf( string, "[ %4d %6d %10d %13e %13e %13e %13e %13e %13e %13e %13e %13e ]", population_index, number_of_generations[population_index], population_sizes[population_index], distribution_multipliers[population_index][0], population_objective_best, population_objective_avg, population_objective_var, population_objective_worst, population_constraint_best, population_constraint_avg, population_constraint_var, population_constraint_worst );
    //fputs( string, file );

    sprintf(string, "\n");
    fputs(string, file);

    fclose(file);

    total_number_of_writes++;
}

/**
 * Writes the solutions to various files. The filenames
 * contain the generation. If the flag final is set
 * (final != 0), the generation number in the filename
 * is replaced with the word "final".
 *
 * all_populations_generation_xxxxx.dat : all populations combined
 * population_xxxxx_generation_xxxxx.dat: the individual populations
 * selection_xxxxx_generation_xxxxx.dat : the individual selections
 */
void rvg_t::writeGenerationalSolutions(short final) {
    int i, j, k;
    char string[1000];
    FILE *file_all, *file_population, *file_selection;

    file_selection = NULL;
    if (final)
        sprintf(string, "all_populations_generation_final.dat");
    else
        sprintf(string, "all_populations_generation_%05ld.dat", populations.size());
    file_all = fopen(string, "w");

    for (i = 0; i < populations.size(); i++) {
        if (final)
            sprintf(string, "population_%05d_generation_final.dat", i);
        else
            sprintf(string, "population_%05d_generation_%05d.dat", i, populations[i]->number_of_generations);
        file_population = fopen(string, "w");

        if( populations[i]->number_of_generations > 0 && !final )
        {
            populations[i]->makeSelection();
            sprintf(string, "selection_%05d_generation_%05d.dat", i, populations[i]->number_of_generations);
            file_selection = fopen(string, "w");
        }

        /* Populations */
        for (j = 0; j < populations[i]->population_size; j++) {
            for (k = 0; k < number_of_parameters; k++) {
                sprintf(string, "%13e", populations[i]->individuals[j]->variables[k]);
                fputs(string, file_all);
                fputs(string, file_population);
                if (k < number_of_parameters - 1) {
                    sprintf(string, " ");
                    fputs(string, file_all);
                    fputs(string, file_population);
                }
            }
            sprintf(string, "     ");
            fputs(string, file_all);
            fputs(string, file_population);
            sprintf(string, "%13e %13e", populations[i]->individuals[j]->objective_value,
                    populations[i]->individuals[j]->constraint_value);
            fputs(string, file_all);
            fputs(string, file_population);
            sprintf(string, "\n");
            fputs(string, file_all);
            fputs(string, file_population);
        }

        fclose(file_population);

        /* Selections */
        if( populations[i]->number_of_generations > 0 && !final )
        {
            for( j = 0; j < populations[i]->selection_size; j++ )
            {
                for( k = 0; k < number_of_parameters; k++ )
                {
                    sprintf( string, "%13e", populations[i]->selection[j]->variables[k] );
                    fputs( string, file_selection );
                    if( k < number_of_parameters-1 )
                    {
                        sprintf( string, " " );
                        fputs( string, file_selection );
                    }
                    sprintf( string, "     " );
                    fputs( string, file_selection );
                }
                sprintf( string, "%13e %13e", populations[i]->selection[j]->objective_value, populations[i]->selection[j]->constraint_value );
                fputs( string, file_selection );
                sprintf( string, "\n" );
                fputs( string, file_selection );
            }
            fclose( file_selection );
        }

          
        if( populations[i]->number_of_generations > 0 && !final )
        {
            sprintf(string, "ams_%05d_generation_%05d.dat", i, populations[i]->number_of_generations);
            FILE * file_ams = fopen(string, "w");
            for( k = 0; k < number_of_parameters; k++ )
            {
                sprintf(string, "%13e ", populations[i]->mean_shift_vector[k]);
                fputs( string, file_ams );
            }
            fclose(file_ams);

            
            for(int l = 0; l < populations[i]->linkage_model->distributions.size(); l++){
                sprintf(string, "cov_%05d_fos_%05d_generation_%05d.dat", i, l, populations[i]->number_of_generations);
                FILE * file_cov = fopen(string, "w");
                for( k = 0; k < populations[i]->linkage_model->distributions[l]->variables.size(); k++ )
                {
                    for( j = 0; j < populations[i]->linkage_model->distributions[l]->variables.size(); j++ )
                    {
                        sprintf(string, "%13e ", ((conditional_distribution_t*) populations[i]->linkage_model->distributions[l])->full_covariance_matrices[0](k,j));
                        fputs( string, file_cov );
                    }
                    sprintf(string, "\n");
                    fputs( string, file_cov );
                }
                fclose(file_cov);
            }
        }

        // fos_t * lm = populations[i]->linkage_model;
        // for(int j = 0; j < lm->getLength(); j++){
        //     auto dist = (conditional_distribution_t*) lm->distributions[j];
        //     for(int k = 0; k < dist->variable_groups.size(); k++){
        //         auto c = dist->variable_groups[k];
        //         int idx = -1;
        //         if (std::find(c.begin(), c.end(), 0) != c.end() && std::find(c.begin(), c.end(), 1) != c.end()) idx = 0;
        //         if (std::find(c.begin(), c.end(), 1) != c.end() && std::find(c.begin(), c.end(), 2) != c.end()) idx = 1;
        //         if (std::find(c.begin(), c.end(), 2) != c.end() && std::find(c.begin(), c.end(), 3) != c.end()) idx = 3;
        //         if(idx < 0) continue;
                
        //         sprintf(string, "ams_%05d_generation_%05d_%05d.dat", i, populations[i]->number_of_generations, idx);
        //         FILE * file_ams = fopen(string, "w");
        //         for(double val : dist->mean_shift_vectors[k]){
        //             sprintf(string, "%13e ", val);
        //         }
        //         fclose(file_ams);
        //     }
        // }
    }

    fclose(file_all);

    writeGenerationalSolutionsBest(final);
}


/**
 * Writes the best solution (measured in the single
 * available objective) to a file named
 * best_generation_xxxxx.dat where xxxxx is the
 * generation number. If the flag final is set
 * (final != 0), the generation number in the filename
 * is replaced with the word "final".The output
 * file contains the solution values with the
 * dimensions separated by a single white space,
 * followed by five white spaces and then the
 * single objective value for that solution
 * and its sum of constraint violations.
 */
void rvg_t::writeGenerationalSolutionsBest(short final) {
    int i, population_index_best, individual_index_best;
    char string[1000];
    FILE *file;
    static int *c = NULL;
    if (c == NULL) {
        c = (int *) Malloc(sizeof(int));
        c[0] = 0;
    }

    /* First find the best of all */
    determineBestSolutionInCurrentPopulations(&population_index_best, &individual_index_best);

    /* Then output it */
    if (final)
        sprintf(string, "best_generation_final.dat");
    else
        sprintf(string, "best_generation_%05d.dat", c[0]);
    file = fopen(string, "w");

    for (i = 0; i < number_of_parameters; i++) {
        sprintf(string, "%13e", populations[population_index_best]->individuals[individual_index_best]->variables[i]);
        fputs(string, file);
        if (i < number_of_parameters - 1) {
            sprintf(string, " ");
            fputs(string, file);
        }
    }
    sprintf(string, "     ");
    fputs(string, file);
    sprintf(string, "%13e %13e",
            populations[population_index_best]->individuals[individual_index_best]->objective_value,
            populations[population_index_best]->individuals[individual_index_best]->constraint_value);
    fputs(string, file);
    sprintf(string, "\n");
    fputs(string, file);
    c[0]++;

    fclose(file);
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Termination -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns 1 if termination should be enforced, 0 otherwise.
 */
short rvg_t::checkTerminationCondition(void) {
    short allTrue;
    int i;

    if (checkNumberOfEvaluationsTerminationCondition())
        return (1);

    if (checkVTRTerminationCondition())
        return (1);

    if (checkTimeLimitTerminationCondition())
        return (1);

    checkAverageFitnessTerminationConditions();

    if (populations.size() < maximum_number_of_populations)
        return (0);

    allTrue = 1;
    for (i = 0; i < populations.size(); i++) {
        if (!populations[i]->population_terminated) {
            allTrue = 0;
            break;
        }
    }
    if (allTrue) {
        restartLargestPopulation();
        allTrue = 0;
    }

    return (allTrue);
}

short rvg_t::checkPopulationTerminationConditions(int population_index) {
    /*if( (FOS_element_size == -3 || FOS_element_size == number_of_parameters || FOS_element_size <= -100000 ) && cholesky_fails > 0 )
    {
		//printf("cholesky failed\n");
	    cholesky_fails = 0;
		return( 1 );
    }*/

    if (checkFitnessVarianceTermination(population_index)) {
        printf("Variance-based restart.\n");
        fflush(stdout);
        return (1);
    }

    if (checkDistributionMultiplierTerminationCondition(population_index)) {
        printf("Distribution-based restart.\n");
        fflush(stdout);
        return (1);
    }

    return (0);
}

short rvg_t::checkSubgenerationTerminationConditions() {
    if (checkNumberOfEvaluationsTerminationCondition())
        return (1);

    if (checkVTRTerminationCondition())
        return (1);

    if (checkTimeLimitTerminationCondition())
        return (1);

    return (0);
}

short rvg_t::checkTimeLimitTerminationCondition(void) {
    return (maximum_number_of_seconds > 0 && getTimer() > maximum_number_of_seconds);
}

/**
 * Returns 1 if the maximum number of evaluations
 * has been reached, 0 otherwise.
 */
short rvg_t::checkNumberOfEvaluationsTerminationCondition(void) {
    if (fitness->number_of_evaluations >= maximum_number_of_evaluations && maximum_number_of_evaluations > 0)
        return (1);

    return (0);
}

/**
 * Returns 1 if the value-to-reach has been reached (in any population).
 */
short rvg_t::checkVTRTerminationCondition(void) {
    return (use_vtr && fitness->vtr_hit_status);
}

void rvg_t::checkAverageFitnessTerminationConditions(void) {
    double *average_objective_values = (double *) Malloc(populations.size() * sizeof(double));
    double *average_constraint_values = (double *) Malloc(populations.size() * sizeof(double));
    for (int i = populations.size() - 1; i >= 0; i--) {
        average_objective_values[i] = 0;
        average_constraint_values[i] = 0;
        for (int j = 0; j < populations[i]->population_size; j++) {
            average_objective_values[i] += populations[i]->individuals[j]->objective_value;
            average_constraint_values[i] += populations[i]->individuals[j]->constraint_value;
        }
        average_objective_values[i] /= populations[i]->population_size;
        average_constraint_values[i] /= populations[i]->population_size;
        if (i < populations.size() - 1 &&
            fitness_t::betterFitness(average_objective_values[i + 1], average_constraint_values[i + 1],
                                     average_objective_values[i], average_constraint_values[i])) {
            for (int j = i; j >= 0; j--)
                populations[j]->population_terminated = 1;
            break;
        }
    }
    free(average_objective_values);
    free(average_constraint_values);
}

/**
 * Determines which solution is the best of all solutions
 * in all current populations.
 */
void rvg_t::determineBestSolutionInCurrentPopulations(int *population_of_best, int *index_of_best) {
    (*population_of_best) = 0;
    (*index_of_best) = 0;
    for (int i = 0; i < populations.size(); i++) {
        for (int j = 0; j < populations[i]->population_size; j++) {
            if (fitness_t::betterFitness(populations[i]->individuals[j]->objective_value,
                                         populations[i]->individuals[j]->constraint_value,
                                         populations[(*population_of_best)]->individuals[(*index_of_best)]->objective_value,
                                         populations[(*population_of_best)]->individuals[(*index_of_best)]->constraint_value)) {
                (*population_of_best) = i;
                (*index_of_best) = j;
            }
        }
    }
}

/**
 * Checks whether the fitness variance in any population
 * has become too small (user-defined tolerance).
 */
short rvg_t::checkFitnessVarianceTermination(int population_index) {

    if (populations[population_index]->getFitnessVariance() <
        fitness_variance_tolerance * populations[population_index]->getFitnessMean())
        return (1);
    return (0);
}


/**
 * Checks whether the distribution multiplier in any population
 * has become too small (1e-10).
 */
short rvg_t::checkDistributionMultiplierTerminationCondition(int population_index) {
    int i = population_index;
    if (!populations[i]->population_terminated) {
        short converged = 1;
        for (int j = 0; j < populations[i]->linkage_model->getLength(); j++) {
            if (populations[i]->linkage_model->getDistributionMultiplier(j) > 1e-10) {
                converged = 0;
                break;
            }
        }

        if (converged)
            return (1);
    }
    return (0);
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Run -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
rvg_t::~rvg_t() {
    for (int i = 0; i < populations.size(); i++)
        delete (populations[i]);
    delete (fitness);
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Run -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

void rvg_t::generationalStepAllPopulations() {
    int population_index_smallest, population_index_biggest;

    population_index_biggest = populations.size() - 1;
    population_index_smallest = 0;
    while (population_index_smallest <= population_index_biggest) {
        if (!populations[population_index_smallest]->population_terminated)
            break;

        population_index_smallest++;
    }

    generationalStepAllPopulationsRecursiveFold(population_index_smallest, population_index_biggest);
}

void rvg_t::generationalStepAllPopulationsRecursiveFold(int population_index_smallest, int population_index_biggest) {
    for (int i = 0; i < number_of_subgenerations_per_population_factor - 1; i++) {
        for (int population_index = population_index_smallest;
             population_index <= population_index_biggest; population_index++) {

            if (!populations[population_index]->population_terminated) {
                populations[population_index]->runGeneration();

                if (populations.size() == 1 && write_generational_statistics)
                    writeGenerationalStatisticsForOnePopulation(0);

                if (populations.size() == 1 && write_generational_solutions)
                    writeGenerationalSolutions(0);

                if (checkSubgenerationTerminationConditions()) {
                    for (int j = 0; j < populations.size(); j++)
                        populations[j]->population_terminated = 1;
                    return;
                }

                if (checkPopulationTerminationConditions(population_index)) {
                    populations[population_index]->population_terminated = 1;
                }
            }
        }

        for (int population_index = population_index_smallest;
             population_index < population_index_biggest; population_index++)
            generationalStepAllPopulationsRecursiveFold(population_index_smallest, population_index);
    }
}

void rvg_t::runAllPopulations() {
    while (!checkTerminationCondition()) {
        if (populations.size() < maximum_number_of_populations) {
            initializeNewPopulation();

            if (populations.size() == 1 && write_generational_statistics)
                writeGenerationalStatisticsForOnePopulation(0);

            if (populations.size() == 1 && write_generational_solutions)
                writeGenerationalSolutions(0);
        }

        generationalStepAllPopulations();

        if (populations.size() > 1 && write_generational_statistics)
            writeGenerationalStatisticsForOnePopulation(populations.size() - 1);

        if (populations.size() > 1 && write_generational_solutions)
            writeGenerationalSolutions(0);
    }
}

/**
 * Runs the IDEA.
 */
void rvg_t::run(void) {
    initialize();

    runAllPopulations();

    printf("evals %f ", fitness->number_of_evaluations);

    printf("obj_val %6.2e ", fitness->elitist_objective_value);

    printf("time %lf ", getTimer());
    printf("generations ");
    if (populations.size() > 1) {
        printf("[ ");
        for (int i = 0; i < populations.size(); i++)
            printf("%d ", populations[i]->number_of_generations);
        printf("]");
    } else
        printf("%d", populations[0]->number_of_generations);
    printf(" gomtime %.3f", gomtime);
    printf(" random_seed %ld", random_seed);
    printf(" chol_fails %d", cholesky_fails);
    printf(" cov_copied %d", cov_matrix_copy);
    printf(" elitist_improved %d", elitist_improved);
    printf("\n");
    if(print_verbose_overview)
        populations[0]->linkage_model->print();
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/**
 * The main function:
 * - interpret parameters on the command line
 * - run the algorithm with the interpreted parameters
 */
int main(int argc, char **argv) {
    rvg_t rvgomea = rvg_t(argc, argv);

    rvgomea.run();

    return (0);
}

