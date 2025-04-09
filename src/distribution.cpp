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

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "distribution.h"



#include <armadillo>
template<class Matrix>
void print_matrix(Matrix matrix) {
    matrix.print(std::cout);
}

//provide explicit instantiations of the template function for 
//every matrix type you use somewhere in your program.
template void print_matrix<arma::mat>(arma::mat matrix);
template void print_matrix<arma::cx_mat>(arma::cx_mat matrix);

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
distribution_t::~distribution_t() {
}

void distribution_t::adaptDistributionMultiplier(partial_solution_t **partial_solutions, int num_solutions) {
    short improvementForFOSElement = 0;
    if ((((double) out_of_bounds_draws) / ((double) samples_drawn)) > 0.9)
        distribution_multiplier *= 0.5;

    improvementForFOSElement = generationalImprovementForOnePopulationForFOSElement(partial_solutions, num_solutions, &this->st_dev_ratio);

    if (improvementForFOSElement) {
        if (distribution_multiplier < 1.0)
            distribution_multiplier = 1.0;

        if (this->st_dev_ratio > st_dev_ratio_threshold)
            distribution_multiplier *= distribution_multiplier_increase;
    } else {
        if (distribution_multiplier > 1.0)
            distribution_multiplier *= distribution_multiplier_decrease;

        if (distribution_multiplier < 1.0)
            distribution_multiplier = 1.0;
    }
}

void
distribution_t::updateConditionals(const std::map<int, std::set<int>> &variable_interaction_graph, int visited[]) {}

void distribution_t::setOrder(const std::vector<int> &order) {
    exit(1);
}

void distribution_t::print() {}

void
distribution_t::adaptDistributionMultiplierMaximumStretch(partial_solution_t **partial_solutions, int num_solutions) {
    short improvementForFOSElement = 0;
    if ((((double) out_of_bounds_draws) / ((double) samples_drawn)) > 0.9)
        distribution_multiplier *= 0.5;

    double st_dev_ratio;
    improvementForFOSElement = generationalImprovementForOnePopulationForFOSElement(partial_solutions, num_solutions,
                                                                                    &st_dev_ratio);

    if (improvementForFOSElement) {
        if (distribution_multiplier < 1.0)
            distribution_multiplier = 1.0;

        if (st_dev_ratio > st_dev_ratio_threshold)
            distribution_multiplier *= distribution_multiplier_increase;
    } else {
        distribution_multiplier *= distribution_multiplier_decrease;
    }
}

double distribution_t::estimateMean(int var, solution_t **selection, int selection_size) {
    double mean = 0.0;
    for (int j = 0; j < selection_size; j++)
        mean += selection[j]->variables[var];
    mean /= (double) selection_size;
    return (mean);
}

double distribution_t::estimateCovariance(int vara, int varb, solution_t **selection, int selection_size) {
    double cov = 0.0;
    double meana = estimateMean(vara, selection, selection_size);
    double meanb = estimateMean(varb, selection, selection_size);
    for (int m = 0; m < selection_size; m++)
        cov += (selection[m]->variables[vara] - meana) * (selection[m]->variables[varb] - meanb);
    cov /= (double) selection_size;
    return (cov);
}

vec distribution_t::estimateMeanVectorML(std::vector<int> variables, solution_t **selection, int selection_size) {
    vec mean_vector = vec(variables.size(), fill::none);
    for (int i = 0; i < variables.size(); i++)
        mean_vector[i] = estimateMean(variables[i], selection, selection_size);
    return (mean_vector);
}

mat distribution_t::estimateUnivariateCovarianceMatrixML(std::vector<int> variables, solution_t **selection,
                                                         int selection_size) {
    /* First do the maximum-likelihood estimate from data */
    mat covariance_matrix = mat(variables.size(), variables.size(), fill::zeros);
    for (int j = 0; j < variables.size(); j++) {
        int vara = variables[j];
        double cov = estimateCovariance(vara, vara, selection, selection_size);
        //covariance_matrix(j,k) = (1-eta_cov)*covariance_matrix(j,k)+ eta_cov*cov;
        covariance_matrix(j, j) = cov * distribution_multiplier;
    }
    return (covariance_matrix);
}

mat
distribution_t::estimateRegularCovarianceMatrixML(std::vector<int> variables, vec mean_vector, solution_t **selection,
                                                  int selection_size) {
    mat covariance_matrix;
    covariance_matrix = estimateCovarianceMatrixML(variables, selection, selection_size);
    /*int n = variables.size();
    if( selection_size < n + 1 )
        covariance_matrix = estimateUnivariateCovarianceMatrixML(variables,selection,selection_size);
    else
    {
        covariance_matrix = estimateCovarianceMatrixML(variables,selection,selection_size);
        if( selection_size < (0.5*n*(n+1))+1 )
            regularizeCovarianceMatrix(covariance_matrix,mean_vector,selection,selection_size);
    }*/
    return (covariance_matrix);
}

mat distribution_t::estimateCovarianceMatrixML(std::vector<int> variables, solution_t **selection, int selection_size) {
    /* First do the maximum-likelihood estimate from data */
    //mat covariance_matrix(variables.size(),variables.size(),fill::none);
    mat covariance_matrix = mat(variables.size(), variables.size());
    for (int j = 0; j < variables.size(); j++) {
        int vara = variables[j];
        for (int k = j; k < variables.size(); k++) {
            int varb = variables[k];
            double cov = estimateCovariance(vara, varb, selection, selection_size);
            //covariance_matrix(j,k) = (1-eta_cov)*covariance_matrix(j,k)+ eta_cov*cov;
            covariance_matrix(j, k) = cov;
            covariance_matrix(k, j) = covariance_matrix(j, k);
        }
    }
    return (covariance_matrix);
}

bool
distribution_t::regularizeCovarianceMatrix(mat &cov_mat, vec &mean_vector, solution_t **selection, int selection_size) {
    // regularization for small populations
    double number_of_samples = (double) selection_size;
    int n = variables.size();

    // either use the univariate matrix as a prior,
    // or a diagonal matrix with the mean variance on all diagonal entries
    bool use_univariate_as_prior = true;

    double meanvar = 0.0;
    if (!use_univariate_as_prior) {
        for (int i = 0; i < n; ++i) {
            meanvar += cov_mat(i, i);
        }
        meanvar /= (double) n;
    }

    double phi = 0.0;
    // y = x.^2
    // phiMat = y'*y/t-sample.^2
    // phi = sum(sum(phiMat))
    mat squared_cov = mat(n, n, fill::none);
    double temp;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            squared_cov(i, j) = 0.0;
            for (int k = 0; k < selection_size; ++k) {
                temp = (selection[k]->variables[variables[i]] - mean_vector[i]) *
                       (selection[k]->variables[variables[j]] - mean_vector[j]);
                squared_cov(i, j) += temp * temp;
            }
            squared_cov(i, j) /= number_of_samples;
        }
    }

    // this can be implemented faster by considering only half this matrix,
    // and we dont need to store square_cov actually.
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            phi += squared_cov(i, j) - cov_mat(i, j) * cov_mat(i, j);

    // Frobenius norm, i.e.,
    // gamma = norm(sample - prior,'fro')^2;
    double gamma = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (use_univariate_as_prior) {
                temp = fabs(cov_mat(i, j) - ((i == j) ? cov_mat(i, i) : 0.0));
            } else {
                temp = fabs(cov_mat(i, j) - ((i == j) ? meanvar : 0.0));
            }
            gamma += temp * temp;
        }
    }

    double kappa = phi / gamma;
    double shrinkage = std::max(0.0, std::min(1.0, kappa / number_of_samples));
    //std::cout << "Shrinkage with factor " << shrinkage << std::endl;
    //shrinkage = fmax(1e-10,shrinkage);

    if (shrinkage == 0.0) {
        return false;
    }

    // shrinking
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (use_univariate_as_prior) {
                cov_mat(i, j) = (1.0 - shrinkage) * cov_mat(i, j) + ((i == j) ? shrinkage * cov_mat(i, i) : 0.0);
            } else {
                cov_mat(i, j) = (1.0 - shrinkage) * cov_mat(i, j) + ((i == j) ? shrinkage * meanvar : 0.0);
            }
        }
    }

    return true;
}

mat distribution_t::pseudoInverse(const mat &matrix) {
    mat result;
    if (!pinv(result, matrix)) {
        cholesky_fails++;
        result = mat(matrix.n_cols, matrix.n_cols, fill::zeros);
        for (int i = 0; i < matrix.n_cols; i++) {
            if (matrix(i, i) == 0) result(i, i) = 1e38;
            else result(i, i) = 1.0 / matrix(i, i);
        }
    }
    return (result);
}

/**
 * BLAS subroutine.
 */
int distribution_t::blasDSWAP(int n, double *dx, int incx, double *dy, int incy) {
    double dtmp;

    if (n > 0) {
        incx *= sizeof(double);
        incy *= sizeof(double);

        dtmp = (*dx);
        *dx = (*dy);
        *dy = dtmp;

        while ((--n) > 0) {
            dx = (double *) ((char *) dx + incx);
            dy = (double *) ((char *) dy + incy);
            dtmp = (*dx);
            *dx = (*dy);
            *dy = dtmp;
        }
    }

    return (0);
}

/**
 * BLAS subroutine.
 */
int distribution_t::blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy) {
    double dtmp0, dtmp, *dx0, *dy0;

    if (n > 0 && da != 0.) {
        incx *= sizeof(double);
        incy *= sizeof(double);
        *dy += da * (*dx);

        if ((n & 1) == 0) {
            dx = (double *) ((char *) dx + incx);
            dy = (double *) ((char *) dy + incy);
            *dy += da * (*dx);
            --n;
        }
        n = n >> 1;
        while (n > 0) {
            dy0 = (double *) ((char *) dy + incy);
            dy = (double *) ((char *) dy0 + incy);
            dtmp0 = (*dy0);
            dtmp = (*dy);
            dx0 = (double *) ((char *) dx + incx);
            dx = (double *) ((char *) dx0 + incx);
            *dy0 = dtmp0 + da * (*dx0);
            *dy = dtmp + da * (*dx);
            --n;
        }
    }

    return (0);
}

/**
 * BLAS subroutine.
 */
void distribution_t::blasDSCAL(int n, double sa, double x[], int incx) {
    int i, ix, m;

    if (n <= 0) {
    } else if (incx == 1) {
        m = n % 5;

        for (i = 0; i < m; i++) {
            x[i] = sa * x[i];
        }

        for (i = m; i < n; i = i + 5) {
            x[i] = sa * x[i];
            x[i + 1] = sa * x[i + 1];
            x[i + 2] = sa * x[i + 2];
            x[i + 3] = sa * x[i + 3];
            x[i + 4] = sa * x[i + 4];
        }
    } else {
        if (0 <= incx) {
            ix = 0;
        } else {
            ix = (-n + 1) * incx;
        }

        for (i = 0; i < n; i++) {
            x[ix] = sa * x[ix];
            ix = ix + incx;
        }
    }
}

/**
 * LINPACK subroutine.
 */
int distribution_t::linpackDCHDC(double a[], int lda, int p, double work[], int ipvt[]) {
    int info, j, jp, k, l, maxl, pl, pu;
    double maxdia, temp;

    pl = 1;
    pu = 0;
    info = p;
    for (k = 1; k <= p; k++) {
        maxdia = a[k - 1 + (k - 1) * lda];
        maxl = k;
        if (pl <= k && k < pu) {
            for (l = k + 1; l <= pu; l++) {
                if (maxdia < a[l - 1 + (l - 1) * lda]) {
                    maxdia = a[l - 1 + (l - 1) * lda];
                    maxl = l;
                }
            }
        }

        if (maxdia <= 0.0) {
            info = k - 1;

            return (info);
        }

        if (k != maxl) {
            blasDSWAP(k - 1, a + 0 + (k - 1) * lda, 1, a + 0 + (maxl - 1) * lda, 1);

            a[maxl - 1 + (maxl - 1) * lda] = a[k - 1 + (k - 1) * lda];
            a[k - 1 + (k - 1) * lda] = maxdia;
            jp = ipvt[maxl - 1];
            ipvt[maxl - 1] = ipvt[k - 1];
            ipvt[k - 1] = jp;
        }
        work[k - 1] = sqrt(a[k - 1 + (k - 1) * lda]);
        a[k - 1 + (k - 1) * lda] = work[k - 1];

        for (j = k + 1; j <= p; j++) {
            if (k != maxl) {
                if (j < maxl) {
                    temp = a[k - 1 + (j - 1) * lda];
                    a[k - 1 + (j - 1) * lda] = a[j - 1 + (maxl - 1) * lda];
                    a[j - 1 + (maxl - 1) * lda] = temp;
                } else if (maxl < j) {
                    temp = a[k - 1 + (j - 1) * lda];
                    a[k - 1 + (j - 1) * lda] = a[maxl - 1 + (j - 1) * lda];
                    a[maxl - 1 + (j - 1) * lda] = temp;
                }
            }
            a[k - 1 + (j - 1) * lda] = a[k - 1 + (j - 1) * lda] / work[k - 1];
            work[j - 1] = a[k - 1 + (j - 1) * lda];
            temp = -a[k - 1 + (j - 1) * lda];

            blasDAXPY(j - k, temp, work + k, 1, a + k + (j - 1) * lda, 1);
        }
    }

    return (info);
}

/**
 * Computes the lower-triangle Cholesky Decomposition
 * of a square, symmetric and positive-definite matrix.
 * Subroutines from LINPACK and BLAS are used.
 */
mat distribution_t::choleskyDecomposition(const mat &matrix) {
    int i, j, k, info, *ipvt;
    double *a, *work;

    int n = matrix.n_rows;
    a = (double *) Malloc(n * n * sizeof(double));
    work = (double *) Malloc(n * sizeof(double));
    ipvt = (int *) Malloc(n * sizeof(int));

    k = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            a[k] = matrix(i, j);
            k++;
        }
        ipvt[i] = 0;
    }

    info = linpackDCHDC(a, n, n, work, ipvt);

    mat result = mat(n, n, fill::none);
    if (info != n) /* Matrix is not positive definite */
    {
        cholesky_fails++;
        k = 0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (i == j) {
                    double matrix_value = matrix(i, j);
                    if (matrix_value <= 0) {
                        matrix_value = 1;
                    }
                    result(i, j) = sqrt(matrix_value);
                } else {
                    result(i, j) = 0.0;
                }

                k++;
            }
        }
    } else {
        k = 0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                result(i, j) = i < j ? 0.0 : a[k];
                k++;
            }
        }
    }

    free(ipvt);
    free(work);
    free(a);

    return (result);
}

normal_distribution_t::normal_distribution_t(std::vector<int> variables) {
    this->variables = variables;
    this->copied_covariance_matrix = false;
    this->distribution_multiplier_decrease = use_incremental ? 0.95 : 0.9;
}


void normal_distribution_t::copyCovariances(const normal_distribution_t * source) {
    if( source->covariance_matrix.size() == 0 ) return;
    int n_vars = variables.size();
    int n_vars_other = variables.size();
    if(!copied_covariance_matrix){
        copied_covariance_matrix = true;
        cov_matrix_copy++;
        covariance_matrix = mat(n_vars, n_vars, fill::zeros);
    }
    // Copy information to correct cell in covariance matrices
    for(int p_i = 0; p_i < n_vars; p_i ++){
        for(int p_j = p_i; p_j < n_vars; p_j ++){
            for(int s_i = 0; s_i < n_vars_other; s_i ++){
                for(int s_j = 0; s_j < n_vars_other; s_j ++){
                    if(source->variables[s_i] == this->variables[p_i] && source->variables[s_j] == this->variables[p_j])
                        covariance_matrix(p_i, p_j) = source->covariance_matrix(s_i, s_j);
                }
            }
        }
    }
}

void normal_distribution_t::estimateDistribution(solution_t **selection, int selection_size,
                                                 vec_t <vec_t<double>> fitness_dependency_matrix){
    mean_vector = estimateMeanVectorML(variables, selection, selection_size);

    /* Change the focus of the search to the best solution */
    if (distribution_multiplier < 1.0)
        for (int j = 0; j < variables.size(); j++)
            mean_vector[j] = selection[0]->variables[variables[j]];

    if(!copied_covariance_matrix){
        if(use_incremental){
            covariance_matrix = mat(variables.size(), variables.size(), fill::zeros);
            for(int i = 0; i < variables.size(); i++)
                covariance_matrix(i, i) = estimateCovariance(variables[i], variables[i], selection, selection_size);
        } else {
            covariance_matrix = estimateRegularCovarianceMatrixML(variables, mean_vector, selection, selection_size);
        }
        copied_covariance_matrix = true; // For possible reuse of this distribution in case of static linkage models
    } else {
        double eta_cov = use_incremental ? 1.0 - exp(alpha_cov[0] * pow(selection_size, alpha_cov[1]) / pow(variables.size(), alpha_cov[2])) : 1.0;
        covariance_matrix = (covariance_matrix * (1.0 - eta_cov)) + (eta_cov * estimateRegularCovarianceMatrixML(variables, mean_vector, selection, selection_size));
    }
    //int c1 = cholesky_fails;
    cholesky_decomposition = choleskyDecomposition(covariance_matrix * distribution_multiplier);
    /*if( c1 != cholesky_fails )
    {
        cholesky_fails--;
        regularizeCovarianceMatrix(covariance_matrix,mean_vector,selection,selection_size);
        cholesky_decomposition = choleskyDecomposition( covariance_matrix );
        //if(cholesky_fails!=c1)printf("BLA\n");
    }*/
}

void normal_distribution_t::estimateDistributionAMS(int selection_size, double * mean_shift_vector){
    this->mean_shift_vector = vec(variables.size(), fill::none);

    /* Change the focus of the search to the best solution */
    for (int j = 0; j < variables.size(); j++)
        this->mean_shift_vector[j] = mean_shift_vector[variables[j]];
}

partial_solution_t *normal_distribution_t::generatePartialSolution(solution_t *parent) {
    std::vector<int> indices = variables;
    int num_indices = variables.size();
    vec result = vec(num_indices, fill::none);

    int times_not_in_bounds = -1;
    out_of_bounds_draws--;

    short ready = 0;
    do {
        times_not_in_bounds++;
        samples_drawn++;
        out_of_bounds_draws++;

        if (times_not_in_bounds >= 100) {
            printf("Sampled out of bounds too many times.\n");
            exit(1);
            /*result = vec(num_indices,fill::none);
            sample_means = vec(num_indices,fill::none);
            sample_zs = zeros<vec>(num_indices);
            for(int i = 0; i < num_indices; i++ )
            {
                result[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]])*randu<double>();
                sample_means[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]])*0.5;
            }*/
        } else {
            vec sample_zs = random1DNormalUnitVector(num_indices);
            result = mean_vector + cholesky_decomposition * sample_zs;
        }

        ready = 1;
        /*for(int i = 0; i < num_indices; i++ )
        {
            if( !fitness->isParameterInRangeBounds( result[i], indices[i] ) )
            {
                ready = 0;
                break;
            }
        }*/ // BLA
    } while (!ready);

    partial_solution_t *res_sol = new partial_solution_t(result, indices);
    res_sol->setSampleMean(mean_vector);
    return (res_sol);
}

void normal_distribution_t::applyPartialAMS(partial_solution_t *solution, double delta_AMS)
{
    short out_of_range = 1;
    double shrink_factor = 2;
    double *result = (double *)Malloc(solution->num_touched_variables * sizeof(double));
    while ((out_of_range == 1) && (shrink_factor > 1e-10))
    {
        shrink_factor *= 0.5;
        out_of_range = 0;
        for (int m = 0; m < solution->num_touched_variables; m++)
        {
            result[m] = solution->touched_variables[m] + shrink_factor * delta_AMS * this->distribution_multiplier * (mean_shift_vector[m]);
            // if (!fitness->isParameterInRangeBounds(result[m], im))
            // {
            //     out_of_range = 1;
            //     break;
            // }
        }
    }
    if (!out_of_range)
    {
        for (int m = 0; m < solution->num_touched_variables; m++)
        {
            int im = solution->touched_indices[m];
            solution->touched_variables[m] = result[m];
        }
    }
    free(result);
}

/*short normal_distribution_t::generationalImprovementForOnePopulationForFOSElement( partial_solution_t** partial_solutions, int num_solutions, double *st_dev_ratio )
{
	short generational_improvement = 0;
    std::vector<int> indices = partial_solutions[0]->touched_indices;
	int num_indices = indices.size();

	vec average_parameters_of_improvements = zeros(num_indices);

	int number_of_improvements  = 0;
	mat cholinv = pinv( trimatl( cholesky_decomposition ) );
	for(int i = 0; i < num_solutions; i++ )
	{
		if( partial_solutions[i]->improves_elitist )
		{
			number_of_improvements++;
			for(int j = 0; j < num_indices; j++ )
				average_parameters_of_improvements[j] += partial_solutions[i]->touched_variables[j] - mean_vector[j];
		}
	}

	// Determine st.dev. ratio
	*st_dev_ratio = 0.0;
	if( number_of_improvements > 0 )
	{
		for(int i = 0; i < num_indices; i++ )
			average_parameters_of_improvements[i] /= (double) number_of_improvements;

		vec zs = cholinv * average_parameters_of_improvements;
		for(int i = 0; i < num_indices; i++ )
			*st_dev_ratio = fmax( *st_dev_ratio, fabs(zs[i]) );

		generational_improvement = 1;
	}

	return( generational_improvement );
}*/

short
normal_distribution_t::generationalImprovementForOnePopulationForFOSElement(partial_solution_t **partial_solutions,
                                                                            int num_solutions, double *st_dev_ratio) {
    short generational_improvement = 0;
    std::vector<int> indices = partial_solutions[0]->touched_indices;
    int num_indices = indices.size();

    double average_z_of_improvements[num_indices];
    for (int i = 0; i < num_indices; i++)
        average_z_of_improvements[i] = 0.0;

    int number_of_improvements = 0;
    mat cholinv = pinv(trimatl(cholesky_decomposition));
    for (int i = 0; i < num_solutions; i++) {
        if (partial_solutions[i]->improves_elitist) {
            number_of_improvements++;
            vec z = cholinv * (partial_solutions[i]->touched_variables - partial_solutions[i]->sample_means);
            for (int j = 0; j < num_indices; j++) {
                average_z_of_improvements[j] += z[j];
            }
        }
    }

    // Determine st.dev. ratio
    *st_dev_ratio = 0.0;
    if (number_of_improvements > 0) {
        for (int i = 0; i < num_indices; i++) {
            average_z_of_improvements[i] /= (double) number_of_improvements;
            *st_dev_ratio = fmax(*st_dev_ratio, fabs(average_z_of_improvements[i]));
        }

        generational_improvement = 1;
    }

    return (generational_improvement);
}

conditional_distribution_t::conditional_distribution_t() {
    this->distribution_multiplier_decrease = use_incremental ? 0.95 : 0.9;
}

conditional_distribution_t::conditional_distribution_t(const std::vector<int> &variables,
                                                       const std::vector<int> &conditioned_variables) {
    addGroupOfVariables(variables, conditioned_variables);
    /*for( int x : variables )
    {
        std::vector<int> vars;
        vars.push_back(x);
        addGroupOfVariables(vars,conditioned_variables);
    }*/
}


conditional_distribution_t::conditional_distribution_t(const std::vector<int> &variables,
                                                       const std::set<int> &conditioned_variables) {
    addGroupOfVariables(variables, conditioned_variables);
    /*for( int x : variables )
    {
        std::vector<int> vars;
        vars.push_back(x);
        addGroupOfVariables(vars,conditioned_variables);
    }*/
}

void conditional_distribution_t::initializeMemory()
{
    variables_copied.resize(variable_groups.size());
    mean_vectors.resize(variable_groups.size());
    mean_shift_vectors.resize(variable_groups.size());
    mean_vectors_conditioned_on.resize(variable_groups.size());
    copied_covariance_matrix.resize(variable_groups.size(), false);
    full_covariance_matrices.resize(variable_groups.size());
    covariance_matrices_multiplied.resize(variable_groups.size());
    rho_matrices.resize(variable_groups.size());
    cholesky_decompositions.resize(variable_groups.size());

    for(int i = 0; i < variable_groups.size(); i++){
        variables_copied[i] = std::vector(variable_groups[i].size() + variables_conditioned_on[i].size(), false);
    }
}

void conditional_distribution_t::addGroupOfVariables(std::vector<int> indices, std::vector<int> indices_cond) {
    std::sort(indices.begin(), indices.end());
    std::sort(indices_cond.begin(), indices_cond.end());
    std::vector<int> indices_map;
    for (int i: indices) {
        indices_map.push_back(variables.size());
        variables.push_back(i);
    }
    for (int i: indices_cond) {
        conditionals.push_back(i);
    }
    //std::sort(variables.begin(),variables.end());
    index_in_var_array.push_back(indices_map);
    variable_groups.push_back(indices);
    variables_conditioned_on.push_back(indices_cond);
    order.push_back(order.size());
}

void
conditional_distribution_t::addGroupOfVariables(const std::vector<int> &indices, const std::set<int> &indices_cond) {
    std::vector<int> cond;
    for (int i: indices_cond)
        cond.push_back(i);
    addGroupOfVariables(indices, cond);
}

void conditional_distribution_t::estimateConditionalAMSML(int variable_group_index, int selection_size, double * mean_shift_vectors)
{
    int n = variable_groups[variable_group_index].size();
    if (copied_covariance_matrix[variable_group_index])
    {
        double eta_ams = use_incremental ? 1.0 - exp(alpha_ams[0] * pow(selection_size, alpha_ams[1]) / pow(n, alpha_ams[2])) : 1.0;
        for(int i = 0; i < n; i++){
            this->mean_shift_vectors[variable_group_index][i] = (1.0 - eta_ams) * this->mean_shift_vectors[variable_group_index][i] + eta_ams * mean_shift_vectors[variable_groups[variable_group_index][i]];
        }
    } else {
        this->mean_shift_vectors[variable_group_index].clear();
        this->mean_shift_vectors[variable_group_index].resize(n);
        for(int i = 0; i < n; i++){
            this->mean_shift_vectors[variable_group_index][i] = mean_shift_vectors[variable_groups[variable_group_index][i]];
        }
    }
}

void conditional_distribution_t::estimateConditionalGaussianML(int variable_group_index, solution_t **selection, int selection_size,
                                                               vec_t<vec_t<double>> fitness_dependency_matrix)
{
    int i = variable_group_index;
    std::vector<int> vars = variable_groups[i];
    int n = vars.size();



    mean_vectors[i] = estimateMeanVectorML(vars, selection, selection_size);

    /* Change the focus of the search to the best solution */
    if (distribution_multiplier < 1.0)
        for (int j = 0; j < n; j++)
            mean_vectors[i][j] = selection[0]->variables[vars[j]];


    std::vector<int> vars_cond = variables_conditioned_on[i];
    int n_cond = vars_cond.size();
    mat conditional_covariance;
    if (n_cond > 0) {
        mean_vectors_conditioned_on[i] = estimateMeanVectorML(vars_cond, selection, selection_size);

        mat A11(n, n, fill::none);
        mat A12(n, n_cond, fill::none);
        mat A22(n_cond, n_cond, fill::none);

        // std::set<int> all_vars_set = std::set<int>(vars.begin(), vars.end());
        // all_vars_set.insert(vars_cond.begin(), vars_cond.end());
        std::vector<int> all_vars = std::vector<int>(vars.begin(), vars.end());
        for(int param_cond : vars_cond){
            if(all_vars.end() == std::find(all_vars.begin(), all_vars.end(), param_cond))
                all_vars.push_back(param_cond);
            else {
                std::cout << "ERROR: variable " << param_cond << " is in both vars and vars_cond" << std::endl;
                assert(false);
            }
        }

        vec mean_vector_all = estimateMeanVectorML(all_vars, selection, selection_size);
        if(use_incremental){
            if(this->copied_covariance_matrix[i]){
                // Regular update of the covariance matrix
                double eta_cov = 1.0 - exp(alpha_cov[0] * pow(selection_size, alpha_cov[1]) / pow(all_vars.size(), alpha_cov[2]));
                full_covariance_matrices[i] = (this->full_covariance_matrices[i] * (1.0 - eta_cov)) + (estimateRegularCovarianceMatrixML(variable_group_index, all_vars, mean_vector_all, selection, selection_size) * eta_cov);

                // Initialize any missing variables to simple variances
                for(int param = 0; param < (n_cond + n); param++){
                    if(!variables_copied[i][param]){
                        for(int param_j = 0; param_j < (n_cond + n); param_j++){
                            if(param_j == param){
                                full_covariance_matrices[i](param,param) = estimateCovariance(all_vars[param], all_vars[param], selection, selection_size);
                            } else {
                                full_covariance_matrices[i](param,param_j) = 0.0f;
                                full_covariance_matrices[i](param_j,param) = 0.0f;
                            }
                        }
                    }
                }
            } else {
                // Initialize the covariance matrix as one containing only simple variances
                full_covariance_matrices[i] = mat(all_vars.size(), all_vars.size(), fill::zeros);
                for(int j = 0; j < all_vars.size(); j++){
                    full_covariance_matrices[i](j,j) = estimateCovariance(all_vars[j], all_vars[j], selection, selection_size);
                }
            }

            // Extract the submatrices from the full covariance matrix
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    A11(j, k) = full_covariance_matrices[i](j,k);
                }
            }

            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n_cond; k++) {
                    A12(j, k) = full_covariance_matrices[i](j,k+n);
                }
            }

            for (int j = 0; j < n_cond; j++) {
                for (int k = 0; k < n_cond; k++) {
                    A22(j, k) = full_covariance_matrices[i](j+n,k+n);
                }
            }
        } else {
            full_covariance_matrices[i] = estimateRegularCovarianceMatrixML(variable_group_index, all_vars, mean_vector_all, selection, selection_size);
            A11 = estimateRegularCovarianceMatrixML(variable_group_index, vars, mean_vectors[i], selection, selection_size);
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n_cond; k++) {
                    A12(j, k) = estimateCovariance( vars[j], vars_cond[k], selection, selection_size );// * distribution_multiplier;
                }
            }

            A22 = estimateRegularCovarianceMatrixML(variable_group_index, vars_cond, mean_vectors_conditioned_on[i], selection, selection_size);
        }

        mat A22inv;
        if (pinv(A22inv, A22)) {
            rho_matrices[i] = A12 * A22inv;
            mat submat = A12 * A22inv * A12.t();
            conditional_covariance = A11 - submat;
        } else {
            printf("pseudo-inverse failed\n");
            fflush(stdout);
        }
    } else {
        if(use_incremental){
            if(this->copied_covariance_matrix[i]){
                double eta_cov = 1.0 - exp(alpha_cov[0] * pow(selection_size, alpha_cov[1]) / pow(vars.size(), alpha_cov[2]));
                full_covariance_matrices[i] = (full_covariance_matrices[i] * (1.0 - eta_cov)) + (estimateRegularCovarianceMatrixML(variable_group_index, vars, mean_vectors[i], selection, selection_size) * eta_cov);
                for(int param = 0; param < n; param++){
                    if(!variables_copied[i][param]){
                        for(int param_j = 0; param_j < n; param_j++){
                            if(param_j == param)
                                full_covariance_matrices[i](param,param) = estimateCovariance(vars[param], vars[param], selection, selection_size);
                            else {
                                full_covariance_matrices[i](param,param_j) = 0.0f;
                                full_covariance_matrices[i](param_j,param) = 0.0f;
                            }
                        }
                    }
                }
            } else {
                full_covariance_matrices[i] = mat(vars.size(), vars.size(), fill::zeros);
                for(int j = 0; j < vars.size(); j++){
                    full_covariance_matrices[i](j,j) = estimateCovariance(vars[j], vars[j], selection, selection_size);
                }
            }
        } else {
            full_covariance_matrices[i] = estimateRegularCovarianceMatrixML(variable_group_index, vars, mean_vectors[i], selection, selection_size);
        }
        conditional_covariance = full_covariance_matrices[i];
    }

    // For possible reuse of this distribution in case of static linkage models
    this->copied_covariance_matrix[i] = true;
    this->variables_copied[i] = std::vector(n, true);

    // Get the final decomposition
    covariance_matrices_multiplied[i] = conditional_covariance * distribution_multiplier;
    cholesky_decompositions[i] = choleskyDecomposition(covariance_matrices_multiplied[i]);
}

void conditional_distribution_t::updateConditionals(const std::map<int, std::set<int>> &variable_interaction_graph,
                                                    int visited[]) {
    const int IS_VISITED = 1;
    const int IN_CLIQUE = 2;

    variables_conditioned_on.clear();
    variables_conditioned_on.resize(variable_groups.size());

    for (int i = 0; i < variable_groups.size(); i++) {
        int ind = i;
        if (order.size() > 0)
            ind = order[i];
        std::vector<int> clique = variable_groups[ind];

        // Add FOS element of all nodes in clique, conditioned on dependent, already visited variables
        std::set<int> cond;
        for (int v: clique)
            visited[v] = IN_CLIQUE;
        for (int v: clique) {
            for (int x: variable_interaction_graph.at(v)) {
                if (visited[x] == IS_VISITED)
                    cond.insert(x);
            }
        }
        for (int v: clique)
            visited[v] = IS_VISITED;

        std::vector<int> cond_vec;
        for (int x: cond)
            cond_vec.push_back(x);
        variables_conditioned_on[ind] = cond_vec;
    }
}

void conditional_distribution_t::copyCovariances(const conditional_distribution_t * source) {
    if( source->full_covariance_matrices.size() == 0 ) return;
    for (int i = 0; i < variable_groups.size(); i++){
        int fit_idx = -1;
        for (int j = 0; j < source->variable_groups.size(); j++){
            if(std::includes(variable_groups[i].begin(), variable_groups[i].end(), source->variable_groups[j].begin(), source->variable_groups[j].end()) &&
               std::includes(variables_conditioned_on[i].begin(), variables_conditioned_on[i].end(), source->variables_conditioned_on[j].begin(), source->variables_conditioned_on[j].end())){
                fit_idx = j;
                break;
            }
        }
        if(fit_idx >= 0){
            int n_vars = variable_groups[i].size();
            int n_cond = variables_conditioned_on[i].size();
            int n_vars_other = source->variable_groups[fit_idx].size();
            int n_cond_other = source->variables_conditioned_on[fit_idx].size();
            if(!copied_covariance_matrix[i]){
                copied_covariance_matrix[i] = true;
                cov_matrix_copy++;
                full_covariance_matrices[i] = mat(n_vars + n_cond, n_vars + n_cond, fill::zeros);
            }
            // Copy information to correct cell in covariance matrices
            for(int p_i = 0; p_i < n_vars; p_i ++){
                for(int p_j = 0; p_j < n_vars; p_j ++){
                    for(int s_i = 0; s_i < n_vars_other; s_i ++){
                        for(int s_j = 0; s_j < n_vars_other; s_j ++){
                            if(source->variable_groups[fit_idx][s_i] == this->variable_groups[i][p_i] && source->variable_groups[fit_idx][s_j] == this->variable_groups[i][p_j]){
                                full_covariance_matrices[i](p_i, p_j) = source->full_covariance_matrices[fit_idx](s_i, s_j);
                                if(p_i == p_j && !variables_copied[i][p_i]){
                                    variables_copied[i][p_i] = true;
                                }
                            }
                        }
                    }
                }
            }
            for(int p_i = 0; p_i < n_vars; p_i ++){
                for(int p_j = 0; p_j < n_cond; p_j ++){
                    for(int s_i = 0; s_i < n_vars_other; s_i ++){
                        for(int s_j = 0; s_j < n_cond_other; s_j ++){
                            if(source->variable_groups[fit_idx][s_i] == this->variable_groups[i][p_i] && source->variables_conditioned_on[fit_idx][s_j] == this->variables_conditioned_on[i][p_j]){
                                full_covariance_matrices[i](p_i, n_vars + p_j) = source->full_covariance_matrices[fit_idx](s_i, n_vars_other + s_j);
                                full_covariance_matrices[i](n_vars + p_j, p_i) = source->full_covariance_matrices[fit_idx](n_vars_other + s_j, s_i);
                            }
                        }
                    }
                }
            }
            for(int p_i = 0; p_i < n_cond; p_i ++){
                for(int p_j = 0; p_j < n_cond; p_j ++){
                    for(int s_i = 0; s_i < n_cond_other; s_i ++){
                        for(int s_j = 0; s_j < n_cond_other; s_j ++){
                            if(source->variables_conditioned_on[fit_idx][s_i] == this->variables_conditioned_on[i][p_i] && source->variables_conditioned_on[fit_idx][s_j] == this->variables_conditioned_on[i][p_j]){
                                full_covariance_matrices[i](n_vars + p_i, n_vars + p_j) = source->full_covariance_matrices[fit_idx](n_vars_other + s_i, n_vars_other + s_j);
                                if(p_i == p_j && !variables_copied[i][n_vars + p_i]){
                                    variables_copied[i][n_vars + p_i] = true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

mat
conditional_distribution_t::estimateRegularCovarianceMatrixML(int group_idx, std::vector<int> variables, vec mean_vector, solution_t **selection,
                                                  int selection_size) {
    mat covariance_matrix;
    covariance_matrix = estimateCovarianceMatrixML(group_idx, variables, selection, selection_size);
    /*int n = variables.size();
    if( selection_size < n + 1 )
        covariance_matrix = estimateUnivariateCovarianceMatrixML(variables,selection,selection_size);
    else
    {
        covariance_matrix = estimateCovarianceMatrixML(variables,selection,selection_size);
        if( selection_size < (0.5*n*(n+1))+1 )
            regularizeCovarianceMatrix(covariance_matrix,mean_vector,selection,selection_size);
    }*/
    return (covariance_matrix);
}

mat conditional_distribution_t::estimateCovarianceMatrixML(int group_idx, std::vector<int> variables, solution_t **selection, int selection_size) {
    /* First do the maximum-likelihood estimate from data */
    //mat covariance_matrix(variables.size(),variables.size(),fill::none);
    mat covariance_matrix = mat(variables.size(), variables.size());
    for (int j = 0; j < variables.size(); j++) {
        int vara = variables[j];
        for (int k = j; k < variables.size(); k++) {
            int varb = variables[k];
            double cov = estimateCovariance(vara, varb, selection, selection_size);
            covariance_matrix(j, k) = cov;
            covariance_matrix(k, j) = cov;
        }
    }
    return (covariance_matrix);
}

void conditional_distribution_t::estimateDistributionAMS(int selection_size, double * mean_shift_vector)
{
    // initializeMemory();
    for (int i = 0; i < variable_groups.size(); i++)
        estimateConditionalAMSML(i, selection_size, mean_shift_vector);
}

void conditional_distribution_t::estimateDistribution(solution_t **selection, int selection_size,
                                                      vec_t<vec_t<double>> fitness_dependency_matrix)
{
    if(this->mean_vectors.size() != variable_groups.size()){
        initializeMemory();
    }
    for (int i = 0; i < variable_groups.size(); i++)
        estimateConditionalGaussianML(i, selection, selection_size, fitness_dependency_matrix);
}

void conditional_distribution_t::applyPartialAMS(partial_solution_t *solution, double delta_AMS)
{
    vec result = vec(variables.size(), fill::none);
    std::map<int, int> sampled_indices;
    assert(order.size() == variable_groups.size());

    for (int k = 0; k < variable_groups.size(); k++)
    {
        int og = order[k];
        std::vector<int> indices = variable_groups[og];
        int num_indices = indices.size();

        int times_not_in_bounds = -1;
        out_of_bounds_draws--;

        vec ams_result;
        short ready = 0;
        do
        {
            times_not_in_bounds++;
            samples_drawn++;
            out_of_bounds_draws++;

            if (times_not_in_bounds >= 100)
            {
                printf("Sampled out of bounds too many times.\n");
                exit(1);
                /*vec sample_result = vec(num_indices, fill::none);
                vec sample_means = vec(num_indices, fill::none);
                for(int i = 0; i < num_indices; i++ )
                {
                    sample_result[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]])*randomRealUniform01();
                    sample_means[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]]) * 0.5;
                }*/
            }
            else
            {
                ams_result = delta_AMS * this->distribution_multiplier * this->mean_shift_vectors[og];
            }

            ready = 1;
            /*for(int i = 0; i < num_indices; i++ )
            {
                if( !fitness->isParameterInRangeBounds( sample_result[i], indices[i] ) )
                {
                    ready = 0;
                    break;
                }
            }*/
            // BLA
        } while (!ready);

        for (int i = 0; i < num_indices; i++)
        {
            assert(indices[i] == variables[index_in_var_array[og][i]]);
            result[index_in_var_array[og][i]] = ams_result[i];
        }
    }

    solution->touched_variables += result;
}

partial_solution_t *conditional_distribution_t::generatePartialSolution(solution_t *solution_conditioned_on) {
    vec result = vec(variables.size(), fill::none);
    vec means = vec(variables.size(), fill::none);
    std::map<int, int> sampled_indices;
    assert(order.size() == variable_groups.size());

    for (int k = 0; k < variable_groups.size(); k++) {
        int og = order[k];
        std::vector<int> indices = variable_groups[og];
        int num_indices = indices.size();

        int times_not_in_bounds = -1;
        out_of_bounds_draws--;

        vec sample_result;
        vec sample_means;
        short ready = 0;
        do {
            times_not_in_bounds++;
            samples_drawn++;
            out_of_bounds_draws++;

            if (times_not_in_bounds >= 100) {
                printf("Sampled out of bounds too many times.\n");
                exit(1);
                /*vec sample_result = vec(num_indices, fill::none);
                vec sample_means = vec(num_indices, fill::none);
                for(int i = 0; i < num_indices; i++ )
                {
                    sample_result[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]])*randomRealUniform01();
                    sample_means[i] = lower_init_ranges[indices[i]] + (upper_init_ranges[indices[i]] - lower_init_ranges[indices[i]]) * 0.5;
                }*/
            } else {
                sample_means = mean_vectors[og];

                std::vector<int> indices_cond = variables_conditioned_on[og];
                int num_indices_cond = indices_cond.size();
                /*printf("means_nc ");
                for( double x : sample_means )
                    printf("%10.3e ",x);
                printf("\n");*/
                if (num_indices_cond > 0) {
                    vec cond = vec(num_indices_cond, fill::none);
                    for (int i = 0; i < num_indices_cond; i++) {
                        auto it = sampled_indices.find(indices_cond[i]);
                        if (it != sampled_indices.end()) {
                            assert(variables[it->second] == indices_cond[i]);
                            cond[i] = result[it->second];
                        } else
                            cond[i] = solution_conditioned_on->variables[indices_cond[i]];
                    }
                    sample_means += rho_matrices[og] * (cond - mean_vectors_conditioned_on[og]);
                    /*printf("means_co ");
                    for( double x : sample_means )
                        printf("%10.3e ",x);
                    printf("\n");*/
                }
                //printf("\n");
                vec sample_zs = random1DNormalUnitVector(num_indices);
                sample_result = sample_means + cholesky_decompositions[og] * sample_zs;
            }

            ready = 1;
            /*for(int i = 0; i < num_indices; i++ )
            {
                if( !fitness->isParameterInRangeBounds( sample_result[i], indices[i] ) )
                {
                    ready = 0;
                    break;
                }
            }*/ // BLA
        } while (!ready);

        for (int i = 0; i < num_indices; i++) {
            assert(indices[i] == variables[index_in_var_array[og][i]]);
            result[index_in_var_array[og][i]] = sample_result[i];
            means[index_in_var_array[og][i]] = sample_means[i];
            sampled_indices[indices[i]] = index_in_var_array[og][i];
        }
    }

    partial_solution_t *sol_res = new partial_solution_t(result, variables);
    sol_res->setSampleMean(means);
    return (sol_res);
}

short
conditional_distribution_t::generationalImprovementForOnePopulationForFOSElement(partial_solution_t **partial_solutions,
                                                                                 int num_solutions,
                                                                                 double *st_dev_ratio) {
    *st_dev_ratio = 0.0;
    short generational_improvement = 0;

    for (int k = 0; k < variable_groups.size(); k++) {
        std::vector<int> indices = variable_groups[k];
        int num_indices = variable_groups[k].size();
        int number_of_improvements = 0;

        double average_z_of_improvements[num_indices];
        for (int i = 0; i < num_indices; i++)
            average_z_of_improvements[i] = 0.0;

        mat cholinv = pseudoInverse(trimatl(cholesky_decompositions[k]));
        vec sample_means(num_indices, fill::none);
        for (int i = 0; i < num_solutions; i++) {
            if (partial_solutions[i]->improves_elitist) {
                number_of_improvements++;
                for (int j = 0; j < num_indices; j++) {
                    int ind = index_in_var_array[k][j];
                    sample_means[j] =
                            partial_solutions[i]->touched_variables[ind] - partial_solutions[i]->sample_means[ind];
                }
                vec z = cholinv * sample_means; //(partial_solutions[i]->touched_variables - partial_solutions[i]->sample_means);
                for (int j = 0; j < num_indices; j++)
                    average_z_of_improvements[j] += z[j];
            }
        }

        // Determine st.dev. ratio
        if (number_of_improvements > 0) {
            for (int i = 0; i < num_indices; i++) {
                average_z_of_improvements[i] /= (double) number_of_improvements;
                *st_dev_ratio = fmax(*st_dev_ratio, fabs(average_z_of_improvements[i]));
            }
            generational_improvement = 1;
        }
    }

    return (generational_improvement);
}

void conditional_distribution_t::setOrder(const std::vector<int> &order) {
    this->order = order;
}

void conditional_distribution_t::print() {
    for (int i = 0; i < variable_groups.size(); i++) {
        int og = order[i];
        printf("[");
        for (int x: variable_groups[og])
            printf(" %d", x);
        printf("]->[");
        for (int x: variables_conditioned_on[og])
            printf(" %d", x);
        printf("],");
    }
    printf("\n");
}
