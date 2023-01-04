#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS    
#endif

#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>

#include <omp.h>
#include <Grishagin/GrishaginProblemFamily.hpp>
#include <Grishagin/GrishaginConstrainedProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <GKLS/GKLSConstrainedProblemFamily.hpp>
#include <mggsa.h>
#include <task.h>

using namespace std;

// #define CALC

const int type = 0; // 0 - count trials, 1 - count points, 2 - c_points / c_trials
const int family_number = 3; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained),
const int number_family = 4;
int current_func[4];

TGrishaginProblemFamily grishaginProblems;
double f_grishagin(vector<double> x, int j) {
    switch (j) {
        case 1: return grishaginProblems[current_func[0]]->ComputeFunction(x);
        default: return numeric_limits<double>::quiet_NaN();
    }
};

TGKLSProblemFamily GKLSProblems;
double f_gkls(vector<double> x, int j) {
    switch (j) {
        case 1: return GKLSProblems[current_func[1]]->ComputeFunction({ x });
        default: return numeric_limits<double>::quiet_NaN();
    }
};

TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
double f_constr_grishagin(vector<double> x, int j) {
    int constr = grishaginConstrainedProblems[current_func[2]]->GetConstraintsNumber();
    if (j >= 1 && j <= constr) {
        return grishaginConstrainedProblems[current_func[2]]->ComputeConstraint(j - 1, x);
    } else if (j - 1 == constr) {
        return grishaginConstrainedProblems[current_func[2]]->ComputeFunction(x);
    } else {
        return numeric_limits<double>::quiet_NaN();
    }
};

TGKLSConstrainedProblemFamily GKLSConstrainedProblems;
double f_constr_gkls(vector<double> x, int j) {
    int constr = GKLSConstrainedProblems[current_func[3]]->GetConstraintsNumber();
    if (j >= 1 && j <= constr) {
        return GKLSConstrainedProblems[current_func[3]]->ComputeConstraint(j - 1, x);
    } else if (j - 1 == constr) {
        return GKLSConstrainedProblems[current_func[3]]->ComputeFunction(x);
    } else {
        return numeric_limits<double>::quiet_NaN();
    }
};

int main() {
#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics_incr_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_operational_characteristics_incr_test_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    int count_func;

    int start_time, end_time;
    double work_time;
    int chunk = 1;
    int num_threads = max(number_family, omp_get_num_procs());

    int count_trials;

    int den = 10, m = 0;
    double eps = 0.01;
    int Nmax = 30000, key = 3;

    vector<double> A, B, X_opt;
    vector<vector<double>> r_array{ {3.0, 2.4, 1.6, 1.0},
                                    {4.2, 3.8, 2.0, 1.0},
                                    {3.5, 2.6, 1.8, 1.0},
                                    {4.7, 3.0, 2.0, 1.0} };
    vector<int> incr_array{0, 20, 40, 60};
    vector<double> d_array{0.0, 0.0, 0.01, 0.01};

    vector<int> count_trials_vec;
    vector<int> count_points_vec;

    vector<vector<vector<int>>> max_count_trials(number_family);
    vector<vector<vector<int>>> max_count_points(number_family);
    for (int i = 0; i < number_family; i++) {
        max_count_trials[i].resize(r_array[i].size());
        max_count_points[i].resize(r_array[i].size());
        for (int j = 0; j < r_array[i].size(); j++) {
            max_count_trials[i][j].resize(incr_array.size());
            max_count_points[i][j].resize(incr_array.size());
        }
    }

    vector<class_problems_fm> problems{ class_problems_fm("GrishaginProblemFamily", &grishaginProblems, type_constraned::NONCONSTR,
                                                          f_grishagin, "Grishagin"),
                                        class_problems_fm("GKLSProblemFamily", &GKLSProblems, type_constraned::NONCONSTR, 
                                                          f_gkls, "GKLS"),
                                        class_problems_fm("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems, 
                                                          type_constraned::CONSTR, f_constr_grishagin, "GrishaginConstrained"),
                                        class_problems_fm("GKLSProblemConstrainedFamily", &GKLSConstrainedProblems, 
                                                          type_constraned::CONSTR, f_constr_gkls, "GKLSConstrained") };

    mggsa_method mggsa(nullptr, -1, -1, A, B, -1.0, -1.0, den, key, eps, Nmax, -1);

#pragma omp parallel for schedule(static, chunk) proc_bind(close) num_threads(num_threads) \
        shared(current_func, problems, d_array, r_array, incr_array, max_count_trials, max_count_points), firstprivate(mggsa) \
        private(A, B, count_func, start_time, end_time, work_time, X_opt, count_trials_vec, count_points_vec, count_trials)
    for (int i = 0; i < number_family; i++) {
        IOptProblemFamily *opt_problem_family;
        IConstrainedOptProblemFamily *constr_opt_problem_family;
        if (problems[i].type == type_constraned::CONSTR) {
            constr_opt_problem_family = static_cast<IConstrainedOptProblemFamily*>(problems[i].problem);
            mggsa.setN((*constr_opt_problem_family)[0]->GetDimension());
            mggsa.setM((*constr_opt_problem_family)[0]->GetConstraintsNumber());
            (*constr_opt_problem_family)[0]->GetBounds(A, B);
            count_trials_vec.resize(constr_opt_problem_family->GetFamilySize());
            count_points_vec.resize(constr_opt_problem_family->GetFamilySize());
        } else {
            opt_problem_family = static_cast<IOptProblemFamily*>(problems[i].problem);
            mggsa.setN((*opt_problem_family)[0]->GetDimension());
            mggsa.setM(0);
            (*opt_problem_family)[0]->GetBounds(A, B);
            count_trials_vec.resize(opt_problem_family->GetFamilySize());
            count_points_vec.resize(opt_problem_family->GetFamilySize());
        }
        mggsa.setAB(A, B);
        mggsa.setF(problems[i].f);
        mggsa.setD(d_array[i]);
        count_func = problems[i].problem->GetFamilySize();

        for (int j = 0; j < r_array[i].size(); j++) {
            mggsa.setR(r_array[i][j]);
            for (int k = 0; k < incr_array.size(); k++) {
                mggsa.setIncr(incr_array[k]);
                start_time = clock();
                for (int l = 0; l < count_func; l++) {
                    current_func[i] = l;
                    if (problems[i].type == type_constraned::CONSTR) {
                        X_opt = (*constr_opt_problem_family)[l]->GetOptimumPoint();
                    } else {
                        X_opt = (*opt_problem_family)[l]->GetOptimumPoint();
                    }
                    if (mggsa.solve_test(X_opt, count_trials, Stop::ACCURNUMBER)) {
                        count_trials_vec[l] = count_trials;
                    } else {
                        count_trials_vec[l] = count_trials + 1;
                    }
                    count_points_vec[l] = mggsa.getCountPoints();
                }
                end_time = clock();
                work_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

                auto iter = max_element(count_trials_vec.begin(), count_trials_vec.end());
                max_count_trials[i][j][k] = *iter;
                iter = max_element(count_points_vec.begin(), count_points_vec.end());
                max_count_points[i][j][k] = *iter;

                string str = problems[i].name + " r = " + to_string(r_array[i][j]) + " incr = " + to_string(incr_array[k]) + 
                             " count trials = " + to_string(max_count_trials[i][j][k]) + " count points = " + 
                             to_string(max_count_points[i][j][k]) + " time: " + to_string(work_time) + " t_num: " + to_string(omp_get_thread_num()) + "\n";
                cout << str;
            }
        }
    }

    for (int i = 0; i < number_family; i++) {
        for (int j = 0; j < r_array[i].size(); j++) {
            for (int k = 0; k < incr_array.size(); k++) {
                ofstr << incr_array[k] << " " << max_count_trials[i][j][k] << " " << max_count_points[i][j][k] << endl;
            }
            ofstr << endl << endl;
        }
        ofstr << endl << endl;
    }
    ofstr.close();

    int size = incr_array.size();
    ofstr_opt << "count_key = " << size << endl;
    ofstr_opt << "array Name[" << number_family << "]" << endl;
    ofstr_opt << "array R[" << size * number_family << "]" << endl;
    for (int i = 0; i < number_family; i++) {
        ofstr_opt << "Name[" << i + 1 << "]=\"" << problems[i].short_name << "\"" << endl;
        for (int j = 0; j < r_array[i].size(); j++) {
            ofstr_opt << "R[" << (i * size) + j + 1 << "]=\"" << r_array[i][j] << "\"" << endl; 
        }
    }
    ofstr_opt.close();
#endif

    // Plotting operational characteristics(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/mggsa_operational_characteristics_incr_test.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/mggsa_operational_characteristics_incr_test.gp %d %d", type, family_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif
	return 0;
}
