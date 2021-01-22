#define ERROR_EACH_LABEL_START_MUST_BE_LESS_THAN_ITS_END 1
#define ERROR_LABELED_NUMBER_OF_CHANGES_MUST_BE_0_OR_1 2
#define ERROR_EACH_LABEL_START_MUST_BE_ON_OR_AFTER_PREVIOUS_END 3
#define ERROR_LABEL_START_MUST_BE_ZERO_OR_LARGER 4
#define ERROR_LABEL_END_MUST_BE_LESS_THAN_N_DATA 5
#define ERROR_PENALTY_MUST_BE_NON_NEGATIVE 6
#define ERROR_NO_DATA 7
#define ERROR_DATA_MUST_BE_FINITE 8

int LOPART
(double *input_data,
 int n_data, //N in paper
 int *input_label_start,
 int *input_label_end,
 int *input_label_changes,
 int n_labels,//M in paper
 double penalty_unlabeled,//lambda.
 double penalty_labeled,//lambda.
 int n_updates,//size of out_ arrays and number of dp updates.
 //inputs above, outputs below.
 double *out_cumsum,//for computing optimal cost of a segment.
 int *out_change_candidates,//T_t
 double *out_cost_candidates,// for visualization.
 double *out_cost, //out_cost[t-1] = W_t in paper, for t=1 to N.
 double *out_mean,
 int *out_last_change //out_last_change[t-1] tau*_t in paper.
 );
