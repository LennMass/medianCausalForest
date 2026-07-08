/*
 * split.Rule = Least Median Square Regression (set temporarily as CT)
 */
#include <math.h>
#include "causalTree.h"
#include "causalTreeproto.h"
#include "Medianfunction4.h"

static double *sums, *wtsums, *treatment_effect;
static double *wts, *trs, *trsums;
static int *countn;
static int *tsplit;
static double *wtsqrsums, *trsqrsums;

int
userinit(int n, double *y[], int maxcat, char **error,
         int *size, int who, double *wt, double *treatment, 
         int bucketnum, int bucketMax, double *train_to_est_ratio)
{
	if (who == 1 && maxcat > 0) {
		graycode_init0(maxcat);
		countn = (int *) ALLOC(2 * maxcat, sizeof(int));
		tsplit = countn + maxcat;
		treatment_effect = (double *) ALLOC(8 * maxcat, sizeof(double));
		wts = treatment_effect + maxcat;
		trs = wts + maxcat;
		sums = trs + maxcat;
		wtsums = sums + maxcat;
		trsums = wtsums + maxcat;
		wtsqrsums = trsums + maxcat;
		trsqrsums = wtsqrsums + maxcat;
	}
	*size = 1;
	*train_to_est_ratio = n * 1.0 / ct.NumHonest;
	return 0;
}


void
userss(int n, double *y[], double *value,  double *con_mean, double *tr_mean, 
       double *risk, double *wt, double *treatment, double max_y,
       double alpha, double train_to_est_ratio)
{
	int i,r=0, k=0;
	double temp0 = 0., temp1 = 0., twt = 0., temp = 0.; /* sum of the weights */ 
	double ttreat = 0.;
	double effect;

	double tr_var, con_var;
	double con_sqr_sum = 0., tr_sqr_sum = 0.;
	double med_array[n];
	for (i = 0; i < n; i++) {
	temp1 += *y[i] * wt[i] * treatment[i]; 
	temp0 += *y[i] * wt[i] * (1 - treatment[i]);
	twt += wt[i];
	ttreat += wt[i] * treatment[i];
	
	tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
	con_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * (1- treatment[i]);
	
	}
	double trmean =  temp1 / ttreat;
	double conmean = temp0 / (twt - ttreat);
	for (i = 0; i < n; i++) {
		double mu_i = treatment[i] ? trmean : conmean;
		double res = (*y[i] - mu_i);
		med_array[i] = res * res; 
	}
	effect = temp1 / ttreat - temp0 / (twt - ttreat);

	tr_var = tr_sqr_sum / ttreat - temp1 * temp1 / (ttreat * ttreat);
	con_var = con_sqr_sum / (twt - ttreat) - temp0 * temp0 / ((twt - ttreat) * (twt - ttreat));

	*tr_mean = temp1 / ttreat;
	*con_mean = temp0 / (twt - ttreat);
	*value = effect;
	*risk = 4 * twt * max_y * max_y + alpha * twt * quick_select5(med_array, n);
}

void user(int n, double *y[], double *x, int nclass, int edge, double *improve, double *split, 
          int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
          double train_to_est_ratio)
{
	int i, j, u; //m=0,p=0, r=0, k=0;
	double temp;
	double left_sum, right_sum;
	double left_tr_sum, right_tr_sum;
	double left_tr, right_tr;
	double left_wt, right_wt;
	int left_n, right_n;
	double best;
	int direction = LEFT;
	int where = 0;
	double node_effect, left_effect, right_effect;
	double left_temp, right_temp;
	int min_node_size = minsize;
	
	int tr_n = 0;
	double tr_var, con_var;
	double right_sqr_sum, right_tr_sqr_sum, left_sqr_sum, left_tr_sqr_sum;
	double left_tr_var, left_con_var, right_tr_var, right_con_var;
	double med_array_right_all[n], med_array_right[n], med_array_left[n];  
	
	right_wt = 0.;
	right_tr = 0.;
	right_sum = 0.;
	right_tr_sum = 0.;
	
	right_sqr_sum = 0.;
	right_tr_sqr_sum = 0.;
	
	right_n = n;
	
	for (i = 0; i < n; i++) {
		right_wt += wt[i];
		right_tr += wt[i] * treatment[i];
		right_sum += *y[i] * wt[i];
		right_tr_sum += *y[i] * wt[i] * treatment[i];
		right_sqr_sum += (*y[i]) * (*y[i]) * wt[i];
		right_tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
	}
	
	double trmean = right_tr_sum / right_tr;
	double conmean = (right_sum - right_tr_sum) / (right_wt - right_tr);
	
	temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
	tr_var = right_tr_sqr_sum / right_tr - right_tr_sum * right_tr_sum / (right_tr * right_tr);
	con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
		- (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
		/ ((right_wt - right_tr) * (right_wt - right_tr));
		
		for (i = 0; i < n; i++) {
			double mu_i = (treatment[i] > 0.5) ? trmean : conmean;
			double res  = (*y[i]) - mu_i;
			med_array_right_all[i] = res * res; 
		}
		node_effect = myrisk;
		if (nclass == 0) {
			/* continuous predictor */
			left_wt = 0;
			left_tr = 0;
			left_n = 0;
			left_sum = 0;
			left_tr_sum = 0;
			best = 0;
			
			for (i = 0; right_n > edge; i++) {
				left_wt += wt[i];
				right_wt -= wt[i];
				
				left_tr += wt[i] * treatment[i];
				right_tr -= wt[i] * treatment[i];
				
				left_n++;
				right_n--;
				
				temp = *y[i] * wt[i] * treatment[i];
				
				left_tr_sum += temp;
				right_tr_sum -= temp;
				
				left_sum += *y[i] * wt[i];
				right_sum -= *y[i] * wt[i];
				////////////////////
				if (x[i + 1] != x[i] && left_n >= edge &&
        (int) left_tr >= min_node_size &&
        (int) left_wt - (int) left_tr >= min_node_size &&
        (int) right_tr >= min_node_size &&
        (int) right_wt - (int) right_tr >= min_node_size) {
					double left_trmean = left_tr_sum / left_tr;
					double left_conmean = (left_sum - left_tr_sum) / (left_wt - left_tr);
					left_temp = left_tr_sum / left_tr - 
						(left_sum - left_tr_sum) / (left_wt - left_tr);
					
					for (u = 0; u <= i; u++) {
						double mu_u = (treatment[u] > 0.5) ? left_trmean : left_conmean;
						double res  = (*y[u]) - mu_u;
						med_array_left[u] = res * res;
					}
					left_effect = alpha * left_wt * quick_select5(med_array_left, left_n);
					double right_trmean = right_tr_sum / right_tr;
					double right_conmean = (right_sum - right_tr_sum) / (right_wt - right_tr);
					right_temp = right_tr_sum / right_tr -
						(right_sum - right_tr_sum) / (right_wt - right_tr);
					
					int m = 0;
					for (u = i+1; u < n; u++) {
						double mu_u = (treatment[u] > 0.5) ? right_trmean : right_conmean;
						double res  = (*y[u]) - mu_u;
						med_array_right[m++] = res * res;
					}
					right_effect = alpha * right_wt * quick_select5(med_array_right, right_n); 
					
					temp = (left_effect + right_effect) - node_effect;
					if (temp > best) {
						best = temp;
						where = i;               
						if (left_temp < right_temp)
							direction = LEFT;
						else
							direction = RIGHT;
					}             
					
				}
			}
			
			*improve = best;
			if (best > 0) {         /* found something */
			csplit[0] = direction;
				*split = (x[where] + x[where + 1]) / 2; 
			}
		}
		
		/*
		 * Categorical predictor
		 */
		else {
			for (i = 0; i < nclass; i++) {
				countn[i] = 0;
				wts[i] = 0;
				trs[i] = 0;
				sums[i] = 0;
				wtsums[i] = 0;
				trsums[i] = 0;
				wtsqrsums[i] = 0;
				trsqrsums[i] = 0;
			}
			
			/* rank the classes by treatment effect */
			for (i = 0; i < n; i++) {
				j = (int) x[i] - 1;
				countn[j]++;
				wts[j] += wt[i];
				trs[j] += wt[i] * treatment[i];
				sums[j] += *y[i];
				wtsums[j] += *y[i] * wt[i];
				trsums[j] += *y[i] * wt[i] * treatment[i];
				wtsqrsums[j] += (*y[i]) * (*y[i]) * wt[i];
				trsqrsums[j] +=  (*y[i]) * (*y[i]) * wt[i] * treatment[i];
			}
			
			for (i = 0; i < nclass; i++) {
				if (countn[i] > 0) {
					tsplit[i] = RIGHT;
					treatment_effect[i] = trsums[i] / trs[i] - (wtsums[i] - trsums[i]) / (wts[i] - trs[i]);
				} else
					tsplit[i] = 0;
			}
			graycode_init2(nclass, countn, treatment_effect);
			
			/*
			 * Now find the split that we want
			 */
			
			left_wt = 0;
			left_tr = 0;
			left_n = 0;
			left_sum = 0;
			left_tr_sum = 0;
			left_sqr_sum = 0.;
			left_tr_sqr_sum = 0.;
			//	r=0;
			//	k=0;
			
			best = 0;
			where = 0;
			while ((j = graycode()) < nclass) {
				tsplit[j] = LEFT;
				left_n += countn[j];
				right_n -= countn[j];
				
				left_wt += wts[j];
				right_wt -= wts[j];
				
				left_tr += trs[j];
				right_tr -= trs[j];
				
				left_sum += wtsums[j];
				right_sum -= wtsums[j];
				
				left_tr_sum += trsums[j];
				right_tr_sum -= trsums[j];
				
				left_sqr_sum += wtsqrsums[j];
				right_sqr_sum -= wtsqrsums[j];
				
				left_tr_sqr_sum += trsqrsums[j];
				right_tr_sqr_sum -= trsqrsums[j];
				
				if (left_n >= edge && right_n >= edge &&
        (int) left_tr >= min_node_size &&
        (int) left_wt - (int) left_tr >= min_node_size &&
        (int) right_tr >= min_node_size &&
        (int) right_wt - (int) right_tr >= min_node_size) {
					double left_trmean = left_tr_sum / left_tr;
					double left_conmean = (left_sum - left_tr_sum) / (left_wt - left_tr);
					left_temp = left_tr_sum / left_tr - (left_sum - left_tr_sum)
						/ (left_wt - left_tr);
					
					for (u = 0; u <= i; u++) {
						double mu_u = (treatment[u] > 0.5) ? left_trmean : left_conmean;
						double res  = (*y[u]) - mu_u;
						med_array_left[u] = res * res;
					}
					left_effect = alpha * left_wt * quick_select5(med_array_left, left_n);
					
					double right_trmean = right_tr_sum / right_tr;
					double right_conmean = (right_sum - right_tr_sum) / (right_wt - right_tr);
					right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) 
						/ (right_wt - right_tr);
					int m = 0;
					for (u = i+1; u < n; u++) {
						double mu_u = (treatment[u] > 0.5) ? right_trmean : right_conmean;
						double res  = (*y[u]) - mu_u;
						med_array_right[m++] = res * res;
					}
					right_effect = alpha * right_wt * quick_select5(med_array_right, right_n); 
					temp = (left_effect + right_effect) - node_effect;
					if (temp > best) {
						best = temp;
						
						if (left_temp > right_temp)
							for (i = 0; i < nclass; i++) csplit[i] = -tsplit[i];
						else
							for (i = 0; i < nclass; i++) csplit[i] = tsplit[i];
					}
				}
			}
			*improve = best;
		}
}


double
	userpred(double *y, double wt, double treatment, double *yhat, double propensity)
	{
		double ystar;
		double temp;
		
		ystar = y[0] * (treatment - propensity) / (propensity * (1 - propensity));
		temp = ystar - *yhat;
		return temp * temp * wt;
	}

// initial user.c file content in htetree package


// /*
//  * split.Rule = user (set temporarily as CT)
//  */
// #include <math.h>
// #include "causalTree.h"
// #include "causalTreeproto.h"
// 
// static double *sums, *wtsums, *treatment_effect;
// static double *wts, *trs, *trsums;
// static int *countn;
// static int *tsplit;
// static double *wtsqrsums, *trsqrsums;
// 
// int
// userinit(int n, double *y[], int maxcat, char **error,
// 		int *size, int who, double *wt, double *treatment, 
// 		int bucketnum, int bucketMax, double *train_to_est_ratio)
// {
// 	if (who == 1 && maxcat > 0) {
// 		graycode_init0(maxcat);
// 		countn = (int *) ALLOC(2 * maxcat, sizeof(int));
// 		tsplit = countn + maxcat;
// 		treatment_effect = (double *) ALLOC(8 * maxcat, sizeof(double));
// 		wts = treatment_effect + maxcat;
// 		trs = wts + maxcat;
// 		sums = trs + maxcat;
// 		wtsums = sums + maxcat;
// 		trsums = wtsums + maxcat;
// 		wtsqrsums = trsums + maxcat;
// 		trsqrsums = wtsqrsums + maxcat;
// 	}
// 	*size = 1;
// 	*train_to_est_ratio = n * 1.0 / ct.NumHonest;
// 	return 0;
// }
// 
// 
// void
// userss(int n, double *y[], double *value,  double *con_mean, double *tr_mean, 
// 		double *risk, double *wt, double *treatment, double max_y,
// 		double alpha, double train_to_est_ratio)
// {
// 	int i;
// 	double temp0 = 0., temp1 = 0., twt = 0.; /* sum of the weights */ 
// 	double ttreat = 0.;
// 	double effect;
// 	double tr_var, con_var;
// 	double con_sqr_sum = 0., tr_sqr_sum = 0.;
// 
// 	for (i = 0; i < n; i++) {
// 		temp1 += *y[i] * wt[i] * treatment[i];
// 		temp0 += *y[i] * wt[i] * (1 - treatment[i]);
// 		twt += wt[i];
// 		ttreat += wt[i] * treatment[i];
// 		tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
// 		con_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * (1- treatment[i]);
// 	}
// 
// 	effect = temp1 / ttreat - temp0 / (twt - ttreat);
// 	tr_var = tr_sqr_sum / ttreat - temp1 * temp1 / (ttreat * ttreat);
// 	con_var = con_sqr_sum / (twt - ttreat) - temp0 * temp0 / ((twt - ttreat) * (twt - ttreat));
// 
// 	*tr_mean = temp1 / ttreat;
// 	*con_mean = temp0 / (twt - ttreat);
// 	*value = effect;
// 	*risk = 4 * twt * max_y * max_y - alpha * twt * effect * effect + 
// 		(1 - alpha) * (1 + train_to_est_ratio) * twt * (tr_var /ttreat  + con_var / (twt - ttreat));
// }
// 
// 
// void user(int n, double *y[], double *x, int nclass, int edge, double *improve, double *split, 
// 		int *csplit, double myrisk, double *wt, double *treatment, int minsize, double alpha,
// 		double train_to_est_ratio)
// {
// 	int i, j;
// 	double temp;
// 	double left_sum, right_sum;
// 	double left_tr_sum, right_tr_sum;
// 	double left_tr, right_tr;
// 	double left_wt, right_wt;
// 	int left_n, right_n;
// 	double best;
// 	int direction = LEFT;
// 	int where = 0;
// 	double node_effect, left_effect, right_effect;
// 	double left_temp, right_temp;
// 	int min_node_size = minsize;
// 
// 	double tr_var, con_var;
// 	double right_sqr_sum, right_tr_sqr_sum, left_sqr_sum, left_tr_sqr_sum;
// 	double left_tr_var, left_con_var, right_tr_var, right_con_var;
// 
// 	right_wt = 0.;
// 	right_tr = 0.;
// 	right_sum = 0.;
// 	right_tr_sum = 0.;
// 	right_sqr_sum = 0.;
// 	right_tr_sqr_sum = 0.;
// 	right_n = n;
// 	for (i = 0; i < n; i++) {
// 		right_wt += wt[i];
// 		right_tr += wt[i] * treatment[i];
// 		right_sum += *y[i] * wt[i];
// 		right_tr_sum += *y[i] * wt[i] * treatment[i];
// 		right_sqr_sum += (*y[i]) * (*y[i]) * wt[i];
// 		right_tr_sqr_sum += (*y[i]) * (*y[i]) * wt[i] * treatment[i];
// 	}
// 
// 	temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) / (right_wt - right_tr);
// 	tr_var = right_tr_sqr_sum / right_tr - right_tr_sum * right_tr_sum / (right_tr * right_tr);
// 	con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
// 		- (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
// 		/ ((right_wt - right_tr) * (right_wt - right_tr));
// 	node_effect = alpha * temp * temp * right_wt - (1 - alpha) * (1 + train_to_est_ratio) 
// 		* right_wt * (tr_var / right_tr  + con_var / (right_wt - right_tr));
// 
// 	if (nclass == 0) {
// 		/* continuous predictor */
// 		left_wt = 0;
// 		left_tr = 0;
// 		left_n = 0;
// 		left_sum = 0;
// 		left_tr_sum = 0;
// 		left_sqr_sum = 0;
// 		left_tr_sqr_sum = 0;
// 		best = 0;
// 
// 		for (i = 0; right_n > edge; i++) {
// 			left_wt += wt[i];
// 			right_wt -= wt[i];
// 			left_tr += wt[i] * treatment[i];
// 			right_tr -= wt[i] * treatment[i];
// 			left_n++;
// 			right_n--;
// 			temp = *y[i] * wt[i] * treatment[i];
// 			left_tr_sum += temp;
// 			right_tr_sum -= temp;
// 			left_sum += *y[i] * wt[i];
// 			right_sum -= *y[i] * wt[i];
// 			temp = (*y[i]) *  (*y[i]) * wt[i];
// 			left_sqr_sum += temp;
// 			right_sqr_sum -= temp;
// 			temp = (*y[i]) * (*y[i]) * wt[i] * treatment[i];
// 			left_tr_sqr_sum += temp;
// 			right_tr_sqr_sum -= temp;
// 
// 
// 			if (x[i + 1] != x[i] && left_n >= edge &&
// 					(int) left_tr >= min_node_size &&
// 					(int) left_wt - (int) left_tr >= min_node_size &&
// 					(int) right_tr >= min_node_size &&
// 					(int) right_wt - (int) right_tr >= min_node_size) {
// 
// 				left_temp = left_tr_sum / left_tr - 
// 					(left_sum - left_tr_sum) / (left_wt - left_tr);
// 				left_tr_var = left_tr_sqr_sum / left_tr - 
// 					left_tr_sum  * left_tr_sum / (left_tr * left_tr);
// 				left_con_var = (left_sqr_sum - left_tr_sqr_sum) / (left_wt - left_tr)  
// 					- (left_sum - left_tr_sum) * (left_sum - left_tr_sum)
// 					/ ((left_wt - left_tr) * (left_wt - left_tr));        
// 				left_effect = alpha * left_temp * left_temp * left_wt
// 					- (1 - alpha) * (1 + train_to_est_ratio) * left_wt 
// 					* (left_tr_var / left_tr + left_con_var / (left_wt - left_tr));
// 
// 				right_temp = right_tr_sum / right_tr -
// 					(right_sum - right_tr_sum) / (right_wt - right_tr);
// 				right_tr_var = right_tr_sqr_sum / right_tr -
// 					right_tr_sum * right_tr_sum / (right_tr * right_tr);
// 				right_con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
// 					- (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
// 					/ ((right_wt - right_tr) * (right_wt - right_tr));
// 				right_effect = alpha * right_temp * right_temp * right_wt
// 					- (1 - alpha) * (1 + train_to_est_ratio) * right_wt * 
// 					(right_tr_var / right_tr + right_con_var / (right_wt - right_tr));
// 
// 				temp = left_effect + right_effect - node_effect;
// 				if (temp > best) {
// 					best = temp;
// 					where = i;               
// 					if (left_temp < right_temp)
// 						direction = LEFT;
// 					else
// 						direction = RIGHT;
// 				}             
// 			}
// 		}
// 
// 		*improve = best;
// 		if (best > 0) {         /* found something */
// 			csplit[0] = direction;
// 			*split = (x[where] + x[where + 1]) / 2; 
// 		}
// 	}
// 
// 	/*
// 	 * Categorical predictor
// 	 */
// 	else {
// 		for (i = 0; i < nclass; i++) {
// 			countn[i] = 0;
// 			wts[i] = 0;
// 			trs[i] = 0;
// 			sums[i] = 0;
// 			wtsums[i] = 0;
// 			trsums[i] = 0;
// 			wtsqrsums[i] = 0;
// 			trsqrsums[i] = 0;
// 		}
// 
// 		/* rank the classes by treatment effect */
// 		for (i = 0; i < n; i++) {
// 			j = (int) x[i] - 1;
// 			countn[j]++;
// 			wts[j] += wt[i];
// 			trs[j] += wt[i] * treatment[i];
// 			sums[j] += *y[i];
// 			wtsums[j] += *y[i] * wt[i];
// 			trsums[j] += *y[i] * wt[i] * treatment[i];
// 			wtsqrsums[j] += (*y[i]) * (*y[i]) * wt[i];
// 			trsqrsums[j] +=  (*y[i]) * (*y[i]) * wt[i] * treatment[i];
// 		}
// 
// 		for (i = 0; i < nclass; i++) {
// 			if (countn[i] > 0) {
// 				tsplit[i] = RIGHT;
// 				treatment_effect[i] = trsums[j] / trs[j] - (wtsums[j] - trsums[j]) / (wts[j] - trs[j]);
// 			} else
// 				tsplit[i] = 0;
// 		}
// 		graycode_init2(nclass, countn, treatment_effect);
// 
// 		/*
// 		 * Now find the split that we want
// 		 */
// 
// 		left_wt = 0;
// 		left_tr = 0;
// 		left_n = 0;
// 		left_sum = 0;
// 		left_tr_sum = 0;
// 		left_sqr_sum = 0.;
// 		left_tr_sqr_sum = 0.;
// 
// 		best = 0;
// 		where = 0;
// 		while ((j = graycode()) < nclass) {
// 			tsplit[j] = LEFT;
// 			left_n += countn[j];
// 			right_n -= countn[j];
// 
// 			left_wt += wts[j];
// 			right_wt -= wts[j];
// 
// 			left_tr += trs[j];
// 			right_tr -= trs[j];
// 
// 			left_sum += wtsums[j];
// 			right_sum -= wtsums[j];
// 
// 			left_tr_sum += trsums[j];
// 			right_tr_sum -= trsums[j];
// 
// 			left_sqr_sum += wtsqrsums[j];
// 			right_sqr_sum -= wtsqrsums[j];
// 
// 			left_tr_sqr_sum += trsqrsums[j];
// 			right_tr_sqr_sum -= trsqrsums[j];
// 
// 			if (left_n >= edge && right_n >= edge &&
// 					(int) left_tr >= min_node_size &&
// 					(int) left_wt - (int) left_tr >= min_node_size &&
// 					(int) right_tr >= min_node_size &&
// 					(int) right_wt - (int) right_tr >= min_node_size) {
// 
// 				left_temp = left_tr_sum / left_tr - (left_sum - left_tr_sum) 
// 					/ (left_wt - left_tr);
// 
// 				left_tr_var = left_tr_sqr_sum / left_tr 
// 					- left_tr_sum  * left_tr_sum / (left_tr * left_tr);
// 				left_con_var = (left_sqr_sum - left_tr_sqr_sum) / (left_wt - left_tr)  
// 					- (left_sum - left_tr_sum) * (left_sum - left_tr_sum)
// 					/ ((left_wt - left_tr) * (left_wt - left_tr));       
// 				left_effect = alpha * left_temp * left_temp * left_wt
// 					- (1 - alpha) * (1 + train_to_est_ratio) * left_wt * 
// 					(left_tr_var / left_tr + left_con_var / (left_wt - left_tr));
// 
// 				right_temp = right_tr_sum / right_tr - (right_sum - right_tr_sum) 
// 					/ (right_wt - right_tr);
// 				right_tr_var = right_tr_sqr_sum / right_tr 
// 					- right_tr_sum * right_tr_sum / (right_tr * right_tr);
// 				right_con_var = (right_sqr_sum - right_tr_sqr_sum) / (right_wt - right_tr)
// 					- (right_sum - right_tr_sum) * (right_sum - right_tr_sum) 
// 					/ ((right_wt - right_tr) * (right_wt - right_tr));
// 				right_effect = alpha * right_temp * right_temp * right_wt
// 					- (1 - alpha) * (1 + train_to_est_ratio) * right_wt *
// 					(right_tr_var / right_tr + right_con_var / (right_wt - right_tr));
// 				temp = left_effect + right_effect - node_effect;
// 
// 
// 				if (temp > best) {
// 					best = temp;
// 
// 					if (left_temp > right_temp)
// 						for (i = 0; i < nclass; i++) csplit[i] = -tsplit[i];
// 					else
// 						for (i = 0; i < nclass; i++) csplit[i] = tsplit[i];
// 				}
// 			}
// 		}
// 		*improve = best;
// 	}
// }
// 
// 
// double
// userpred(double *y, double wt, double treatment, double *yhat, double propensity)
// {
// 	double ystar;
// 	double temp;
// 
// 	ystar = y[0] * (treatment - propensity) / (propensity * (1 - propensity));
// 	temp = ystar - *yhat;
// 	return temp * temp * wt;
// }

