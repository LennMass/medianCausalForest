/*
 * rundown function fot tot
 */
#include "causalTree.h"
#include "node.h"
#include "causalTreeproto.h"

void
totrundown(pNode tree, int obs, double *cp, double *xpred, double *xtemp)
{
    int i, obs2 = (obs < 0) ? -(1 + obs) : obs;
    pNode otree =  tree;

    int opnumber = 0;

    for (i = 0; i < ct.num_unique_cp; i++) {
        while (cp[i] < tree->complexity) {
	        tree = branch(tree, obs);
	        if (tree == 0)
		    goto oops;
	        otree = tree;
	    }

        xpred[i] = tree->response_est[0];

        xtemp[i] = (*ct_xeval)(ct.ydata[obs2], tree->response_est,ct.wt[obs2], ct.treatment[obs2],
                ct.propensity);
        // xtemp[i] = (*ct_xeval)(ct.ydata[obs2],0,0,tree->response_est,ct.wt[obs2],0,0,0,0,0,0,ct.propensity,0,ct.treatment[obs2],0,0,0,0,0,ct.wt[obs2],0,0,0);
    }
    return;

oops:;
    if (ct.usesurrogate < 2) {  /* must have hit a missing value */
	for (i = 0; i < ct.num_unique_cp; i++)
	    xpred[i] = otree->response_est[0];

	xtemp[i] = (*ct_xeval)(ct.ydata[obs2], ct.wt[obs2], ct.treatment[obs2], tree->response_est,
             ct.propensity);
	// xtemp[i] = (*ct_xeval)(ct.ydata[obs2],0,0,tree->response_est,0,0,0,0,0,0,propensity,0,ct.treatment[obs2],0,0,0,0,0,ct.wt[obs2],0,0,0);
	Rprintf("oops number %d.\n", opnumber++);
  return;
    }
    warning("Warning message--see rundown.c");
}
