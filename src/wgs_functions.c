//
//  wgs_functions.c
//  
//
//  Created by Colin Worby on 08/05/2015.
//
//

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "wgsmcmc_head.h"

extern int num_patients;
extern int studyspan;
extern int TP;
extern int num_sequences;
extern int allneg;

// total infected days avoiding colonisation
int avoideddays(struct augdata *a_data, int from, int to) {
    int i = 0, tot=0;
    for (i=from;i<=to;i++) {
        tot += col_pop(a_data, i);
    }
    return(tot);
}

// colonised population on a particular day
int col_pop(struct augdata *a_data, int day) {
    int i=0, cp=0;
    for (i=0;i<num_patients;i++) {
        if (a_data -> col_t <= day && a_data -> col_t !=0 && a_data -> day_d >= day) { // present and has col time prior to day
            if (a_data -> col_t == day && a_data -> onadm == 1) { // pos on admission today
                cp++;
            } else if (a_data -> col_t < day) { // col time earlier than today
                cp++;
            }
        }
        a_data++;
    }
    return(cp);
}

// Number of imported cases already
int out_pop(struct augdata *a_data, int day) {
    int op=0;
    int i;
    for(i=0;i<num_patients;i++){
       if (a_data -> day_a < day && a_data -> onadm == 1) { //already admitted, pos on adm
          //if (a_data -> obspos == 1) { // observed strain
             op++;
          //}
       }
       a_data++;
    }
    return(op);
}

// test result loglikelihood contribution
double likelihoodimport(double p, int num_onadm) {
    double cont = 0.0;
    cont = num_onadm*log(p)+(num_patients-num_onadm)*log(1.0-p);
    return(cont);
}

// test observation loglikelihood contribution
double likelihoodtests(struct augdata *a_data, double z, int *resultsmatrix) {
    double cont=0.0;
    cont = TP*log(z) + Num_FN(a_data,resultsmatrix)*log(1-z);
    return(cont);
    
}

// transmission dynamic loglikelihood contribution
double likelihoodtrans(struct augdata *a_data, double beta) {
    int i=0;
    double ll=0.0;
    
    for (i=0;i<num_patients;i++) { //  cycle through patients
        if (a_data -> psource!=0 && a_data -> onadm == 0) { // acquisition
            ll += log(1-exp(-beta*col_pop(a_data-i, a_data -> col_t))) - log(col_pop(a_data-i, a_data -> col_t)); // col day
            if (col_pop(a_data-i, a_data -> col_t)==0) {
                Rprintf("Zero col pop at col time %d (pid %d)\n", a_data -> col_t, a_data -> patID);
            }
            if (a_data -> col_t > a_data -> day_a) {
                ll -= beta*avoideddays(a_data-i,a_data -> day_a, a_data -> col_t-1); // before col day
            }
        } else if (a_data -> onadm == 0) { // avoided infection
            ll -= beta*avoideddays(a_data-i,a_data -> day_a, a_data -> day_d);
        }
        a_data++;
    }
    return(ll);
}

// genetic distance loglikelihood contribution: importation structure model

double likelihoodgen_IS(struct augdata *a_data, int dmat[][num_sequences], double gamma, double gamma_G, double clust_par, int *patseqID) {
    int i=0, j=0, k=0, grps=0;
    double gll = 0.0;
    
    for (i=0;i<num_patients;i++) { // cycle through patients
        if (a_data -> col_t>0) { // positive?
            for (j=0; j<num_sequences; j++) { // look for corresponding sequence
                if (patseqID[j] == a_data -> patID) { // found a sequence for patient i
                    for (k=0; k<j; k++) { // go through all gen distances to previous sequences
                        if (dmat[j][k]>10000) {
                            Rprintf("Extreme distance: j=%d, k=%d, seqID=%d\n", j, k, patseqID[j]);
                        }
                        if (samegroup(a_data -> patID, patseqID[k], a_data-i)==1) { // belong to same group
                            gll += log(gamma) + dmat[j][k]*log(1.0-gamma);
                        } else {
                            gll += log(gamma_G) + dmat[j][k]*log(1.0-gamma_G);
                        }
                    }
                }
            }
        }
        //if (a_data -> onadm == 1 && a_data -> group == a_data -> patID) { //unclustered
        //  gll += log(1-clust_par);
        //} else if (a_data -> onadm == 1 && a_data -> group != a_data -> patID) { // clustered
        //  gll += log(clust_par) - log(out_pop(a_data-i, a_data -> day_a));
        //}
        a_data++;
    }
    
    grps = tot_groups(a_data-num_patients+1);
    gll += grps*log(clust_par) + (num_imports(a_data-num_patients+1)-grps)*log(1-clust_par);
    return(gll);
}

// genetic distance loglikelihood contribution: transmission diversity model

double likelihoodgen_TD(struct augdata *a_data, int dmat[][num_sequences], double gamma, double gamma_G, double trans_par, int *patseqID) {
    int i=0, j=0, k=0;
    double gll = 0.0;
    int ge=0;
    for (i=0;i<num_patients;i++) { // cycle through patients
        if (a_data -> col_t>0) { // positive?
            for (j=0; j<num_sequences; j++) { // look for corresponding sequence
                if (patseqID[j] == a_data -> patID) { // found a sequence for patient i
                    for (k=0; k<j; k++) { // go through all gen distances to previous sequences
                        if (dmat[j][k]>10000) {
                            Rprintf("Extreme distance: j=%d, k=%d, seqID=%d\n", j, k, patseqID[j]);
                        }
                        if (samegroup(a_data -> patID, patseqID[k], a_data-i)==0) { // belong to different groups (unlinked)
                            gll += log(gamma_G) + dmat[j][k]*log(1.0-gamma_G);
                        } else if (translinks(a_data -> patID, patseqID[k], a_data-i)==0) { // same host
                            gll += log(gamma) + dmat[j][k]*log(1.0-gamma);
                        } else if (translinks(a_data -> patID, patseqID[k], a_data-i)!=-1) { // same chain
                            ge = translinks(a_data -> patID, patseqID[k], a_data-i);
                            gll += log(gamma)+ge*log(trans_par) + dmat[j][k]*log(1.0-gamma*pow(trans_par,ge));
                        } else {
                            Rprintf("no trans link, pats %d and %d\n", a_data -> patID, patseqID[k]);
                        }
                    }
                }
            }
        }
        a_data++;
    }
    return(gll);

}


// are patients pat1 and pat2 in the same group?
int samegroup(int pat1, int pat2, struct augdata *a_data) {
    int i=0, grp1=0, grp2=0, mark=0;
    for (i=0; i<num_patients; i++) {
        if (a_data -> patID == pat1) {
            grp1 = a_data -> group;
        }
        if (a_data -> patID == pat2) {
            grp2 = a_data -> group;
        }
        if (grp1!=0 && grp2!=0) {
            break;
        }
        a_data++;
    }
    if (grp1==0||grp2==0) {
        Rprintf("Group not found for %d, %d (%d, %d)\n", pat1, pat2, grp1, grp2);
    }
    if (grp1==grp2) {
        mark = 1;
    } else {
        mark = 0;
    }
    return(mark);
}

double loglikelihood(double *param, struct augdata *a_data, int dmat[][num_sequences], int *patseqID, int *resultsmatrix, int num_onadm, int model) {
    double ll=0.0;
    ll = likelihoodimport(param[0],num_onadm)+likelihoodtests(a_data,param[1], resultsmatrix)+likelihoodtrans(a_data, param[2]);
    if (model==1) {
        ll += likelihoodgen_IS(a_data, dmat, param[3], param[4], param[5], patseqID);
    } else if (model==2) {
        ll += likelihoodgen_TD(a_data, dmat, param[3], param[4], param[5], patseqID);
    }
    return(ll);
}

// number of transmission links between pat1 and pat2
int translinks(int pat1, int pat2, struct augdata *a_data) {
    int curpat = pat1;
    int j=0;
    int gen=1;
    if (pat1==pat2) {
        gen=0;
    } else {
        while (0) { // not reached end of chain
            if (a_data -> patID == curpat) { // got current patient?
                curpat = a_data -> psource; // set current patient to the source
                if (curpat == pat2) { // if the source is import
                    break; // end the process
                } else if (curpat == -1){ // got to start of chain before pat2
                    gen=0;
                    break;
                } else {
                    gen++; // otherwise, increase generation
                    a_data -= j; // and put pointer back to the start.
                    j = 0;
                }
            } else { // not found patient yet
                a_data++; // so continue combing data structure
                j++;
                if (j>=num_patients) { // didn't find anything
                    gen=0;
                    break;
                }
            }
        }
        if (gen==0) {
            a_data -= j;
            curpat = pat2;
            j=0;
            gen=1;
            while (0) { // not reached end of chain
                if (a_data -> patID == curpat) { // got current patient?
                    curpat = a_data -> psource; // set current patient to the source
                    if (curpat == pat1) { // if the source is import
                        break; // end the process
                    } else if (curpat == -1){ // got to start of chain before pat2
                        gen=0;
                        break;
                    } else {
                        gen++; // otherwise, increase generation
                        a_data -= j; // and put pointer back to the start.
                        j = 0;
                    }
                } else { // not found patient yet
                    a_data++; // so continue combing data structure
                    j++;
                    if (j>=num_patients) { // didn't find anything
                        gen=-1;
                        break;
                    }
                }
            }          
        }
    }
    return(gen);
}

// copy data structures
void copyaugdata(struct augdata *data1, struct augdata *data2) {
    int w;
    for (w=0; w<num_patients; w++) {
        data1 -> psource = data2 -> psource;
        data1 -> group = data2 -> group;
        data1 -> onadm = data2 -> onadm;
        data1 -> col_t = data2 -> col_t;
        data1 -> day_a = data2 -> day_a;
        data1 -> day_d = data2 -> day_d;
        data1 -> patID = data2 -> patID;
        data1 -> obspos = data2 -> obspos;
        data1++;
        data2++;
    }
}

// fetch test result for patient on given day
int result(int pID, int day, int *resultsmatrix) {
    // What test result for PID on this day?
    // resultsmatrix =1 if pos, 0 if neg, -1 if NA
    int i,res=0;
    for (i=0;i<num_patients;i++) {
        //first column of resultsmatrix is PID values
        if(resultsmatrix[i*(1+studyspan)]==pID) {
            res=resultsmatrix[(i*(1+studyspan))+day];
        }
    }
    if(res<-1||res>1) {
        Rprintf("WARNING: result(Pid=%d,day=%d)=%d\n", pID, day, res);
    }
    return(res);
}

// number of true positive results
int Num_TP(struct augdata *a_data, int *resultsmatrix) {
    int i,j, TP=0;
    for (i=0; i<num_patients; i++) {
        if (a_data -> col_t !=0) {
            for (j=a_data -> col_t;j<=a_data -> day_d; j++) {
                if (result(a_data -> patID, j, resultsmatrix)==1) TP++; // RESULT 1 == POSITIVE
            }
        }
        a_data++;
    }
    return(TP);
}

// number of false negative results
int Num_FN(struct augdata *a_data, int *resultsmatrix) {
    int i,j, FN=0;
    for (i=0; i<num_patients; i++) {
        if (a_data -> col_t !=0) {
            for (j=a_data -> col_t;j<=a_data -> day_d; j++) {
                if (result(a_data -> patID, j, resultsmatrix)==0) FN++; // RESULT 0 == NEGATIVE
            }
        }
        a_data++;
    }
    return(FN);
}

// log prob ratio: acq -> imp
double lpr_acq2imp(struct augdata *a_data_can, struct augdata *a_data_cur, double acq_const, double grp_const, int lastday, int j, int md) {
    double lp =0.0;
    lp = log(1-acq_const)-log(acq_const)-log(col_pop(a_data_cur-j, a_data_cur -> col_t))-log(lastday - a_data_cur -> day_a + 1);
    if (col_pop(a_data_cur-j, a_data_cur -> col_t)==0) {
        Rprintf("Warning: old colpop = %d\n", col_pop(a_data_cur-j, a_data_cur -> col_t));
    }
    if (md==1) {
        if (a_data_can -> group == a_data_cur -> patID) { // forming new group
            lp -= log(grp_const);
        } else {
            lp -= log(1-grp_const) + log(out_pop(a_data_cur-j, a_data_cur -> day_a));
            if (out_pop(a_data_cur-j, a_data_cur -> day_a)==0) {
                Rprintf("Warning (A->I): outpop = %d\n", out_pop(a_data_cur-j, a_data_cur -> day_a));
            }
        }
    }
    return(lp);
}

// log prob ratio: imp -> acq
double lpr_imp2acq(struct augdata *a_data_can, struct augdata *a_data_cur, double acq_const, double grp_const, int lastday, int j, int md) {
    double lp = 0.0;
    lp = log(acq_const)-log(1-acq_const)+log(col_pop(a_data_cur-j, a_data_can -> col_t)-1)+log(lastday - a_data_cur -> day_a + 1);
    if (col_pop(a_data_cur-j, a_data_can -> col_t)<=1) {
        Rprintf("Warning: old colpop = %d\n", col_pop(a_data_cur-j, a_data_can -> col_t));
    }
    if (md==1) {
        if (a_data_cur -> group == a_data_cur -> patID) { // previously group originator
            lp += log(grp_const);
        } else {
            lp += log(1-grp_const) + log(out_pop(a_data_cur-j, a_data_cur -> day_a));
            if (out_pop(a_data_cur-j, a_data_cur -> day_a)==0) {
                Rprintf("Warning (I -> A): outpop = %d\n", out_pop(a_data_cur-j, a_data_cur -> day_a));
            }
        }
    }
    return(lp);
    
}

// log prob ratio: acq -> acq
double lpr_acq2acq(struct augdata *a_data_can, struct augdata *a_data_cur, int j) {
    double lp = 0.0;
    if (a_data_cur -> col_t < a_data_can -> col_t) {
        lp = log(col_pop(a_data_cur-j, a_data_can -> col_t)-1)-log(col_pop(a_data_cur-j, a_data_can -> col_t));
    } else if (a_data_cur -> col_t > a_data_can -> col_t){
        lp = log(col_pop(a_data_cur-j, a_data_can -> col_t))-log(col_pop(a_data_cur-j, a_data_cur -> col_t));
    } else {
        lp = 0.0;
    }
    return(lp);
}

double lpr_imp2imp(struct augdata *a_data_can, struct augdata *a_data_cur, double grp_const, int j) {
    double lp=0.0;
    if (a_data_can -> group != a_data_cur -> group) { // new group
        if (a_data_can -> group == a_data_cur -> patID) { // creating own group
            lp = log(1-grp_const)-log(out_pop(a_data_cur-j, a_data_cur -> day_a)) - log(grp_const);
        } else if (a_data_cur -> group == a_data_cur -> patID) { // joining other group
            lp = log(grp_const)-log(1-grp_const) + log(out_pop(a_data_cur-j, a_data_cur -> day_a));
        }
    }
    return(lp);
}

double lpr_addacq(struct augdata *a_data_can, double acq_const, int col_add, int j) {
    double lp=0.0;
    lp = log(allneg-col_add)+log(a_data_can -> day_d - a_data_can -> day_a + 1)+log(col_pop(a_data_can-j, a_data_can -> col_t))-log(1-acq_const)-log(col_add+1);
    return(lp);
}

double lpr_addimp(struct augdata *a_data_can, double acq_const, double grp_const, int col_add, int md, int j) {
    double lp=0.0;
    lp = log(allneg-col_add)-log(col_add+1)-log(acq_const);
    if (md==1) {
        if (a_data_can -> patID == a_data_can -> group) { // own group
            lp -= log(grp_const);
        } else {
            lp += log(out_pop(a_data_can-j, a_data_can -> day_a)-1) - log(1-grp_const); // MINUS 1 HERE
            //lp += log(tot_groups(a_data_can-j)) - log(1-grp_const);
        }
    }
    return(lp);
}

double lpr_removeacq(struct augdata *a_data_cur, double acq_const, int col_add, int j) {
    double lp=0.0;
    lp = log(1-acq_const)+log(col_add)-log(a_data_cur -> day_d - a_data_cur -> day_a + 1)-log(col_pop(a_data_cur-j, a_data_cur -> col_t))-log(allneg-col_add-1);
    return(lp);
}

double lpr_removeimp(struct augdata *a_data_cur, double acq_const, double grp_const, int col_add, int md, int j) {
    double lp=0.0;
    lp = log(col_add)+log(acq_const)-log(allneg-col_add-1);
    if (md==1) {
        if (a_data_cur -> patID == a_data_cur -> group) { // own group
            lp += log(grp_const);
        } else {
            lp += log(1-grp_const) - log(out_pop(a_data_cur-j, a_data_cur -> day_a)-1); // MINUS 1 HERE
            //lp +=  log(1-grp_const) - log(tot_groups(a_data_cur-j));
        }
    }
    return(lp);
}

int offspringlimit(struct augdata *a_data, int pID) {
    int i=0, lim=99999;
    for (i=0;i<num_patients;i++) {
        if (a_data -> col_t != 0 && a_data -> patID != pID) {
            if (a_data -> psource == pID) {
                if (a_data -> col_t < lim) {
                    lim = a_data -> col_t - 1;
                }
            }
            if (a_data -> patID == pID) {
                if (a_data -> day_d < lim) {
                    lim = a_data -> day_d;
                }
            }
        }
        a_data++;
    }
    return(lim);
}


int directdescendant(struct augdata *a_data, int pID1, int pID2) {
    int curpat=pID1;
    int breaker = 0, mark=0, k=0, g=0;
    //Rprintf("%d", curpat);
    while (breaker==0) {
        if (g>10000) {
            Rprintf("DD function overload, PID %d %d\n", pID1, pID2);
            break;
        }
        if (a_data -> patID == curpat) {
            if (a_data -> psource == pID2) {
                mark = 1;
                breaker = 1;
            } else if (a_data -> psource == a_data -> patID) {
                mark = 0;
                breaker = 1;
            } else {
                curpat = a_data -> psource;
                //Rprintf("-%d", curpat);
                a_data -= k+1;
                k = -1;
            }
        }
        if (k==num_patients) {
            Rprintf("Didn't find psource %d (%d, %d))\n", curpat, pID1, pID2);
            a_data -= k;
            break;
        }
        k++;
        g++;
        a_data++;
    }
    //Rprintf("\n");
    return(mark);
}

int tot_groups(struct augdata *a_data) {
    int i=0, j=0;
    int gps=0;
    int groupvec[num_patients];
    for (i=0;i<num_patients;i++) {
        groupvec[i] = 0;
    }
    for (i=0;i<num_patients;i++) {
        if (a_data -> onadm == 1) {
            for (j=0; j<num_patients; j++) {
                if (groupvec[j] == a_data -> group) {
                    break;
                } else if (groupvec[j] == 0) {
                    groupvec[j] = a_data -> group;
                    break;
                }
            }
        }
        a_data++;
    }
    for (i=0;i<num_patients;i++) {
        if (groupvec[i]>0) {
            gps++;
        } else if (groupvec[i] == 0) {
            break;
        }
    }
    return(gps);
}

int num_imports(struct augdata *a_data) {
    int i=0;
    int imps=0;
    for (i=0;i<num_patients;i++) {
        if (a_data -> onadm == 1) {
            imps++;
        }
        a_data++;
    }
    return(imps);
}



