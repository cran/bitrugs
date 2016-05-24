//
//  wgs_mcmc.c
//  
//
//  Created by Colin Worby on 20/05/2015.
//
//

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <time.h>
#include <string.h>
#include "wgsmcmc_head.h"

extern int num_patients;
extern int studyspan;
extern int num_sequences;
extern int TP;
extern int allneg;

void network_mcmc(int *n_patients, int *pID, int *day_adm, int *day_dis, int *col_times, int *psrc, int *grp, int *resultsmatrix, int *model,
                  int *n_seqs, int *seqID, int *distmat, int *iterations, int *augmoves, int *study_span, double *sigma,
                  double *priors, int *printno, char **path, double *outp, int *noaug) {
    
    int i=0, j=0, mc=0, mv=0, a=0, c=0, w=0, ol=0;
    double ll1=0.0, ll2=0.0, ll3=0.0, ll4=0.0;
    int src=0, host=0, colpool=0, outpool=0, lastday=0;
    num_patients = *n_patients;
    studyspan = *study_span;
    num_sequences = *n_seqs;
    int md = *model;
    
    Rprintf("Initialise sequences\n");
    
    int sqID_cur[num_sequences]; // vector of sequence IDs
    int sqID_can[num_sequences];
    int dmat_cur[num_sequences][num_sequences]; // distance matrix
    int dmat_can[num_sequences][num_sequences];
    for (i=0; i < num_sequences;i++) {
        for (j=0; j < num_sequences; j++) {
            dmat_cur[i][j] = distmat[num_sequences*i + j];
            //Rprintf("%d ", dmat_cur[i][j]);
            dmat_can[i][j] = dmat_cur[i][j];

        }
        //Rprintf("\n");
    }
    for (i=0;i<num_sequences;i++) {
        sqID_cur[i] = seqID[i];
        sqID_can[i] = sqID_cur[i];
    }
    
    int day_col[num_patients-1];
    int psource[num_patients-1];
    int group[num_patients-1];
    int onadm[num_patients-1];
    int op[num_patients-1]; // 1 if an added col time (no pos results)
    
    // parameter vectors
    double param_cur[6] = {0.05, 0.75, 0.005, 0.1, 0.05, 0.8};
    double param_can[6] = {0.05, 0.75, 0.005, 0.1, 0.05, 0.8};
    
    int num_acq = 0, num_posadm = 0, c_num_acq=0, c_num_posadm=0;
    int num_groups = 0, c_num_groups = 0;
    
    double cur_likelihood=0.0, can_likelihood=0.0;
    
    GetRNGstate();
    
    //initialize data
    
    
    for (i=0; i<num_patients; i++) {
        day_col[i] = col_times[i];
        psource[i] = psrc[i];
        group[i] = grp[i];
        if (psrc[i]==pID[i]) {
            onadm[i]=1;
        } else {
            onadm[i]=0;
        }
        if (onadm[i]==1) {
            num_posadm++;
        } else if (day_col[i]>0) {
            num_acq++;
        }
        if (group[i]==pID[i]) {
            num_groups++;
        }
        op[i] = 0;
        //Rprintf("%d\t%d\t%d\t%d\t%d\t%d\n", pID[i], day_adm[i], day_dis[i], day_col[i], onadm[i], group[i]);
    }

    c_num_acq = num_acq;
    c_num_posadm = num_posadm;
    c_num_groups = num_groups;

    allneg = num_patients-num_acq-num_posadm;
    
    struct augdata aug_cur[num_patients];               //create structure
    struct augdata aug_can[num_patients];

    
    for (a=0; a < num_patients; a++) {
        aug_cur[a].day_a = day_adm[a];
        aug_cur[a].day_d = day_dis[a];
        aug_cur[a].psource = psource[a];
        aug_cur[a].group = group[a];
        aug_cur[a].onadm = onadm[a];
        aug_cur[a].col_t = day_col[a];
        aug_cur[a].patID = pID[a];
        aug_cur[a].obspos = op[a];
        
        aug_can[a].day_a = day_adm[a];
        aug_can[a].day_d = day_dis[a];
        aug_can[a].psource = psource[a];
        aug_can[a].group = group[a];
        aug_can[a].onadm = onadm[a];
        aug_can[a].col_t = day_col[a];
        aug_can[a].patID = pID[a];
        aug_can[a].obspos = op[a];
        //Rprintf("%d %d %d %d %d %d\n", pID[a], day_adm[a], day_dis[a], psource[a], day_col[a], group[a]);
    }


    //Rprintf("No. pos on adm: %d\n", num_posadm);
    //Rprintf("No. acquired: %d\n", num_acq);
    //Rprintf("No. groups: %d\n", num_groups);
    
    struct augdata *acur_ptr;                //set pointers
    acur_ptr = &aug_cur[0];
    
    struct augdata *acan_ptr;
    acan_ptr = &aug_can[0];
    
    struct augdata *src_ptr;
    src_ptr = &aug_cur[0];


    //for (i=0;i<studyspan;i++) {
     //   Rprintf("%d ", col_pop(acur_ptr, i));
       // if (i%20==0) {
         //   Rprintf("\n");
       // }
    //}
    
    double logaccept_num=0.0, logaccept_denom=0.0;
    //Number of MCMC accepts
    int acc_lo = 200;
    int acc_hi = 400;
    //double log_prop_ratio = 0.0;

    //starting loglikelihood values
    TP = Num_TP(acur_ptr,resultsmatrix);
    Rprintf("\nTP: %d\n", TP);
    Rprintf("FN: %d\n", Num_FN(acur_ptr,resultsmatrix));
                            
    ll1 = likelihoodimport(param_cur[0],num_posadm);
    ll2 = likelihoodtests(acur_ptr,param_cur[1], resultsmatrix);
    ll3 = likelihoodtrans(acur_ptr, param_cur[2]);
    if (md==1) {
        ll4 = likelihoodgen_IS(acur_ptr, dmat_cur, param_cur[3], param_cur[4], param_cur[5], seqID);
    } else if (md==2) {
        ll4 = likelihoodgen_TD(acur_ptr, dmat_cur, param_cur[3], param_cur[4], param_cur[5], seqID);
    }
    
    cur_likelihood = loglikelihood(param_cur, acur_ptr, dmat_cur, seqID, resultsmatrix, num_posadm, md);
    Rprintf("ll1 %.4f\nll2 %.4f\nll3 %.4f\nll4 %.4f\nLL %.4f\n", ll1, ll2, ll3, ll4, cur_likelihood);
    can_likelihood = cur_likelihood;

    //open writeout files
    char var_file[100];
    //char aug_file[100];
    
    sprintf(var_file, "%sWGS_mcmc%d.txt", path[0], *printno);
    //sprintf(aug_file, "%sWGS_mcmc_detail%d.txt", path[0], *printno);
    Rprintf("File: %s\n", var_file);
    
    FILE *V0 = fopen(var_file, "w");
    
    //FILE *S1 = fopen(aug_file, "w");
    
    int updatedgroup = 0;
    //double probcarry = 0.0;
    //double probcarry2 = 0.0;
    //int mv;
    //int h=0;
    //int tempcolpop = 0;//, clst=0;
    //int tempID = 0;
    //int tempday = 0;
    //int clustered = 0;
    //int c_clustered = 0;
    int col_add = 0, col_add_c = 0;
    double choosemove = 0.0;
    double lpr = 0.0;
    double acq_const = 0.3;
    double grp_const = 0.4;

    int trans_acc[4] = {0, 0, 0, 0};
    //begin MCMC loop
    j=0;

    for (mc=1; mc<=*iterations; mc++) {                        //start MCMC loops
        
        ////fprintf(S1, "Move %d\n", mc);
        // Reset pointers
        acur_ptr = &aug_cur[0];
        acan_ptr = &aug_can[0];
        src_ptr = &aug_cur[0];
        
        if (mc%(int)pow(10.0,*outp)==0) {
           Rprintf(".");
        }
        // Check acceptance rates
        if (mc%(int)pow(10.0,*outp+1.0)==0) {
           Rprintf("\n%d\n", mc);
           Rprintf("Pos adm %d, num acq %d, no groups %d\n", num_posadm, num_acq, num_groups);
           if (*noaug>1) {
                Rprintf("Col add: %d, TP: %d, FN: %d\n", col_add, TP, Num_FN(acur_ptr, resultsmatrix));
           }

           ll1 = likelihoodimport(param_cur[0],num_posadm);
           ll2 = likelihoodtests(acur_ptr,param_cur[1], resultsmatrix);
           ll3 = likelihoodtrans(acur_ptr, param_cur[2]);
           if (md==1) {
                ll4 = likelihoodgen_IS(acur_ptr, dmat_cur, param_cur[3], param_cur[4], param_cur[5], seqID);
           } else if (md==2) {
                ll4 = likelihoodgen_TD(acur_ptr, dmat_cur, param_cur[3], param_cur[4], param_cur[5], seqID);
           }

           Rprintf("LL=%.3f %.3f %.3f %.3f\n", ll1, ll2, ll3, ll4);
           Rprintf("p\tz\tbeta\tgamma\tgamma_G\ttpar\n");
           for (i=0;i<6;i++) {
                Rprintf("%.5f\t", param_cur[i]);
            }
            Rprintf("\n");
           //Rprintf("col_add %d, fixed %d, clustered %d\n", col_add, add_offsp, clustered);
        }
        if (mc%1000==0) {
           for (i=0; i<=3; i++) {
               Rprintf("s%d=%.6f", i+1, sigma[i]);
               if (trans_acc[i] < acc_lo) {
                  sigma[i] *= 0.8;
                  Rprintf("-\t");	   
               } else if (trans_acc[i] > acc_hi) {
                  sigma[i] *= 1.25;
                  Rprintf("+\t");
               } else {
                  Rprintf("\t");
               }
               Rprintf("\tAcc: %d\n", trans_acc[i]);
               trans_acc[i] = 0;
           }
        }
        //fprintf(S1, "\nIteration %d:\np=%.3f, z=%.3f\n",mc,param_cur[0],param_cur[1]);
        //fprintf(S1,"Pos adm, num acq, num groups: %d\t%d\t%d\n", num_posadm, num_acq, num_groups);
        //fprintf(S1, "Col add: %d\n", col_add);
        
        // Parameter updates - Gibbs updates
        param_cur[0] = rbeta(num_posadm + priors[0], num_patients - num_posadm + priors[1]);
        param_cur[1] = rbeta(TP + priors[2], Num_FN(acur_ptr, resultsmatrix) + priors[3]);
        
        param_can[0] = param_cur[0];
        param_can[1] = param_cur[1];

        cur_likelihood = loglikelihood(param_cur, acur_ptr, dmat_cur, seqID, resultsmatrix, num_posadm, md);
        can_likelihood = cur_likelihood;
        
        // Parameter updates - MH updates
        
        // beta
        param_can[2] = param_cur[2] + rnorm(0.0,sigma[0]);
        if (param_can[2] > 0.0) {
        
           can_likelihood = loglikelihood(param_can, acur_ptr, dmat_cur, seqID, resultsmatrix, num_posadm, md);
           logaccept_num =  can_likelihood + dexp(param_can[2], priors[4],1);
           logaccept_denom = cur_likelihood + dexp(param_cur[2], priors[4],1);
           double u = runif(0.0, 1.0);
           if (log(u) < logaccept_num - logaccept_denom) {
              trans_acc[0]++;
              //fprintf(S1, "beta: %.6f -> %.6f acc\n",param_cur[2],param_can[2]);
              param_cur[2] = param_can[2];
              cur_likelihood = can_likelihood;
           } else {
              param_can[2] = param_cur[2];
              //fprintf(S1, "beta: %.6f -> %.6f rej\n",param_cur[2],param_can[2]);
           }
        } else {
           //fprintf(S1, "beta: %.6f -> %.6f rej\n",param_cur[2],param_can[2]);
        }
        
        // gamma
        param_can[3] = param_cur[3] + rnorm(0.0,sigma[1]);
        if (param_can[3] > 0.0) {
        
           can_likelihood = loglikelihood(param_can, acur_ptr, dmat_cur, seqID, resultsmatrix, num_posadm, md);
           logaccept_num =  can_likelihood + dexp(param_can[3], priors[5],1);
           logaccept_denom = cur_likelihood + dexp(param_cur[3], priors[5],1);
           //Rprintf("Cur: %.6f + %.6f", cur_likelihood, dexp(param_can[3], priors[5],1))
           double u = runif(0.0, 1.0);
           if (log(u) < logaccept_num - logaccept_denom) {
              trans_acc[1]++;
              //fprintf(S1, "gamma: %.6f -> %.6f acc\n",param_cur[3],param_can[3]);
              param_cur[3] = param_can[3];
              cur_likelihood = can_likelihood;
           } else {
              param_can[3] = param_cur[3];
              //fprintf(S1, "gamma: %.6f -> %.6f rej\n",param_cur[3],param_can[3]);
           }
        } else {
           //fprintf(S1, "gamma: %.6f -> %.6f rej\n",param_cur[3],param_can[3]);
        }
        
        // gamma_G
        param_can[4] = param_cur[4] + rnorm(0.0,sigma[2]);
        if (param_can[4] > 0.0) {
        
           can_likelihood = loglikelihood(param_can, acur_ptr, dmat_cur, seqID, resultsmatrix, num_posadm, md);
           logaccept_num =  can_likelihood + dexp(param_can[4], priors[6],1);
           logaccept_denom = cur_likelihood + dexp(param_cur[4], priors[6],1);
           double u = runif(0.0, 1.0);
           if (log(u) < logaccept_num - logaccept_denom) {
              trans_acc[2]++;
              //fprintf(S1, "gamma_G: %.6f -> %.6f acc\n",param_cur[4],param_can[4]);
              param_cur[4] = param_can[4];
              cur_likelihood = can_likelihood;
           } else {
              param_can[4] = param_cur[4];
              //fprintf(S1, "gamma_G: %.6f -> %.6f rej\n",param_cur[4],param_can[4]);
           }
        } else {
           //fprintf(S1, "gamma_G: %.6f -> %.6f rej\n",param_cur[4],param_can[4]);
        }
        
       
        // transpar
        param_can[5] = param_cur[5] + rnorm(0.0,sigma[3]);
        if (param_can[5] > 0.0) {
        
           can_likelihood = loglikelihood(param_can, acur_ptr, dmat_cur, seqID, resultsmatrix, num_posadm, md);
           logaccept_num =  can_likelihood + dexp(param_can[5], priors[7],1);
           logaccept_denom = cur_likelihood + dexp(param_cur[5], priors[7],1);
           double u = runif(0.0, 1.0);
           if (log(u) < logaccept_num - logaccept_denom) {
              trans_acc[3]++;
              //fprintf(S1, "tpar: %.6f -> %.6f acc\n",param_cur[5],param_can[5]);
              param_cur[5] = param_can[5];
              cur_likelihood = can_likelihood;
           } else {
              param_can[5] = param_cur[5];
              //fprintf(S1, "tpar: %.6f -> %.6f rej\n",param_cur[5],param_can[5]);
           }
        } else {
          //fprintf(S1, "tpar: %.6f -> %.6f rej\n",param_cur[5],param_can[5]);
        }
        // Sample augmented data
        if (*noaug>0) {
            for (mv=1;mv<=*augmoves;mv++) {
            
                if (*noaug == 1 || md==1) {
                    choosemove = 0.0; // only move
                } else {
                    if (col_add == 0) { // no cols added so only move, add
                        choosemove = runif(0.0, 2.0);
                    } else { // move, add, remove
                        choosemove = runif(0.0, 3.0);
                    }
                }
                // Move a colonisation time - pick individual to move
                c=0;
                acur_ptr = &aug_cur[0];
                acan_ptr = &aug_can[0]; 
                src_ptr = &aug_cur[0];
                
                if (choosemove < 1.0 && num_acq+num_posadm>0) {
                    src = (int) floor(runif(1.0, num_acq+num_posadm+1.0));
                    //fprintf(S1, "\nMove %dth col host\n", src);
                    for (i=0;i<num_patients;i++) {
                        if (acur_ptr -> col_t > 0) {
                            c++;
                        }
                        if (c==src) {
                            //fprintf(S1,"pID: %d adm: %d dis: %d onadm: %d\n", acur_ptr->patID, acur_ptr-> day_a, acur_ptr-> day_d, acur_ptr -> onadm);
                            //fprintf(S1,"col time: %d, psource: %d, group: %d\n", acur_ptr-> col_t, acur_ptr-> psource, acur_ptr-> group);
                            if (acur_ptr -> onadm == 1) {
                                if (runif(0.0,1.0)<acq_const && md==1) {
                                    //fprintf(S1, "Remain importation\n");
                                    //stay importation
                                    if (runif(0.0,1.0)<grp_const) {
                                        // Form new group
                                        //fprintf(S1, "New group\n");
                                        acan_ptr -> group = acur_ptr -> patID;
                                        c_num_groups = tot_groups(acan_ptr-i);
                                    } else {
                                        //Join existing group
                                        outpool = out_pop(acur_ptr-i, acan_ptr -> col_t);
                                        if (outpool>0) {
                                            w=0;
                                            host = (int) floor(runif(1.0, outpool+1.0));
                                            for (j=0;j<num_patients;j++) {
                                                if (src_ptr -> onadm == 1 && src_ptr -> day_a < acan_ptr -> col_t) {
                                                    w++;
                                                }
                                                if (w==host) {
                                                    // Found group
                                                    acan_ptr -> group = src_ptr -> group;
                                                    //fprintf(S1, "Join group %d (pID %d)\n", acan_ptr -> group, src_ptr -> patID);
                                                    c_num_groups = tot_groups(acan_ptr-i);
                                                    break;
                                                }
                                                if (j==(num_patients-1)) {
                                                    Rprintf("Error: %dth group for p%d not found (day %d, outpop %d)\n", host, acur_ptr -> patID, acan_ptr -> col_t, outpool);
                                                }
                                                src_ptr++;
                                            }
                                        } else {
                                            // If no previous imports
                                            //fprintf(S1, "No previous imports, day %d\n", acan_ptr -> col_t);
                                            acan_ptr -> group = acur_ptr -> patID;
                                            c_num_groups = tot_groups(acan_ptr-i);
                                        }
                                    }
                                    if (acur_ptr -> group == acan_ptr -> group) {
                                        lpr = 0.0;
                                    } else {
                                        lpr = lpr_imp2imp(acan_ptr, acur_ptr, grp_const, i);
                                    }
                                    //fprintf(S1, "LPR=%.4f\n", lpr);
                                    
                                } else {
                                    // reassign import to acquisition
                                    // 1. Look for first positive result
                                    lastday = acur_ptr -> day_d;
                                    for (j=acur_ptr -> day_a; j<=acur_ptr -> day_d; j++) {
                                        if (result(acur_ptr -> patID, j, resultsmatrix)==1) {
                                            //fprintf(S1, "Last day=%d\n", lastday);
                                            lastday = j;
                                            break;
                                        }
                                    }
                                    // What if there this host is the source for other individuals?
                                    ol = offspringlimit(acur_ptr-i, acur_ptr -> patID);
                                    if (ol < lastday) {
                                        lastday = ol;
                                    }
                                    if (lastday<acur_ptr -> day_a) {
                                        lastday = acur_ptr -> day_a;
                                    }
                                    // 2. Propose day of colonisation
                                    acan_ptr -> col_t = (int) floor(runif(acur_ptr -> day_a, lastday+1.0));
                                    //fprintf(S1, "Propose colonisation on day %d\n", acan_ptr -> col_t);
                                    colpool = col_pop(acur_ptr-i, acan_ptr -> col_t)-1; // Don't include self in pool of sources.

                                    if (colpool>0) {
                                        w=0;
                                        host = (int) floor(runif(1.0, colpool+1.0));
                                        for (j=0;j<num_patients;j++) {
                                            if (src_ptr -> col_t!=0 && src_ptr -> day_d >= acan_ptr -> col_t && src_ptr -> patID != acur_ptr -> patID) {
                                                if (src_ptr -> col_t < acan_ptr -> col_t) {
                                                    w++;
                                                } else if (src_ptr -> col_t==acan_ptr -> col_t && src_ptr -> onadm == 1) {
                                                    w++;
                                                }
                                                if (w==host) {
                                                    // 3. Found host, take source, group
                                                    //fprintf(S1, "Proposed host: ID %d, a %d d %d c %d p %d oa %d g %d\n", src_ptr -> patID, src_ptr -> day_a, src_ptr -> day_d, src_ptr -> col_t, src_ptr -> psource, src_ptr -> onadm, src_ptr -> group);
                                                    acan_ptr -> psource = src_ptr -> patID;
                                                    acan_ptr -> onadm = 0;
                                                    acan_ptr -> group = src_ptr -> group;
                                                    c_num_acq++;
                                                    c_num_posadm--;
                                                    c_num_groups = tot_groups(acan_ptr-i);
                                                    break;
                                                }
                                            }
                                            if (j==(num_patients-1)) {
                                                Rprintf("Error: %dth source for p%d not found (day %d, colpop %d)\n", host, acan_ptr -> col_t, acur_ptr -> patID, colpool);
                                            }
                                            src_ptr++;
                                        }
                                    } else {
                                        // If no other colonised host on proposed day
                                        //fprintf(S1, "No hosts on day %d (colpop=%d)\n", acan_ptr -> col_t, col_pop(acur_ptr-i, acan_ptr -> col_t));
                                        acan_ptr -> col_t = acur_ptr -> col_t;
                                        acan_ptr -> onadm = 1;
                                        acan_ptr -> group = acur_ptr -> group;
                                    }
                                    if (acan_ptr -> onadm != 1) {
                                        lpr = lpr_imp2acq(acan_ptr, acur_ptr, acq_const, grp_const, lastday, i, md);
                                    } else {
                                        lpr = 0.0;
                                    }
                                    //fprintf(S1, "LPR=%.4f\n", lpr);
                                }
                                    
                            } else if (runif(0.0,1.0)<acq_const) {
                                // reassign acquisition to an importation
                                //fprintf(S1, "Propose importation\n");
                                acan_ptr -> onadm = 1;
                                acan_ptr -> psource = acur_ptr -> patID;
                                acan_ptr -> col_t = acur_ptr -> day_a;
                                // Find possible days on which acquisition could have occured
                                lastday = acur_ptr -> day_d;
                                for (j=acur_ptr -> day_a; j<=acur_ptr -> day_d; j++) {
                                    if (result(acur_ptr -> patID, j, resultsmatrix)==1) {
                                        lastday = j;
                                        break;
                                    }
                                }
                                // To deal with potential offspring
                                ol = offspringlimit(acur_ptr-i, acur_ptr -> patID);
                                if (ol < lastday) {
                                    lastday = ol;
                                }
                                if (lastday<acur_ptr -> day_a) {
                                    lastday = acur_ptr -> day_a;
                                }
                                if (md==2) {
                                    acan_ptr -> group = acur_ptr -> patID;
                                    c_num_groups++;
                                } else {
                                    
                                    if (runif(0.0,1.0)<grp_const && md==1) {
                                        // Form new group
                                        //fprintf(S1, "New group\n");
                                        acan_ptr -> group = acur_ptr -> patID;
                                        c_num_groups = tot_groups(acan_ptr-i);
                                        
                                    } else {
                                        //Join existing group
                                        outpool = out_pop(acur_ptr-i, acan_ptr -> col_t);
                                        if (outpool>0) {
                                            w=0;
                                            host = (int) floor(runif(1.0, outpool+1.0));
                                            for (j=0;j<num_patients;j++) {
                                                if (src_ptr -> onadm == 1 && src_ptr -> day_a < acan_ptr -> col_t) {
                                                    w++;
                                                }
                                                if (w==host) {
                                                    // Found group
                                                    acan_ptr -> group = src_ptr -> group;
                                                    //fprintf(S1, "Join group %d (pID %d)\n", acan_ptr -> group, src_ptr -> patID);
                                                    break;
                                                }
                                                if (j==(num_patients-1)) {
                                                    Rprintf("Error: %dth group for p%d not found (day %d, outpop %d)\n", host, acur_ptr -> patID, acan_ptr -> col_t, outpool);
                                                }
                                                src_ptr++;
                                            }
                                        } else {
                                            // If no previous imports
                                            //fprintf(S1, "No previous imports, day %d\n", acan_ptr -> col_t);
                                            acan_ptr -> group = acur_ptr -> patID;
                                            c_num_groups = tot_groups(acan_ptr-i);
                                        }
                                    }
                                }
                                c_num_posadm++;
                                c_num_acq--;
                                lpr = lpr_acq2imp(acan_ptr, acur_ptr, acq_const, grp_const, lastday, i, md);
                                //fprintf(S1, "LPR=%.4f\n", lpr);
                                
                            } else {
                                // acquisition remains acquisition
                                // 1. Look for first positive result
                                lastday = acur_ptr -> day_d;
                                for (j=acur_ptr -> day_a; j<=acur_ptr -> day_d; j++) {
                                    if (result(acur_ptr -> patID, j, resultsmatrix)==1) {
                                        lastday = j;
                                        break;
                                    }
                                }
                                // Account for potential offspring limits
                                ol = offspringlimit(acur_ptr-i, acur_ptr -> patID);
                                if (ol < lastday) {
                                    lastday = ol;
                                }
                                if (lastday<acur_ptr -> day_a) {
                                    lastday = acur_ptr -> day_a;
                                }
                                // 2. Propose day of colonisation
                                acan_ptr -> col_t = (int) floor(runif(acur_ptr -> day_a, lastday+1.0));
                                //fprintf(S1, "Propose colonisation on day %d\n", acan_ptr -> col_t);
                                colpool = col_pop(acur_ptr-i, acan_ptr -> col_t);
                                if (acan_ptr -> col_t > acur_ptr -> col_t) {
                                    colpool--; // Don't include self in pool of sources.
                                }
                                if (colpool>0) {
                                    w=0;
                                    host = (int) floor(runif(1.0, colpool+1.0));
                                    for (j=0;j<num_patients;j++) {
                                        if (src_ptr -> col_t!=0 && src_ptr -> day_d >= acan_ptr -> col_t) {
                                            if (src_ptr -> col_t < acan_ptr -> col_t && src_ptr -> patID != acur_ptr -> patID) {
                                                w++;
                                            } else if (src_ptr -> col_t==acan_ptr -> col_t && src_ptr -> onadm == 1 && src_ptr -> patID != acur_ptr -> patID) {
                                                w++;
                                            }
                                            if (w==host) {
                                                // 3. Found host, take source, group
                                                //fprintf(S1, "Proposed host: ID %d, a %d d %d c %d p %d oa %d g %d\n", src_ptr -> patID, src_ptr -> day_a, src_ptr -> day_d, src_ptr -> col_t, src_ptr -> psource, src_ptr -> onadm, src_ptr -> group);
                                                acan_ptr -> psource = src_ptr -> patID;
                                                acan_ptr -> onadm = 0;
                                                acan_ptr -> group = src_ptr -> group;
                                                
                                                break;
                                            }
                                        }
                                        if (j==(num_patients-1)) {
                                            Rprintf("Error: %dth source for p%d not found (day %d, colpop %d)\n", host, acur_ptr -> patID, colpool);
                                        }
                                        src_ptr++;
                                    }
                                } else {
                                    // If no other colonised host on proposed day
                                    //fprintf(S1, "No hosts on day %d\n", acan_ptr -> col_t);
                                    acan_ptr -> col_t = acur_ptr -> col_t;
                                    acan_ptr -> onadm = 0;
                                }
                                lpr = lpr_acq2acq(acan_ptr, acur_ptr, i);
                                //fprintf(S1, "LPR=%.4f\n", lpr);
                            }
                            
                            break;
                        }
                        if (i==(num_patients-1)) {
                            Rprintf("Error: %dth col host to move not found\n", c);
                        }
                        acur_ptr++;
                        acan_ptr++;
                    }
                    
                } else if (choosemove < 2.0) {
                    // Add colonisation time
                    //fprintf(S1, "\nAdd col. time\n");
                    // 1. Pick a negative host
                    src = (int) floor(runif(1.0, num_patients-num_acq-num_posadm+1.0));
                    //fprintf(S1, "Pick %dth negative host\n", src);
                    for (i=0;i<num_patients;i++) {
                        if (acur_ptr -> col_t == 0) {
                            c++;
                            if (c==src) {
                                //fprintf(S1, "Pick: ID %d, a %d d %d c %d oa %d p %d g %d\n", acur_ptr -> patID, acur_ptr -> day_a, acur_ptr -> day_d, acur_ptr -> col_t, acur_ptr -> onadm, acur_ptr -> psource, acur_ptr -> group);
                                // 2. Choose importation / acquisition
                                if (runif(0.0,1.0)<acq_const) {
                                    //fprintf(S1, "Add importation\n");
                                    acan_ptr -> onadm = 1;
                                    acan_ptr -> col_t = acur_ptr -> day_a;
                                    acan_ptr -> psource = acur_ptr -> patID;
                                    acan_ptr -> obspos = 1;
                                    // New group?
                                    if (runif(0.0,1.0)< grp_const || md == 2) { //grpconst
                                        // Form new group
                                        //fprintf(S1, "New group\n");
                                        acan_ptr -> group = acur_ptr -> patID;
                                        c_num_groups = tot_groups(acan_ptr-i);
                                        c_num_posadm++;
                                        col_add_c++;
                                        lpr = lpr_addimp(acan_ptr, acq_const, grp_const, col_add, md, i);
                                    } else {
                                        // Add existing group
                                        outpool = out_pop(acur_ptr-i, acan_ptr -> col_t);
                                        //fprintf(S1, "Outpop: %d\n", outpool);
                                        if (outpool>0) {
                                            w=0;
                                            host = (int) floor(runif(1.0, outpool+1.0));
                                            for (j=0;j<num_patients;j++) {
                                                if (src_ptr -> onadm == 1 && src_ptr -> day_a < acan_ptr -> col_t) {
                                                    w++;
                                                }
                                                if (w==host) {
                                                    // Found group
                                                    acan_ptr -> group = src_ptr -> group;
                                                    //fprintf(S1, "Join group %d (pID %d)\n", acan_ptr -> group, src_ptr -> patID);
                                                    c_num_posadm++;
                                                    col_add_c++;
                                                    lpr = lpr_addimp(acan_ptr, acq_const, grp_const, col_add, md, i);
                                                    break;
                                                }
                                                if (j==(num_patients-1)) {
                                                    Rprintf("Error (add imp): %dth group for p%d not found (day %d, outpop %d)\n", host, acur_ptr -> patID, acan_ptr -> col_t,  outpool);
                                                }
                                                src_ptr++;
                                            }
                                        } else {
                                            // If no previous imports
                                            // make own group
                                            //fprintf(S1, "No previous imports, day %d\n", acan_ptr -> col_t);
                                            acan_ptr -> group = acur_ptr -> patID;
                                            c_num_groups = tot_groups(acan_ptr-i);
                                            c_num_posadm++;
                                            col_add_c++;
                                            lpr = lpr_addimp(acan_ptr, acq_const, grp_const, col_add, md, i);
                                        }

                                    }
                                    //fprintf(S1, "LPR=%.4f\n", lpr);
                                
                                } else {
                                    //fprintf(S1, "Add acquisition\n");
                                    // Propose day of colonisation
                                    acan_ptr -> col_t = (int) floor(runif(acur_ptr -> day_a, acur_ptr -> day_d+1.0));
                                    colpool = col_pop(acur_ptr-i, acan_ptr -> col_t);
                                    //fprintf(S1, "Propose colonisation on day %d, colpop %d\n", acan_ptr -> col_t, colpool);

                                    if (colpool>0) {
                                        w=0;
                                        host = (int) floor(runif(1.0, colpool+1.0));
                                        for (j=0;j<num_patients;j++) {
                                            if (src_ptr -> col_t!=0 && src_ptr -> day_d >= acan_ptr -> col_t && src_ptr -> patID != acur_ptr -> patID) {
                                                if (src_ptr -> col_t < acan_ptr -> col_t) {
                                                    w++;
                                                } else if (src_ptr -> col_t==acan_ptr -> col_t && src_ptr -> onadm == 1) {
                                                    w++;
                                                }
                                                if (w==host) {
                                                    // 3. Found host, take source, group
                                                    //fprintf(S1, "Proposed host: ID %d, a %d d %d c %d p %d oa %d g %d\n", src_ptr -> patID, src_ptr -> day_a, src_ptr -> day_d, src_ptr -> col_t, src_ptr -> psource, src_ptr -> onadm, src_ptr -> group);
                                                    acan_ptr -> psource = src_ptr -> patID;
                                                    acan_ptr -> onadm = 0;
                                                    acan_ptr -> group = src_ptr -> group;
                                                    acan_ptr -> obspos = 1;
                                                    c_num_acq++;
                                                    col_add_c++;
                                                    lpr = lpr_addacq(acan_ptr, acq_const, col_add, i);
                                                    break;
                                                }
                                            }
                                            if (j==(num_patients-1)) {
                                                Rprintf("Error (add acq): %dth source for p%d not found (day %d, colpop %d)\n", host, acan_ptr -> col_t, acur_ptr -> patID, colpool);
                                            }
                                            src_ptr++;
                                        }
                                    } else {
                                        // If no other colonised host on proposed day
                                        //fprintf(S1, "No hosts on day %d (colpop=%d)\n", acan_ptr -> col_t, col_pop(acur_ptr-i, acan_ptr -> col_t));
                                        acan_ptr -> col_t = acur_ptr -> col_t;
                                        acan_ptr -> onadm = acur_ptr -> onadm;
                                        acan_ptr -> group = acur_ptr -> group;
                                        acan_ptr -> obspos = acur_ptr -> obspos;
                                        lpr = 0.0;
                                    }

                                    //fprintf(S1, "LPR=%.4f\n", lpr);

                                }
                                break;
                            }
                        }
                        if (i==(num_patients-1)) {
                            Rprintf("Error: %dth col host to move not found\n", c);
                        }
                        acur_ptr++;
                        acan_ptr++;
                    }
                    
                } else {
                    // Remove colonisation time
                    //fprintf(S1, "\nRemove col. time\n");
                    src = (int) floor(runif(1.0, col_add+1.0));
                    //fprintf(S1, "Pick %dth added host\n", src);
                    for (i=0;i<num_patients;i++) {
                        if (acur_ptr -> obspos == 1) {
                            c++;
                        }
                        if (c==src) { // found person to remove
                            //fprintf(S1, "Pick: ID %d, a %d d %d c %d oa %d p %d g %d\n", acur_ptr -> patID, acur_ptr -> day_a, acur_ptr -> day_d, acur_ptr -> col_t, acur_ptr -> onadm, acur_ptr -> psource, acur_ptr -> group);
                            if (offspringlimit(acur_ptr-i, acur_ptr -> patID)==99999) { // no offspring
                            
                                if (acur_ptr -> onadm == 1) {
                                    //fprintf(S1, "Propose removing import\n");
                                    // what if group originator? Shouldn't matter - sorted later (reconcile groupings)
                                    outpool = out_pop(acur_ptr-i, acan_ptr -> col_t);
                                    //fprintf(S1, "Outpop: %d\n", outpool);
                                    acan_ptr -> col_t = 0;
                                    acan_ptr -> psource = 0;
                                    acan_ptr -> group = 0;
                                    acan_ptr -> obspos = 0;
                                    acan_ptr -> onadm = 0;
                                    col_add_c--;
                                    c_num_posadm--;
                                    c_num_groups = tot_groups(acan_ptr-i);
                                    lpr = lpr_removeimp(acur_ptr, acq_const, grp_const, col_add, md, i);
                                    for (j=0;j<num_patients;j++) {
                                        if (src_ptr -> group == acur_ptr -> group) {
                                            //fprintf(S1, "Same group %d (p %d)\n", src_ptr -> patID, src_ptr-> psource);
                                        }
                                        src_ptr++;
                                    }
                                    
                                } else {
                                    //fprintf(S1, "Propose removing acquisition\n");
                                    acan_ptr -> col_t = 0;
                                    acan_ptr -> psource = 0;
                                    acan_ptr -> group = 0;
                                    acan_ptr -> obspos = 0;
                                    acan_ptr -> onadm = 0;
                                    col_add_c--;
                                    c_num_acq--;
                                    c_num_groups = tot_groups(acan_ptr-i);
                                    lpr = lpr_removeacq(acur_ptr, acq_const, col_add, i);
                                    for (j=0;j<num_patients;j++) {
                                        if (src_ptr -> group == acur_ptr -> group) {
                                            //fprintf(S1, "Same group %d (p %d)\n", src_ptr -> patID, src_ptr-> psource);
                                        }
                                        src_ptr++;
                                    }
                                }
                            } else {
                                //fprintf(S1, "Host %d has offspring\n", acur_ptr -> patID);
                                lpr = 0.0;
                            }
                            //fprintf(S1, "LPR: %.4f\n", lpr);
                            break;
                        }
                        if (i==num_patients-1) {
                            Rprintf("Error (remove): %dth col host to remove not found\n", src);
                        }
                        acur_ptr++;
                        acan_ptr++;
                    }
                    
                }
                //fprintf(S1, "Cur: ID %d, a %d d %d c %d oa %d p %d g %d\n", acur_ptr -> patID, acur_ptr -> day_a, acur_ptr -> day_d, acur_ptr -> col_t, acur_ptr -> onadm, acur_ptr -> psource, acur_ptr -> group);
                //fprintf(S1, "Can: ID %d, a %d d %d c %d oa %d p %d g %d\n", acan_ptr -> patID, acan_ptr -> day_a, acan_ptr -> day_d, acan_ptr -> col_t, acan_ptr -> onadm, acan_ptr -> psource, acan_ptr -> group);
                //fprintf(S1, "New col pop: %d\n", col_pop(acan_ptr-i, acan_ptr -> col_t));
                
                // Reconcile groupings
                updatedgroup = acan_ptr -> group;
                
                // group originator (re)moved?
                if (acur_ptr -> onadm==1 && acan_ptr -> onadm==0) {
                    acan_ptr = &aug_can[0];
                    int earliestimport=99999;
                    int potgroup=0;
                    for (j=0;j<num_patients;j++) {
                        if (acan_ptr -> group == acur_ptr -> group && acan_ptr -> onadm == 1 && acan_ptr -> patID != acur_ptr -> patID) { // found candidate import
                            if (acan_ptr -> day_a < earliestimport) {
                                earliestimport = acan_ptr -> day_a;
                                potgroup = acan_ptr -> patID;
                            }
                        }
                        acan_ptr++;
                    }
                    acan_ptr = &aug_can[0];
                    if (earliestimport<99999 && acan_ptr -> group != potgroup) {
                        for (j=0;j<num_patients;j++) {
                            if (acan_ptr -> group == acur_ptr -> group) {
                                //fprintf(S1, "Update group (orig) for %d (psource %d) from %d to %d\n", acan_ptr -> patID, acan_ptr -> psource, acan_ptr -> group, potgroup);
                                acan_ptr -> group = potgroup;
                            }
                            acan_ptr++;
                        }
                    }
                } else if (acur_ptr -> group != updatedgroup && choosemove < 1.0) { // groups need reconciling after a move
                    acan_ptr = &aug_can[0];
                    //fprintf(S1, "updatedgroup %d\n", updatedgroup);
                    for (j=0;j<num_patients;j++) {
                        if (acan_ptr -> group == acur_ptr -> group && acan_ptr -> patID != acur_ptr -> patID) { // if belongs to the old group
                            if (directdescendant(acur_ptr-i, acan_ptr -> patID, acur_ptr -> patID)==1) {
                                //fprintf(S1, "Update group (trans) for %d (psource %d) from %d to %d\n", acan_ptr -> patID, acan_ptr -> psource, acan_ptr -> group, updatedgroup);
                                acan_ptr -> group = updatedgroup;
                            }
                        }
                        acan_ptr++;
                    }
                }

                
                acur_ptr = &aug_cur[0];
                acan_ptr = &aug_can[0];
                src_ptr = &aug_cur[0];
                
                // Evaluate and accept/reject proposed update
            
                can_likelihood = loglikelihood(param_cur, acan_ptr, dmat_cur, seqID, resultsmatrix, c_num_posadm, md);
                //fprintf(S1, "Cur LL: %.4f, Can LL: %.4f\n", cur_likelihood, can_likelihood);

//                ll1 = likelihoodimport(param_cur[0],num_posadm);
//                ll2 = likelihoodtests(acur_ptr,param_cur[1], resultsmatrix);
//                ll3 = likelihoodtrans(acur_ptr, param_cur[2]);
//                if (md==1) {
//                    ll4 = likelihoodgen_IS(acur_ptr, dmat_cur, param_cur[3], param_cur[4], param_cur[5], seqID);
//                } else if (md==2) {
//                    ll4 = likelihoodgen_TD(acur_ptr, dmat_cur, param_cur[3], param_cur[4], param_cur[5], seqID);
//                }
//                //fprintf(S1, "Cur LL1 %.4f, LL2 %.4f, LL3 %.4f, LL4 %.4f = %.4f (%.4f)\n", ll1, ll2, ll3, ll4, ll1+ll2+ll3+ll4, cur_likelihood);
//                
//                ll1 = likelihoodimport(param_cur[0],c_num_posadm);
//                ll2 = likelihoodtests(acan_ptr,param_cur[1], resultsmatrix);
//                ll3 = likelihoodtrans(acan_ptr, param_cur[2]);
//                if (md==1) {
//                    ll4 = likelihoodgen_IS(acan_ptr, dmat_cur, param_cur[3], param_cur[4], param_cur[5], seqID);
//                } else if (md==2) {
//                    ll4 = likelihoodgen_TD(acan_ptr, dmat_cur, param_cur[3], param_cur[4], param_cur[5], seqID);
//                }
//                //fprintf(S1, "Can LL1 %.4f, LL2 %.4f, LL3 %.4f, LL4 %.4f = %.4f\n", ll1, ll2, ll3, ll4, ll1+ll2+ll3+ll4);
//                can_likelihood = ll1+ll2+ll3+ll4;
                
                ////fprintf(S1, "Accept with P=%.4f\n", exp(can_likelihood - cur_likelihood));
                double ran_acc = runif(0.0,1.0);
                if (log(ran_acc)< can_likelihood - cur_likelihood + lpr) {
                    
                    cur_likelihood = can_likelihood;
                    copyaugdata(acur_ptr, acan_ptr);
                    num_acq = c_num_acq;
                    num_posadm = c_num_posadm;
                    num_groups = c_num_groups;
                    col_add = col_add_c;
                    //fprintf(S1, "Accept\n\n");
                    
                } else {
                    copyaugdata(acan_ptr, acur_ptr);
                    c_num_acq = num_acq;
                    c_num_posadm = num_posadm;
                    c_num_groups =num_groups;
                    col_add_c = col_add;
                    //fprintf(S1, "Reject\n\n");
                }
                //fprintf(S1,"Num_posadm: %d Num_acq: %d Num groups: %d\n", num_posadm, num_acq, num_groups);
                //fprintf(S1,"TP:%d\tFN:%d\n", TP, Num_FN(acur_ptr, resultsmatrix));
            } // End augdata moves
        }
        
        // Output parameters
        for (i=0;i<6;i++) {
           fprintf(V0, "%.8f\t", param_cur[i]);
        }
        fprintf(V0,"%d\t%d\t%d\t%.6f", num_posadm, num_acq, num_groups, cur_likelihood);
        for (i=0;i<num_patients;i++) {
            fprintf(V0, "\t%d", acur_ptr -> psource);
            acur_ptr++;
        }
        acur_ptr = &aug_cur[0];
        for (i=0;i<num_patients;i++) {
            fprintf(V0, "\t%d", acur_ptr -> group);
            acur_ptr++;
        }
        acur_ptr = &aug_cur[0];
        fprintf(V0, "\n");

    }

    
     //fclose(S1);
     fclose(V0);
     PutRNGstate();
    
}
    
