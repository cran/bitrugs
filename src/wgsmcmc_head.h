int num_patients;
int studyspan;
int TP;
int num_sequences;
int allneg;


struct augdata {
       int day_a;
       int day_d;
       int psource;
       int group;
       int onadm;
       int col_t;
       int patID;
       int obspos;
};

// total infected days avoiding colonisation
int avoideddays(struct augdata*, int from, int to);

// colonised population on a particular day
int col_pop(struct augdata*, int day);

// Number of importation groups
int out_pop(struct augdata*, int day);

// test result loglikelihood contribution
double likelihoodimport(double p, int num_onadm);

// test observation loglikelihood contribution
double likelihoodtests(struct augdata*, double z, int *resultsmatrix);

// transmission dynamic loglikelihood contribution
double likelihoodtrans(struct augdata*, double beta);

// genetic distance loglikelihood contribution: importation structure model

double likelihoodgen_IS(struct augdata*, int dmat[][num_sequences], double gamma, double gamma_G, double clust_par, int *patseqID);

// genetic distance loglikelihood contribution: transmission diversity model

double likelihoodgen_TD(struct augdata*, int dmat[][num_sequences], double gamma, double gamma_G, double trans_par, int *patseqID);

// are patients pat1 and pat2 in the same group?
int samegroup(int pat1, int pat2, struct augdata*);

// full likelihood
double loglikelihood(double *param, struct augdata*, int dmat[][num_sequences], int *patseqID, int *resultsmatrix, int num_onadm, int model);

// number of transmission links between pat1 and pat2
int translinks(int pat1, int pat2, struct augdata*);
// copy data structures
void copyaugdata(struct augdata*, struct augdata *data2);

// fetch test result for patient on given day
int result(int pID, int day, int *resultsmatrix);

// number of true positive results
int Num_TP(struct augdata*, int *resultsmatrix);

// number of false negative results
int Num_FN(struct augdata*, int *resultsmatrix);

double lpr_acq2imp(struct augdata*, struct augdata *data2, double acq_const, double grp_const, int lastday, int j, int md);

double lpr_imp2acq(struct augdata*, struct augdata *data2, double acq_const, double grp_const, int lastday, int j, int md);

double lpr_acq2acq(struct augdata*, struct augdata *data2, int j);

double lpr_imp2imp(struct augdata*, struct augdata *data2, double grp_const, int j);

int offspringlimit(struct augdata*, int pID);

int directdescendant(struct augdata*, int pID1, int pID2);

int tot_groups(struct augdata *a_data);

int num_imports(struct augdata *a_data);

double lpr_addacq(struct augdata*, double acq_const, int col_add, int j);

double lpr_addimp(struct augdata*, double acq_const, double grp_const, int col_add, int md, int j);

double lpr_removeacq(struct augdata*, double acq_const, int col_add, int j);

double lpr_removeimp(struct augdata*, double acq_const, double grp_const, int col_add, int md, int j);
