#include <stdio.h>
#include "time.h"
#include "math.h"
#include "string.h"
#include "unistd.h"
#include "getopt.h"
#include "global.h"
#include "util.h"
#include "string.h"
double rho = 0.9;
clock_t dur = 0;
int main (int argc, char *argv[])
{
    char Time_Buffer[64];
    double REP_NUM=0.0;
  double DELTA=0.0;
    double TOTAL_READ_COUNT=1; 
    double FDR_CUT=1.00;
  char *input_geneExpre = (char*)malloc(sizeof(char)*MAX_CHAR), 
    *outputf = (char*)malloc(sizeof(char)*MAX_CHAR),
    *read_counts_input = (char*)malloc(sizeof(char)*MAX_CHAR),
    *output_inf = (char*)malloc(sizeof(char)*MAX_CHAR);
    char str_line[MAX_LINE], *title_element_list[100];
    int nthread = 1, batch_size = 1, opt, row_num = 0, i = 0;

   //time recorder
    clock_t start,end,runtime;


    while ((opt = getopt(argc, argv, "hp:r:f:d:n:i:o:")) != -1) {
        switch (opt) {
        case 'h':
            printf("usage: -p: multiple processes number\n");
            printf("-r: min readCounts\n");
            printf("-f: fdr cutoff\n");
            printf("-d: delta cutoff\n");
            printf("-n: sample size\n");
            exit(1);
            break; 
        case 'p':
            nthread = atoi(optarg);
            break;           
        case 'r':
           TOTAL_READ_COUNT=atof(optarg);
           break;
        case 'f':
            FDR_CUT=atof(optarg);
            break;
        case 'd':
            DELTA=atof(optarg);
            break;
        case 'n':
           REP_NUM=atoi(optarg);
            break;
        case 'i':
            input_geneExpre = optarg;
            // printf("optind  %d\n",optind);
            // printf("argv_1   %s\n",argv[1]);
            // printf("argv_2   %s\n",argv[1]);
            optarg+=1+strlen(input_geneExpre); 
            read_counts_input = optarg;

            break;
        case 'o':
            outputf = optarg;
            break;
        default: /* '?' */
            fprintf(stderr, " 'h' for  Usage help: %s \n", argv[0]);
            exit(EXIT_FAILURE);
        }
        
    }
        printf("nthread  %d\nTOTAL_READ_COUNT  %f\nFDR_CUT  %f\nDELTA  %f\nREP_NUM  %f\n",nthread,TOTAL_READ_COUNT,FDR_CUT,DELTA,REP_NUM);
        printf("input_geneExpre   %s\n",input_geneExpre );
        printf("read_counts_input   %s\n",read_counts_input );
        printf("outputf   %s\n",outputf );


    //read and write
    FILE *ifp_geneExpre = NULL, *ifp_readCount = NULL,*ofp = NULL,*ofp_inf = NULL;
    if ((ifp_geneExpre = fopen(input_geneExpre, "r")) == NULL) {
        printf("Fail to open %s!\n", input_geneExpre);
        return -1;
    }
    if ((ifp_readCount = fopen(read_counts_input, "r")) == NULL) {
        printf("Fail to open %s!\n", read_counts_input);
        return -1;
    }    
    if ((ofp = fopen(outputf, "w")) == NULL) {
        printf("Fail to open %s!\n", outputf);
        return -1;
    }
    if ((ofp_inf = fopen("output_inf", "w")) == NULL) {
        printf("Fail to open output_inf.txt!\n");
        return -1;
    }
    start=clock();
    struct tm *tm_now;
    char *datetime;
    time(&start);
    tm_now = localtime(&start);
    datetime = asctime(tm_now);
    fprintf(ofp_inf,"start time: %s", datetime);

    //do something here

    end=clock();
    time(&end);
    tm_now = localtime(&end);
    datetime = asctime(tm_now);
    fprintf(ofp_inf,"end time: %s", datetime);


    int msec = runtime * 1000 / CLOCKS_PER_SEC;

    fprintf(ofp_inf,"Total Wallclock time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    fprintf(ofp_inf,"Wallclock time per thread taken %d seconds %d milliseconds\n", msec/1000/nthread, (msec%1000)/nthread);
    fprintf(ofp_inf,"Time for func(single thread): %d seconds %d milliseconds\n", msec/1000, msec%1000);
  return 0;
 }