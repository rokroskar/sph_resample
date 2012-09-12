#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <fenv.h>
#include "kd.h"
#include "smooth1.h"
//#include "tipsydefs.h"

#define PRUNE_STEPS		5
#define MICRO_STEP		0.1

struct dump {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
} ;


void usage(void)
{
	fprintf(stderr,"BASIC USAGE:\n");
	fprintf(stderr,"sph_rsmpl -nSplitting <N> -radius <R> -s <nSmooth> -I <filename> \n");
	fprintf(stderr,"	 reads TIPSY standard or checkpoint file from stdin\n");
	fprintf(stderr,"-diag  : print diagnostic output\n");
        
        fprintf(stderr,"\nOTHER OPTIONS:\n");
	fprintf(stderr,"-smbh  : this is a run with SMBHs (applies to determining com)\n");
        fprintf(stderr,"-ckpnt : the input file is a checkpoint (default is tipsy standard file)\n");
        fprintf(stderr,"-Lc    : specifies max. angular momentum instead of radius for particle selection\n");
        
	fprintf(stderr,"\nOUTPUT filenames are \n");
	fprintf(stderr,"\t INPUT_FILE.subsample.tipsy and \n");
	fprintf(stderr,"\t INPUT_FILE.subsample.chk\n");
	exit(1);
	}

int main(int argc,char **argv)
{
        KD kdg, kds, kdd;
	SMX smx;
	/*
	 ** Input argument variables and control.
	 */
	int nSmooth,bCkpnt=0, bSmbh;
	int nMaxMembers;
	int bPeriodic;
	int bStandard;
	int bAux;
	float fPeriod[3],fCenter[3];
	char achName[256];
	/*
	 ** Working variables for SKID mainline.
	 */
	float fStep;
	int nBucket,i,j,nActive,nIttr;
	int sec1,usec1;
	char achFile[256];
	int iExt;
	int nScat;
	int bOutDiag;

	int nSplitting;
	CHK_HEADER chk_header;
	struct dump tipsy_header;
	float radius, Lc = 0.0, alpha = 1.0;

	feenableexcept(FE_OVERFLOW | FE_DIVBYZERO | FE_INVALID);

	printf("SPH_RSMPL: Rok Roskar 2009/2010 - based entirely on SKID by Joachim Stadel\n");
	/*
	 ** Bucket size set to 16, user cannot affect this!
	 */
   	nBucket = 16;

	/*
	 * Default tipsy native format.
	 */
	bStandard = 0;

	/*
	 * Default number of nearest neighbors
	 */
	nSmooth = 32;

	/* 
	** Assume it is *not* an SMBH run - this applies to determination of COFM
	*/
	bSmbh = 0;
	bPeriodic = 0;
	bAux = 0;
	/*
	 ** Default output parameters.
	 */
	strcpy(achName,"split");
	
	/*
	 ** Now get the command line arguments!
	 */
	
	//	if (argc == 1) usage();
	i=1;

        if (argc == 1) usage();

	while (i < argc) {
	        if (!strcmp(argv[i],"-s")) {
			++i;
			if (i >= argc) usage();
			nSmooth = atoi(argv[i]);
			++i;
			}
	        else if (!strcmp(argv[i],"-diag")) {
		        bOutDiag = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-std")) {
		        bStandard = 1;
			++i;
		        }
		else if (!strcmp(argv[i],"-nSplitting")) {
		  ++i;
		  if (i >= argc) usage();
		  nSplitting = atof(argv[i]);
		  ++i;
		}
		else if (!strcmp(argv[i],"-radius")) {
		  ++i;
		  if (i >= argc) usage();
		  radius = atof(argv[i]);
		  ++i;
		}
		else if (!strcmp(argv[i],"-Lc")) {
		  ++i;
		  if (i >= argc) usage();
		  Lc = atof(argv[i]);
		  ++i;
		}
		else if (!strcmp(argv[i],"-alpha")) {
		  ++i;
		  if (i >= argc) usage();
		  alpha = atof(argv[i]);
		  ++i;
		}
		else if (!strcmp(argv[i], "-ckpnt")) {
		  bCkpnt = 1;
		  ++i;
		}

		else if (!strcmp(argv[i], "-smbh")) {
		  bSmbh = 1;
		  ++i;
		}
		else if (!strcmp(argv[i], "-I")) {
		  strcpy(achName, argv[++i]);
		  i++;
		}
		else if(!strcmp(argv[i], "-aux")) {
		  bAux = 1;
		  ++i;
		}
		else usage();
	}

	fprintf(stderr, "nSmooth = %d\n", nSmooth);

	/* 
	** Initialize a KD context for each type of particle
	*/

	kdInit(&kdg,nBucket,fPeriod,fCenter,bOutDiag);
	kdInit(&kds,nBucket,fPeriod,fCenter,bOutDiag);
	//kdInit(&kds_old,nBucket,fPeriod,fCenter,bOutDiag);
	kdInit(&kdd,nBucket,fPeriod,fCenter,bOutDiag);

	if(bCkpnt) chk_header   = kdReadTipsyCheckpoint(kdg, kds, kdd, achName, bSmbh);
	else       tipsy_header = kdReadTipsy(kdg, kds, kdd, achName, 1, bSmbh, bAux);
	//else kdReadTipsy(kd,stdin);

	/*
	** Split the gas
	*/

	kdScatterActive(kdg,radius);
	kdBuildTree(kdg);
	smInit(&smx,kdg,nSmooth);
	kdTime(kdg,&sec1,&usec1);
	smDensityInit(smx,bPeriodic);
	kdTime(kdg,&sec1,&usec1);
	kdInitResample(kdg, nSplitting, 0, radius, Lc, alpha);
	kdBuildMoveTree(kdg);
	smTemperature(smx);
	smFinish(smx);
       

	/*
	** We don't split the dark, but only set the active particles
	*/

	if(kdd->nParticles > 0) {
	  kdScatterActive(kdd, radius);
          if (bSmbh) kdSmbhSoft(kdd, nSplitting);
        }
	/*
	** Split the stars
	*/

	// first, separate out the stars we don't want to split

	if(kds->nParticles > 0) {
	  
	  kdScatterActive(kds,radius);
	  kdBuildTree(kds);
	  smInit(&smx,kds,nSmooth);
	  kdTime(kds,&sec1,&usec1);
	  smDensityInit(smx,bPeriodic);
	  kdTime(kds,&sec1,&usec1);
	  kdInitResample(kds, nSplitting, kdg->nMove + kdd->nParticles, radius, Lc, alpha);
	  kdBuildMoveTree(kds);
	  smTemperature(smx);
	  smFinish(smx);

	}

	/* strcpy(achFile, achName); */
/* 	strcat(achFile, ".temps"); */
/* 	kdOutTemperature(kds, achFile); */
	
	
	strcpy(achFile, achName);
	strcat(achFile, ".subsample.tipsy");

	kdWriteTipsyStd(kdg, kdd, kds, achFile, radius, nSplitting);

	kdSetIord(kdg, kdd, kds, radius);

	strcpy(achFile, achName);
	strcat(achFile, ".subsample.chk");
	
	if(!bCkpnt) {
	  /*
	     if we started from a checkpoint, keep the same header,
	     otherwise set appropriate values
	  */
	  
	  chk_header.version = 8;
	  chk_header.not_corrupt_flag = 1;
	  chk_header.current_timestep = 0;
	  chk_header.current_time = tipsy_header.time;	
	}

	kdWriteTipsyCheckpoint(kdg, kdd, kds, chk_header, achFile, radius);

	fflush(stdout);
	kdFinish(kdg);
	kdFinish(kds);
	kdFinish(kdd);
	return 0;
}
	

