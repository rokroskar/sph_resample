#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <fenv.h>
#include "kd.h"
#include "smooth1.h"


#define PRUNE_STEPS		5
#define MICRO_STEP		0.1


void usage(void)
{
	fprintf(stderr,"USAGE:\n");
	fprintf(stderr,"sph_rsmpl -nSplitting <N> -radius <R> -s <nSmooth> INPUT_FILE \n");
	fprintf(stderr,"	 reads TIPSY CHECKPOINT input file from stdin\n");
	fprintf(stderr,"-diag : print diagnostic output\n");
	fprintf(stderr,"OUTPUT filenames are \n");
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
	int nSmooth,bCkpnt, bSmbh;
	int nMaxMembers;
	int bPeriodic;
	int bStandard;
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
	CHK_HEADER header;
	float radius;

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
	
	/*
	 ** Default output parameters.
	 */
	strcpy(achName,"split");
	
	/*
	 ** Now get the command line arguments!
	 */
	i = 1;
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
		else if (!strcmp(argv[i], "-ckpnt")) {
		  bCkpnt = 1;
		  ++i;
		}

		else if (!strcmp(argv[i], "-smbh")) {
		  bSmbh = 1;
		  ++i;
		}
		
		else usage();
		}

	/* 
	** Initialize a KD context for each type of particle
	*/
	
	kdInit(&kdg,nBucket,fPeriod,fCenter,bOutDiag);
	kdInit(&kds,nBucket,fPeriod,fCenter,bOutDiag);
	//kdInit(&kds_old,nBucket,fPeriod,fCenter,bOutDiag);
	kdInit(&kdd,nBucket,fPeriod,fCenter,bOutDiag);

	if(bCkpnt) header = kdReadTipsyCheckpoint(kdg, kds, kdd, stdin, bSmbh);
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
	kdInitResample(kdg, nSplitting, 0, 0, radius);
	kdBuildMoveTree(kdg);
	smTemperature(smx);
	smFinish(smx);
       

	/*
	** We don't split the dark, but only set the active particles
	*/
	
	kdScatterActive(kdd, radius);


	/*
	** Split the stars
	*/

	// first, separate out the stars we don't want to split

	kdScatterActive(kds,radius);
	kdBuildTree(kds);
	smInit(&smx,kds,nSmooth);
	kdTime(kds,&sec1,&usec1);
	smDensityInit(smx,bPeriodic);
	kdTime(kds,&sec1,&usec1);
	kdInitResample(kds, nSplitting, kdg->nMove + kdd->nParticles, 1, radius);
	kdBuildMoveTree(kds);
	smTemperature(smx);
	smFinish(smx);

	/* strcpy(achFile, achName); */
/* 	strcat(achFile, ".temps"); */
/* 	kdOutTemperature(kds, achFile); */
	
	
	strcpy(achFile, achName);
	strcat(achFile, ".subsample.tipsy");

	kdWriteTipsyStd(kdg, kdd, kds, achFile, radius, nSplitting);

	kdSetIord(kdg, kdd, kds, radius);

	strcpy(achFile, achName);
	strcat(achFile, ".subsample.chk");
	
	if (bCkpnt) kdWriteTipsyCheckpoint(kdg, kdd, kds, header, achFile, radius);

	fflush(stdout);
	kdFinish(kdg);
	kdFinish(kds);
	kdFinish(kdd);
	return 0;
}
	

