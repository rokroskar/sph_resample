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
	fprintf(stderr,"skid -tau <fLinkLength> [OPTIONAL ARGUMENTS]\n");
	fprintf(stderr,"	 reads TIPSY BINARY input file from stdin\n");
	fprintf(stderr,"     [-std]\n");
	fprintf(stderr,"COSMOLOGY and UNITS arguments:\n");
	fprintf(stderr,"     [-z <fRedShift>] [-O <fOmega>]\n");
	fprintf(stderr,"     [-G <fGravConst>] [-H <fHubble>]\n");
	fprintf(stderr,"     [-Lambda <fLambda>] [-Q <fQuintessense]\n");
	fprintf(stderr,"GROUP FINDING arguments (see man page!):\n");
	fprintf(stderr,"     [-s <nSmooth>] [-d <fMinDensity>] [-t <fMaxTemp>]\n");
	fprintf(stderr,"     [-cvg <fConvergeRadius>] [-scoop <fScoopRadius>]\n");
	fprintf(stderr,"     [-m <nMinMembers>] [-nu] [-gd] [-unbind <GroupName>[.grp]]\n");
	fprintf(stderr,"     [-M <fMaxMass>] [-fic] [-go] [-maxgroup nMaxMembers]\n");
	fprintf(stderr,"GRAVITATIONAL SOFTENING arguments:\n");
	fprintf(stderr,"     [-spline] [-plummer] [-e <fSoft>]\n"); 
	fprintf(stderr,"PERIODIC BOX specification:\n");
	fprintf(stderr,"     [-p <xyzPeriod>]\n");
	fprintf(stderr,"     [-c <xyzCenter>]\n");
	fprintf(stderr,"     [-cx <xCenter>] [-cy <yCenter>] [-cz <zCenter>]\n");
	fprintf(stderr,"OUTPUT arguments:\n");
	fprintf(stderr,"     [-o <Output Name>] [-ray] [-den] [-stats] [-diag]\n");
	fprintf(stderr,"\nSee man page skid(1).\n");
	exit(1);
	}

int main(int argc,char **argv)
{
        KD kdg, kds, kdd;
	SMX smx;
	/*
	 ** Input argument variables and control.
	 */
	int bTau,bCvg,bScoop,nSmooth,nMembers,bNoUnbind,bUnbindOnly,iSoftType,bCkpnt;
	int bEps,bOutRay,bOutDens,bGasAndDark,bGasOnly,bOutStats,bForceInitialCut;
	int nMaxMembers;
	int bPeriodic;
	int bStandard;
	float fTau,z,Omega0,G,H0,fDensMin,fTempMax,fMassMax,fCvg,fScoop,fEps;
	float Lambda;		/* Cosmological constant */
	float fQuintess;	/* Omega_quintessence */
	float fPeriod[3],fCenter[3];
	char achGroup[256],achName[256];
	/*
	 ** Working variables for SKID mainline.
	 */
	float fStep;
	int nBucket,i,j,nActive,nIttr;
	int sec1,usec1;
	int sec2,usec2;
	int sec3,usec3;
	int sec4,usec4;
	int sec5,usec5;
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
	 ** bTau flag to make sure user has at least specified the -tau argument.
	 */
	bTau = 0;
	/*
	 * Default tipsy native format.
	 */
	bStandard = 0;
	/*
	 ** Default Cosmological parameters.
	 */
	z = 0.0;
	Omega0 = 1.0;
	Lambda = 0.0;
	fQuintess = 0.0;
	G = 1.0;
	H0 = 0.0;
	/*
	 ** Default group finding parameters (those not dependent on fTau).
	 */
	nSmooth = 64;
	fDensMin = 0.0;
	fTempMax = HUGE;
	fMassMax = HUGE;
	bCvg = 0;
	bScoop = 0;
	nMembers = 8;
	nMaxMembers = INT_MAX;
	bNoUnbind = 0;
	bGasAndDark = 0;
	bGasOnly = 0;
	bUnbindOnly = 0;
	bForceInitialCut = 0;
	bPeriodic = 0;
	/*
	 ** Default gravitational parameters.
	 */
	bEps = 0;
	iSoftType = SPLINE;
	/*
	 ** Default periodic box parameters.
	 */
	for (j=0;j<3;++j) {
		fPeriod[j] = HUGE;
		fCenter[j] = 0.0;
		}
	/*
	 ** Default output parameters.
	 */
	strcpy(achName,"skid");
	bOutRay = 0;
	bOutDens = 0;
	bOutStats = 0;
	bOutDiag = 0;
	/*
	 ** Now get the command line arguments!
	 */
	i = 1;
	while (i < argc) {
		if (!strcmp(argv[i],"-tau")) {
			++i;
			if (i >= argc) usage();
			fTau = atof(argv[i]);
			bTau = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-z")) {
			++i;
			if (i >= argc) usage();
			z = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-O")) {
			++i;
			if (i >= argc) usage();
			Omega0 = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-Lambda")) {
			++i;
			if (i >= argc) usage();
			Lambda = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-Q")) {
			++i;
			if (i >= argc) usage();
			fQuintess = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-G")) {
			++i;
			if (i >= argc) usage();
			G = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-H")) {
			++i;
			if (i >= argc) usage();
			H0 = atof(argv[i]);
			++i;
			}
	    else if (!strcmp(argv[i],"-s")) {
			++i;
			if (i >= argc) usage();
			nSmooth = atoi(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-d")) {
			++i;
			if (i >= argc) usage();
			fDensMin = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-t")) {
			++i;
			if (i >= argc) usage();
			fTempMax = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-M")) {
			++i;
			if (i >= argc) usage();
			fMassMax = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-fic")) {
			bForceInitialCut = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-cvg")) {
			++i;
			if (i >= argc) usage();
			fCvg = atof(argv[i]);
			bCvg = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-scoop")) {
			++i;
			if (i >= argc) usage();
			fScoop = atof(argv[i]);
			bScoop = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-m")) {
			++i;
			if (i >= argc) usage();
			nMembers = atoi(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-maxgroup")) {
			++i;
			if (i >= argc) usage();
			nMaxMembers = atoi(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-nu")) {
			bNoUnbind = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-gd")) {
			bGasAndDark = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-go")) {
			bGasOnly = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-unbind")) {
			++i;
			if (i >= argc) usage();
			strcpy(achGroup,argv[i]);
			bUnbindOnly = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-spline")) {
			iSoftType = SPLINE;
			++i;
			}
		else if (!strcmp(argv[i],"-plummer")) {
			iSoftType = PLUMMER;
			++i;
			}
		else if (!strcmp(argv[i],"-e")) {
			++i;
			if (i >= argc) usage();
			fEps = atof(argv[i]);
			bEps = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-p")) {
			++i;
			if (i >= argc) usage();
			fPeriod[0] = atof(argv[i]);
			fPeriod[1] = atof(argv[i]);
			fPeriod[2] = atof(argv[i]);
			bPeriodic = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-c")) {
			++i;
			if (i >= argc) usage();
			fCenter[0] = atof(argv[i]);
			fCenter[1] = atof(argv[i]);
			fCenter[2] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-cx")) {
			++i;
			if (i >= argc) usage();
			fCenter[0] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-cy")) {
			++i;
			if (i >= argc) usage();
			fCenter[1] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-cz")) {
			++i;
			if (i >= argc) usage();
			fCenter[2] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-o")) {
			++i;
			if (i >= argc) usage();
			strcpy(achName,argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-ray")) {
			bOutRay = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-den")) {
			bOutDens = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-stats")) {
			bOutStats = 1;
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
		  ++i;
		  bCkpnt = 1;
		  ++i;
		}

		else usage();
		}
	/*
	 ** Make sure user has specified the -tau argument.
	 */
	//	if (!bTau) usage();
	/*
	 ** Default other parameters if needed.
	 */
	if (!bCvg) fCvg = 0.5*fTau;
	if (!bScoop) fScoop = 2.0*fTau;
	fStep = 0.5*fCvg;

	kdInit(&kdg,nBucket,fPeriod,fCenter,bOutDiag);
	kdInit(&kds,nBucket,fPeriod,fCenter,bOutDiag);
	//kdInit(&kds_old,nBucket,fPeriod,fCenter,bOutDiag);
	kdInit(&kdd,nBucket,fPeriod,fCenter,bOutDiag);

	if(bCkpnt) header = kdReadTipsyCheckpoint(kdg, kds, kdd, stdin);
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
	

