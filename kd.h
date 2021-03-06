#ifndef KD_HINCLUDED
#define KD_HINCLUDED

#include "cosmo.h"


#define ROOT		1
#define LOWER(i)	(i<<1)
#define UPPER(i)	((i<<1)+1)
#define PARENT(i)	(i>>1)
#define SIBLING(i) 	((i&1)?i-1:i+1)
#define SETNEXT(i)\
{\
	while (i&1) i=i>>1;\
	++i;\
	}

#define DARK	1
#define GAS		2
#define STAR	4

/*
 ** Softening types!
 */
#define PLUMMER		1
#define SPLINE		2



typedef struct CoolingParticleStruct {
  double Y_HI,Y_HeI,Y_HeII;        /* Abundance of ions */
} COOLPARTICLE;


typedef struct chk_particle {
  int iOrder;
  double fMass;
  double fSoft;
  double r[3];
  double v[3];
  double u;
  double fMetals;
  COOLPARTICLE CoolParticle;
  double fTimeCoolIsOffUntil;
  double fTimeForm;
  double fMassForm;
  double fNSN;
  double fMFracOxygen;
  double fMFracIron;
  int iGasOrder;
} CHK_PART;

typedef struct chkHeader {
  int version;  /* 8 */
  int not_corrupt_flag;
  int number_of_particles;
  int number_of_gas_particles;
  int number_of_dark_particles;
  int number_of_star_particles;
  int max_order;
  int max_order_gas;
  int max_order_dark;
  int current_timestep;
  double current_time;
  double current_ecosmo;
  double old_time;
  double old_potentiale;
  int opening_type;
  double opening_criterion;
  int number_of_outs;
  //  double redout_array;
  int a[16];
  double b[20];
  int c[9];
  double d[4];
  //  int e[3];
  double f[16];
} CHK_HEADER;


typedef struct pInitial {
  int iOrder;
  double fMass;
  double fSoft;
  float r[3];
  float v[3];
  double u;
  double fMetals;
  COOLPARTICLE CoolParticle;
  double fTimeCoolIsOffUntil;
  double fTimeForm;
  double fMassForm;
  double fNSN;
  double fMFracOxygen;
  double fMFracIron;
  double fDensity;
  float fBall2;
  float fTemp;
  int iGasOrder;
} PINIT;

typedef struct pMoved {
  int iOrder;
  double fMass;
  double fSoft;
  float r[3];
  float v[3];
  double u;
  double fMetals;
  COOLPARTICLE CoolParticle;
  double fTimeCoolIsOffUntil;
  double fTimeForm;
  double fMassForm;
  double fNSN;
  double fMFracOxygen;
  double fMFracIron;
  double fDensity;
  float fBall2;
  float fTemp;
  float fWeights;
  float rOld[3];
  float a[3];
  int iGasOrder;
  int iParentOrder;
} PMOVE;


/* typedef struct pInitial { */
/* 	float r[3]; */
/* 	float v[3]; */
/* 	float fMass; */
/* 	float fSoft; */
/* 	float fTemp; */
/* 	float fBall2; */
/* 	float fDensity; */
/*         float fMetals; */
/* 	int iOrder; */
/* } PINIT; */

/* /\* typedef struct pMoved { *\/ */
/* /\* 	float r[3]; *\/ */
/* /\* 	float rOld[3]; *\/ */
/* /\* 	float a[3]; *\/ */
/* /\* 	int iOrder; *\/ */
/* /\* 	} PMOVE; *\/ */

/* typedef struct pMoved { */
/* 	float r[3]; */
/* 	float v[3]; */
/* 	float fMass; */
/* 	float fSoft; */
/* 	float fTemp; */
/* 	float fBall2; */
/* 	float fDensity; */
/*         float fWeights; */
/*         float fMetals; */
/* 	int iOrder; */
/*         int iParentOrder; */
/*         float rOld[3]; */
/*         float a[3]; */
/* 	} PMOVE; */


typedef struct pGroup {
	float rel[3];
	float rCenter[3];
	float rBound[3];
	float vcm[3];
	float fMass;
	float fRadius;
	int nMembers;
	int pStart;
	int pCurr;
	} PGROUP;

typedef struct bndBound {
	float fMin[3];
	float fMax[3];
	} BND;

typedef struct kdNode {
	float fSplit;
	BND bnd;
	int iDim;
	int pLower;
	int pUpper;
	} KDN;


typedef struct kdContext {
	int nBucket;
	float fPeriod[3];
	float fCenter[3];
	float G;
	CSM csm;
	float z;
	int nParticles;
	int nDark;
	int nGas;
	int nStar;
	int inType;
	float fTime;
	int nLevels;
	int nNodes;
	int nSplit;
	int nMove;
	int nActive;
	int nInitActive;
	PINIT *pInit;
	PMOVE *pMove;
	PGROUP *pGroup;
	KDN *kdNodes;
	int nGroup;
	int *piGroup;
	int uSecond;
	int uMicro;
	int bOutDiag;
      	} * KD;


#define INTERSECTNP(pkdn,x,y,z,fDist2)\
{\
	float INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;\
	float INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1;\
	INTRSCT_dx = (pkdn)->bnd.fMin[0] - x;\
	INTRSCT_dx1 = x - (pkdn)->bnd.fMax[0];\
	INTRSCT_dy = (pkdn)->bnd.fMin[1] - y;\
	INTRSCT_dy1 = y - (pkdn)->bnd.fMax[1];\
	INTRSCT_dz = (pkdn)->bnd.fMin[2] - z;\
	INTRSCT_dz1 = z - (pkdn)->bnd.fMax[2];\
	if (INTRSCT_dx > 0.0) fDist2 = INTRSCT_dx*INTRSCT_dx;\
	else if (INTRSCT_dx1 > 0.0) fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
	else fDist2 = 0.0;\
	if (INTRSCT_dy > 0.0) fDist2 += INTRSCT_dy*INTRSCT_dy;\
	else if (INTRSCT_dy1 > 0.0) fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
	if (INTRSCT_dz > 0.0) fDist2 += INTRSCT_dz*INTRSCT_dz;\
	else if (INTRSCT_dz1 > 0.0) fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
	}


#define INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz)\
{\
	float INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;\
	float INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1,INTRSCT_fDist2;\
	INTRSCT_dx = c[cp].bnd.fMin[0]-x;\
	INTRSCT_dx1 = x-c[cp].bnd.fMax[0];\
	INTRSCT_dy = c[cp].bnd.fMin[1]-y;\
	INTRSCT_dy1 = y-c[cp].bnd.fMax[1];\
	INTRSCT_dz = c[cp].bnd.fMin[2]-z;\
	INTRSCT_dz1 = z-c[cp].bnd.fMax[2];\
	if (INTRSCT_dx > 0.0) {\
		INTRSCT_dx1 += lx;\
		if (INTRSCT_dx1 < INTRSCT_dx) {\
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
			sx = x+lx;\
			}\
		else {\
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
			sx = x;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (INTRSCT_dx1 > 0.0) {\
		INTRSCT_dx += lx;\
		if (INTRSCT_dx < INTRSCT_dx1) {\
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
			sx = x-lx;\
			}\
		else {\
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
			sx = x;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		INTRSCT_fDist2 = 0.0;\
		sx = x;\
		}\
	if (INTRSCT_dy > 0.0) {\
		INTRSCT_dy1 += ly;\
		if (INTRSCT_dy1 < INTRSCT_dy) {\
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
			sy = y+ly;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
			sy = y;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (INTRSCT_dy1 > 0.0) {\
		INTRSCT_dy += ly;\
		if (INTRSCT_dy < INTRSCT_dy1) {\
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
			sy = y-ly;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
			sy = y;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sy = y;\
		}\
	if (INTRSCT_dz > 0.0) {\
		INTRSCT_dz1 += lz;\
		if (INTRSCT_dz1 < INTRSCT_dz) {\
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
			sz = z+lz;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
			sz = z;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (INTRSCT_dz1 > 0.0) {\
		INTRSCT_dz += lz;\
		if (INTRSCT_dz < INTRSCT_dz1) {\
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
			sz = z-lz;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
			sz = z;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sz = z;\
		}\
	}


#define INTERCONT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz)\
{\
	float INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;\
	float INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1;\
	float INTRSCT_fDist2,INTRSCT_fMax2;\
	INTRSCT_dx = c[cp].bnd.fMin[0]-x;\
	INTRSCT_dx1 = x-c[cp].bnd.fMax[0];\
	INTRSCT_dy = c[cp].bnd.fMin[1]-y;\
	INTRSCT_dy1 = y-c[cp].bnd.fMax[1];\
	INTRSCT_dz = c[cp].bnd.fMin[2]-z;\
	INTRSCT_dz1 = z-c[cp].bnd.fMax[2];\
	if (INTRSCT_dx > 0.0) {\
		if (INTRSCT_dx1+lx < INTRSCT_dx) {\
			INTRSCT_dx1 += lx;\
			INTRSCT_dx -= lx;\
			sx = x+lx;\
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
			INTRSCT_fMax2 = INTRSCT_dx*INTRSCT_dx;\
			}\
		else {\
			sx = x;\
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
			INTRSCT_fMax2 = INTRSCT_dx1*INTRSCT_dx1;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (INTRSCT_dx1 > 0.0) {\
		if (INTRSCT_dx+lx < INTRSCT_dx1) {\
		    INTRSCT_dx += lx;\
			INTRSCT_dx1 -= lx;\
			sx = x-lx;\
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
			INTRSCT_fMax2 = INTRSCT_dx1*INTRSCT_dx1;\
			}\
		else {\
			sx = x;\
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
			INTRSCT_fMax2 = INTRSCT_dx*INTRSCT_dx;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sx = x;\
		INTRSCT_fDist2 = 0.0;\
		if (INTRSCT_dx < INTRSCT_dx1) INTRSCT_fMax2 = INTRSCT_dx*INTRSCT_dx;\
		else INTRSCT_fMax2 = INTRSCT_dx1*INTRSCT_dx1;\
		}\
	if (INTRSCT_dy > 0.0) {\
		if (INTRSCT_dy1+ly < INTRSCT_dy) {\
		    INTRSCT_dy1 += ly;\
			INTRSCT_dy -= ly;\
			sy = y+ly;\
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
			INTRSCT_fMax2 += INTRSCT_dy*INTRSCT_dy;\
			}\
		else {\
			sy = y;\
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
			INTRSCT_fMax2 += INTRSCT_dy1*INTRSCT_dy1;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (INTRSCT_dy1 > 0.0) {\
		if (INTRSCT_dy+ly < INTRSCT_dy1) {\
		    INTRSCT_dy += ly;\
			INTRSCT_dy1 -= ly;\
			sy = y-ly;\
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
			INTRSCT_fMax2 += INTRSCT_dy1*INTRSCT_dy1;\
			}\
		else {\
			sy = y;\
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
			INTRSCT_fMax2 += INTRSCT_dy*INTRSCT_dy;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sy = y;\
		if (INTRSCT_dy < INTRSCT_dy1) INTRSCT_fMax2 += INTRSCT_dy*INTRSCT_dy;\
		else INTRSCT_fMax2 += INTRSCT_dy1*INTRSCT_dy1;\
		}\
	if (INTRSCT_dz > 0.0) {\
		if (INTRSCT_dz1+lz < INTRSCT_dz) {\
		    INTRSCT_dz1 += lz;\
            INTRSCT_dz -= lz;\
			sz = z+lz;\
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
			INTRSCT_fMax2 += INTRSCT_dz*INTRSCT_dz;\
			}\
		else {\
			sz = z;\
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
			INTRSCT_fMax2 += INTRSCT_dz1*INTRSCT_dz1;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (INTRSCT_dz1 > 0.0) {\
		if (INTRSCT_dz+lz < INTRSCT_dz1) {\
			INTRSCT_dz += lz;\
		    INTRSCT_dz1 -= lz;\
			sz = z-lz;\
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
			INTRSCT_fMax2 += INTRSCT_dz1*INTRSCT_dz1;\
			}\
		else {\
			sz = z;\
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
			INTRSCT_fMax2 += INTRSCT_dz*INTRSCT_dz;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sz = z;\
		if (INTRSCT_dz < INTRSCT_dz1) INTRSCT_fMax2 += INTRSCT_dz*INTRSCT_dz;\
		else INTRSCT_fMax2 += INTRSCT_dz1*INTRSCT_dz1;\
		}\
	if (INTRSCT_fMax2 < fBall2) goto ContainedCell;\
	}


void kdTime(KD,int *,int *);
int kdInit(KD *,int,float *,float *,int);
void kdSetSoft(KD,float);
void kdSetUniverse(KD,float,float,float,float,float,float);
int kdParticleType(KD,int);
struct chkHeader kdReadTipsyCheckpoint(KD, KD, KD, char *, int);
void kdcofm(KD, KD, KD, int);
struct dump kdReadTipsy(KD, KD, KD, char *,int,int,int);
int kdBuildTree(KD);
int kdBuildMoveTree(KD);
int kdInitMove(KD,float,float,float,float,int,int);
int kdInitResample(KD, int, int, float, float, float);
int kdScatterActive(KD,float);
int bAllowInitialCut(KD);
void kdMoveParticles(KD,float);
int kdPruneInactive(KD,float);
void kdReactivateMove(KD);
void kdFoF(KD,float);
void kdInGroup(KD,char *);
void kdInitpGroup(KD);
void kdCalcCenter(KD);
void kdReadCenter(KD,char *,int);
void kdGroupOrder(KD);
void kdTooSmall(KD,int);
void kdUnbind(KD,int,float,int,int,int);
void kdOutGroup(KD,char *);
void kdOutDensity(KD,char *);
void kdOutVector(KD,char *);
void kdWriteGroup(KD,char *,int);
void kdOutStats(KD kd,char *,float,float);
void kdFinish(KD);
void kdOutTemperature(KD, char *);
void kdWriteTipsyStd(KD, KD, KD, char *, float, int);
void kdWriteTipsyCheckpoint(KD, KD, KD, struct chkHeader, char *, float);
void kdSetIord(KD, KD, KD, float);
void kdWriteMovedParticle(CHK_PART *, int, PMOVE);
void kdWriteInitParticle(CHK_PART *, int, PINIT);
float *kdReadFloatArray(char *, int, char *);
long *kdReadLongArray(char *, int, char *);

#endif






