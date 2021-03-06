#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include "kd.h"
#include "grav.h"
#include "tipsydefs.h"


#define MSOLG 1.99e33
/* G in cgs */
#define GCGS 6.67e-8
/* kiloparsec in centimeters */
#define KPCCM 3.085678e21
/* Gas constant */
#define Rgas 8.2494e7

int xdrHeader(XDR *pxdrs,struct dump *ph)
{
	int pad = 0;
	
	if (!xdr_double(pxdrs,&ph->time)) return 0;
	if (!xdr_int(pxdrs,&ph->nbodies)) return 0;
	if (!xdr_int(pxdrs,&ph->ndim)) return 0;
	if (!xdr_int(pxdrs,&ph->nsph)) return 0;
	if (!xdr_int(pxdrs,&ph->ndark)) return 0;
	if (!xdr_int(pxdrs,&ph->nstar)) return 0;
	if (!xdr_int(pxdrs,&pad)) return 0;
	return 1;
	}

void kdTime(KD kd,int *puSecond,int *puMicro)
{
	struct rusage ru;

	getrusage(0,&ru);
	*puMicro = ru.ru_utime.tv_usec - kd->uMicro;
	*puSecond = ru.ru_utime.tv_sec - kd->uSecond;
	if (*puMicro < 0) {
		*puMicro += 1000000;
		*puSecond -= 1;
		}
	kd->uSecond = ru.ru_utime.tv_sec;
	kd->uMicro = ru.ru_utime.tv_usec;
	}


int kdInit(KD *pkd,int nBucket,float *fPeriod,float *fCenter,int bOutDiag)
{
	KD kd;
	int j;

	kd = (KD)malloc(sizeof(struct kdContext));
	assert(kd != NULL);
	kd->nBucket = nBucket;
	for (j=0;j<3;++j) {
		kd->fPeriod[j] = fPeriod[j];
		kd->fCenter[j] = fCenter[j];
		}
	kd->bOutDiag = bOutDiag;
	kd->G = 1.0;
	kd->csm = NULL;
	kd->z = 0.0;
	kd->nParticles = 0;
	kd->nDark = 0;
	kd->nGas = 0;
	kd->nStar = 0;
	kd->inType = 0;
	kd->fTime = 0.0;
	kd->nLevels = 0;
	kd->nNodes = 0;
	kd->nSplit = 0;
	kd->nMove = 0;
	kd->nActive = 0;
	kd->nInitActive = 0;
	kd->pMove = NULL;
	kd->pInit = NULL;
	kd->pGroup = NULL;
	kd->kdNodes = NULL;
	kd->piGroup = NULL;
	kd->nGroup = 0;
	*pkd = kd;
	return(1);
	}

int kdParticleType(KD kd,int iOrder)
{
	if (iOrder < kd->nGas) return(GAS);
	else if (iOrder < (kd->nGas + kd->nDark)) return(DARK);
	else if (iOrder < kd->nParticles) return(STAR);
	else return(0);
	}


struct dump kdReadTipsy(KD kdg, KD kds, KD kdd, char *filename, int bStandard, int bSmbh, int bAux)
{
        PINIT *pg, *pd, *ps;
	int i,j;
	struct dump h;
	struct gas_particle gp;
	struct dark_particle dp;
	struct star_particle sp;
	float *fTimeForm, *fMassForm, *fTimeCoolIsOffUntil, *fNSN, *fMFracOxygen, *fMFracIron, *HI, *HeI, *HeII;
	long *iGasOrder;
	XDR xdrs;
	FILE *fp;
	
	if (kdg->bOutDiag) puts(">> kdReadTipsy()");
	fflush(stdout);
	
	if(fp = fopen(filename,"r"));
	else {
	  fprintf(stderr, "file \'%s\' not found - exiting promptly\n", filename);
	  exit(-1);
	}
	   

	if (bStandard) {
	    assert(sizeof(Real)==sizeof(float)); /* Otherwise, this XDR stuff
						    ain't gonna work */
	    xdrstdio_create(&xdrs, fp, XDR_DECODE);
	    xdrHeader(&xdrs,&h);
	    }
	else {
	    fread(&h,sizeof(struct dump),1,fp);
	    }

	/* 
	** Initialize and allocate arrays for the gas kd
	*/

	kdg->nDark = h.ndark;
	kdg->nGas = h.nsph;
	kdg->nStar = h.nstar;
	kdg->fTime = h.time;
	kdg->inType = GAS;

	kdg->nParticles = kdg->nGas;   // ONLY CONSIDERING GAS IN THIS CASE
	kdg->nInitActive = kdg->nParticles; 
	
	kdg->pInit = (PINIT *)malloc(kdg->nParticles*sizeof(PINIT));
	assert(kdg->pInit != NULL);
	
	/* 
	** Initialize and allocate arrays for the dark matter kd
	*/

	kdd->nDark = h.ndark;
	kdd->nGas = h.nsph;
	kdd->nStar = h.nstar;
	kdd->fTime = h.time;
	kdd->inType = DARK;
	
	kdd->nParticles = kdd->nDark;   // ONLY CONSIDERING DM IN THIS CASE
	kdd->nInitActive = kdd->nParticles; 
	
	kdd->pInit = (PINIT *)malloc(kdd->nParticles*sizeof(PINIT));
	assert(kdd->pInit != NULL);
	
	/* 
	** Initialize and allocate arrays for the stars kd
	*/
	
	kds->nDark = h.ndark;
	kds->nGas = h.nsph;
	kds->nStar = h.nstar;
	kds->fTime = h.time;
	kds->inType = STAR;
	
	kds->nParticles = kds->nStar;   // ONLY CONSIDERING STARS IN THIS CASE
	kds->nInitActive = kds->nParticles; 
	
	kds->pInit = (PINIT *)malloc(kds->nParticles*sizeof(PINIT));
	assert(kds->pInit != NULL);
	
	
	printf("nDark:%d nGas:%d nStar:%d\n",kdg->nDark,kdg->nGas,kdg->nStar);
	fflush(stdout);
	
	pg = kdg->pInit;
	ps = kds->pInit;
	pd = kdd->pInit;
  

	/*
	 ** Read Stuff!
	 */


	/* Read in the auxiliary files */

	if(bAux) {
	  fTimeCoolIsOffUntil = kdReadFloatArray(filename,0,"coolontime");
	  fTimeForm           = kdReadFloatArray(filename, 0, "timeform");
	  fMassForm           = kdReadFloatArray(filename, 0, "massform");
	  fNSN                = kdReadFloatArray(filename, 0, "ESNRate");
	  fMFracOxygen        = kdReadFloatArray(filename, 0, "OxMassFrac");
	  fMFracIron          = kdReadFloatArray(filename, 0, "FeMassFrac");
	  iGasOrder           = kdReadLongArray(filename, 0, "igasorder");
	  HI                  = kdReadFloatArray(filename, 0, "HI");
	  HeI                 = kdReadFloatArray(filename, 0, "HeI");
          HeII                = kdReadFloatArray(filename, 0, "HeII");
	}
	  
	for(i=0;i<kdg->nParticles;i++){
	  pg[i].iOrder = i;
	  pg[i].fDensity = 0.0;
	  if (bStandard) {
	    xdr_vector(&xdrs, (char *) &gp, 12,
		       sizeof(Real), (xdrproc_t) xdr_float);
	  }
	  else {
	    fread(&gp,sizeof(struct gas_particle),1,fp);
	  }
	  pg[i].fMass = gp.mass;
	  pg[i].fSoft = gp.hsmooth;
	  pg[i].fTemp = gp.temp;
	  pg[i].fMetals = gp.metals;
	  for (j=0;j<3;++j) {
	    pg[i].r[j] = gp.pos[j];
	    pg[i].v[j] = gp.vel[j];
	  }
	  
	  if(bAux) {
	    pg[i].fTimeCoolIsOffUntil = fTimeCoolIsOffUntil[i];
	    pg[i].fNSN                = fNSN[i];
	    pg[i].fMFracOxygen        = fMFracOxygen[i];
	    pg[i].fMFracIron          = fMFracIron[i];
	    pg[i].iGasOrder           = iGasOrder[i];
	    pg[i].CoolParticle.Y_HI   = HI[i];
	    pg[i].CoolParticle.Y_HeI  = HeI[i];
	    pg[i].CoolParticle.Y_HeII = HeII[i];
	  }
	}

	for(i=0;i<kdd->nParticles;i++){
	  if (bStandard) {
	    xdr_vector(&xdrs, (char *) &dp, 9,
		       sizeof(Real), (xdrproc_t) xdr_float);
	  }
	  else {
	    fread(&dp,sizeof(struct dark_particle),1,fp);
	  }
	  pd[i].fMass = dp.mass;
	  pd[i].fSoft = dp.eps;
	  pd[i].fTemp = 0.0;
	  for (j=0;j<3;++j) {
	    pd[i].r[j] = dp.pos[j];
	    pd[i].v[j] = dp.vel[j];
	  }
	}
	
	for(i=0;i<kds->nParticles;i++){
	  if (bStandard) {
	    xdr_vector(&xdrs, (char *) &sp, 11,
		       sizeof(Real), (xdrproc_t) xdr_float);
	  }
	  else {
	    fread(&sp,sizeof(struct star_particle),1,fp);
	  }
	  ps[i].fMass = sp.mass;
	  ps[i].fSoft = sp.eps;
	  ps[i].fTemp = 0.0;
	  ps[i].fMetals = sp.metals;
	  ps[i].fTimeForm = sp.tform;
	  for (j=0;j<3;++j) {
	    ps[i].r[j] = sp.pos[j];
	    ps[i].v[j] = sp.vel[j];
	  }
	  
	}
	
	if (bStandard) xdr_destroy(&xdrs);
	
	
	
	

	kdcofm(kdg, kdd, kds, bSmbh);
	
	if (kdg->bOutDiag) puts("<< kdReadTipsy()");
	fflush(stdout);
	fclose(fp);
	return h;
}

float *kdReadFloatArray(char *filename, int bAscii, char *arrayName)
// Reads a tipsy array file
{
  XDR xdrs;
  FILE *fp;
  long np;
  float *arr, temp;
  int i;
  char arrFile[256];

  strcpy(arrFile, filename);
  strcat(arrFile, ".");
  strcat(arrFile, arrayName);

  fprintf(stderr, "array = %s\n", arrFile);

  if (!bAscii) {
    assert(sizeof(Real)==sizeof(float)); /* Otherwise, this XDR stuff
					    ain't gonna work */
    
    fp = fopen(arrFile, "r");
    xdrstdio_create(&xdrs, fp, XDR_DECODE);
    xdr_long(&xdrs, &np);
    arr = malloc(sizeof(float)*np);
    for(i=0;i<np;i++) xdr_float(&xdrs,&temp);
  }
  fclose(fp);
  
  return arr; 
}

long *kdReadLongArray(char *filename, int bAscii, char *arrayName)
// Reads a tipsy array file
{
  XDR xdrs;
  FILE *fp;
  long np;
  long *arr, temp;
  int i;
  char arrFile[256];

  strcpy(arrFile, filename);
  strcat(arrFile, ".");
  strcat(arrFile, arrayName);

  fprintf(stderr, "array = %s\n", arrFile);

  if (!bAscii) {
    assert(sizeof(Real)==sizeof(float)); /* Otherwise, this XDR stuff
					    ain't gonna work */
    
    fp = fopen(arrFile, "r");
    xdrstdio_create(&xdrs, fp, XDR_DECODE);
    xdr_long(&xdrs, &np);
    arr = malloc(sizeof(float)*np);
    for(i=0;i<np;i++) xdr_long(&xdrs,&temp);
  }
  fclose(fp);
  
  return arr; 
}
    


struct chkHeader kdReadTipsyCheckpoint(KD kdg, KD kds, KD kdd, char *filename, int bSmbh)
{

  PINIT *pg, *ps, *pd;
  CHK_PART cp;
  CHK_HEADER h;
  int i,j;
  long offset;
  char name[100],fdlname[500],c[10000];
  FILE *IN, *FDL, *OUT, *fp;
 
  
  if (kdg->bOutDiag) puts(">> kdReadTipsyCheckpoint()");
  fflush(stdout);

  if(getenv("PKDGRAV_CHECKPOINT_FDL") == NULL) {
    fprintf(stderr,"PKDGRAV_CHECKPOINT_FDL environment variable not set\n");
    return;
  }
  else  FDL = fopen(getenv("PKDGRAV_CHECKPOINT_FDL"),"r");
  if(FDL == NULL) {
    fprintf(stderr,"error opening %s\n",getenv("PKDGRAV_CHECKPOINT_FDL"));
    return;
  }

  fp = fopen(filename, "r");

  while(!feof(FDL)) {
    fread(&c[0],sizeof(char),1,FDL);
    fscanf(fp,"%c",&c[0]);
  }

  fseek(fp,-2,SEEK_CUR);
 
  fread(&h,sizeof(CHK_HEADER),1,fp);

  fprintf(stderr, "Output time = %f\n", h.current_time);
  fprintf(stderr, "max order = %d\nmax gas order = %d\nmax dark order = %d\n", 
	  h.max_order, h.max_order_gas, h.max_order_dark);


  /* 
  ** Initialize and allocate arrays for the gas kd
  */

  kdg->nDark = h.number_of_dark_particles;
  kdg->nGas = h.number_of_gas_particles;
  kdg->nStar = h.number_of_star_particles;
  kdg->fTime = h.current_time;
  kdg->inType = GAS;

  kdg->nParticles = kdg->nGas;   // ONLY CONSIDERING GAS IN THIS CASE
  kdg->nInitActive = kdg->nParticles; 
  
  kdg->pInit = (PINIT *)malloc(kdg->nParticles*sizeof(PINIT));
  assert(kdg->pInit != NULL);

  /* 
  ** Initialize and allocate arrays for the dark matter kd
  */

  kdd->nDark = h.number_of_dark_particles;
  kdd->nGas = h.number_of_gas_particles;
  kdd->nStar = h.number_of_star_particles;
  kdd->fTime = h.current_time;
  kdd->inType = DARK;

  kdd->nParticles = kdd->nDark;   // ONLY CONSIDERING DM IN THIS CASE
  kdd->nInitActive = kdd->nParticles; 
  
  kdd->pInit = (PINIT *)malloc(kdd->nParticles*sizeof(PINIT));
  assert(kdd->pInit != NULL);

  /* 
  ** Initialize and allocate arrays for the stars kd
  */

  kds->nDark = h.number_of_dark_particles;
  kds->nGas = h.number_of_gas_particles;
  kds->nStar = h.number_of_star_particles;
  kds->fTime = h.current_time;
  kds->inType = STAR;

  kds->nParticles = kds->nStar;   // ONLY CONSIDERING STARS IN THIS CASE
  kds->nInitActive = kds->nParticles; 
  
  kds->pInit = (PINIT *)malloc(kds->nParticles*sizeof(PINIT));
  assert(kds->pInit != NULL);

  
  printf("nDark:%d nGas:%d nStar:%d\n",kdg->nDark,kdg->nGas,kdg->nStar);
  fflush(stdout);
  
  pg = kdg->pInit;
  ps = kds->pInit;
  pd = kdd->pInit;
  

  /*
  ** First set up the gas 
  */

  for (i=0;i<kdg->nParticles;i++) {

    fread(&cp, sizeof(CHK_PART), 1, fp);
    
    pg[i].iOrder = i;
    pg[i].fMass = cp.fMass;
    pg[i].fSoft = cp.fSoft;
    pg[i].u = cp.u;
    pg[i].fMetals = cp.fMetals;
    pg[i].CoolParticle = cp.CoolParticle;
    pg[i].fTimeCoolIsOffUntil = cp.fTimeCoolIsOffUntil;
    pg[i].fTimeForm = cp.fTimeForm;
    pg[i].fMassForm = cp.fMassForm;
    pg[i].fNSN = cp.fNSN;
    pg[i].fMFracOxygen = cp.fMFracOxygen;
    pg[i].fMFracIron = cp.fMFracIron;
    pg[i].iGasOrder = cp.iGasOrder;

    pg[i].fTimeForm = 0.0;
    pg[i].fMassForm = 0.0;
    //    pg[i].fNSN = 0.0;

    pg[i].fTemp = cp.u*(GCGS*2.3e5*MSOLG)/(1.0*KPCCM)/	     \
		       (2*(1.-0.25) - cp.CoolParticle.Y_HI + \
			3*0.25/4.0 - 2*cp.CoolParticle.Y_HeI - cp.CoolParticle.Y_HeII) \
			/ Rgas / 1.5;

    //pg[i].fTemp = cp.u*1.5*GCGS*2.3e5*MSOLG/(1.0*KPCCM)*1.38e-16;

    for (j = 0; j <= 2; j++) {
      pg[i].r[j] = cp.r[j];
      pg[i].v[j] = cp.v[j];
    }
      
  }

  /*
  ** Set up the dark
  */

  for (i=0;i<kdd->nParticles;i++) {

    fread(&cp, sizeof(CHK_PART), 1, fp);

    pd[i].iOrder = i+kdg->nGas;
    pd[i].fMass = cp.fMass;
    pd[i].fSoft = cp.fSoft;
    pd[i].u = cp.u;
    pd[i].fMetals = cp.fMetals;
    pd[i].CoolParticle = cp.CoolParticle;
    pd[i].fTimeCoolIsOffUntil = cp.fTimeCoolIsOffUntil;
    pd[i].fTimeForm = cp.fTimeForm;
    pd[i].fMassForm = cp.fMassForm;
    pd[i].fNSN = cp.fNSN;
    pd[i].fMFracOxygen = cp.fMFracOxygen;
    pd[i].fMFracIron = cp.fMFracIron;
    pd[i].iGasOrder = cp.iGasOrder;

    pd[i].fTimeForm = 0.0;
    pd[i].fMassForm = 0.0;
    pd[i].fNSN = 0.0;

    /* p[i].fTemp = cp.u*(GCGS*2.3e5*MSOLG)/(1.0*KPCCM)/	     \ */
/* 		       (2*(1.-0.25) - cp.CoolParticle.Y_HI + \ */
/* 			3*0.25/4.0 - 2*cp.CoolParticle.Y_HeI - cp.CoolParticle.Y_HeII \ */
/* 			* Rgas * 1.5); */

    //pd[i].fTemp = cp.u*1.5*GCGS*2.3e5*MSOLG/(1.0*KPCCM)*1.38e-16;

    for (j = 0; j <= 2; j++) {
      pd[i].r[j] = cp.r[j];
      pd[i].v[j] = cp.v[j];
    }
  }

  /* 
  ** Set up the stars
  */

  for (i=0;i<kds->nParticles;i++) {
    
    fread(&cp, sizeof(CHK_PART), 1, fp);

    ps[i].iOrder = i+kdg->nGas+kdg->nDark;
    ps[i].fMass = cp.fMass;
    ps[i].fSoft = cp.fSoft;
    ps[i].u = cp.u;
    ps[i].fMetals = cp.fMetals;
    ps[i].CoolParticle = cp.CoolParticle;
    ps[i].fTimeCoolIsOffUntil = cp.fTimeCoolIsOffUntil;
    ps[i].fTimeForm = cp.fTimeForm;
    ps[i].fMassForm = cp.fMassForm;
    if(cp.fTimeForm > 0)
      assert(cp.fMassForm > 0);
    ps[i].fNSN = cp.fNSN;
    ps[i].fMFracOxygen = cp.fMFracOxygen;
    ps[i].fMFracIron = cp.fMFracIron;
    ps[i].iGasOrder = cp.iGasOrder;

    
    for (j = 0; j <= 2; j++) {
      ps[i].r[j] = cp.r[j];
      ps[i].v[j] = cp.v[j];
    }

  }
  
  kdcofm(kdg, kdd, kds, bSmbh);

  if(kdg->bOutDiag) puts("<< kdReadTipsyCheckpoint()");
  fflush(stdout);

  return(h);
}

void kdcofm(KD kdg, KD kdd, KD kds, int bSmbh)
{

   /*
   ** Center the particle positions and velocities based on the 
   ** center of mass -> take the center between the two BHs and center all particles
   */
  PINIT *pg, *pd, *ps;
  float comr[3], comv[3];
  double totmass = 0.0;
  int i, j;

  pg = kdg->pInit;
  pd = kdd->pInit;
  ps = kds->pInit;

  if(bSmbh) 
    {
      
      /*for(i=0;i<3;i++) 
	{
	  comr[i] = (pd[0].r[i]+pd[1000001].r[i])/2.0;
	  comv[i] = (pd[0].v[i]+pd[1000001].v[i])/2.0;
	}
      */
      // set the cofm by hand
      comr[0] = 1.38185727;
      comr[1] = -2.50419233;
      comr[2] = -0.15508558;
      comv[0] = -0.64014463;
      comv[1] = -0.67666935;
      comv[2] = 0.77240984;
      

      // set the locations and velocities of SMBHs by hand to match the cofm pos/vel of progenitors

      pd[0].r[0] = -0.99154286; 
      pd[0].r[1] = -4.98014241;
      pd[0].r[2] = -0.19191479;
      pd[1000001].r[0] = 3.65599845;
      pd[1000001].r[1] = -0.13178999;
      pd[1000001].r[2] = -0.11979661;

      pd[0].v[0] = 203.81394998;
      pd[0].v[1] = 24.8927977;
      pd[0].v[2] = 0.81545728;
      pd[1000001].v[0] = -196.54367953;
      pd[1000001].v[1] = -25.17678508;
      pd[1000001].v[2] = 0.73116272;
    }
  else 
    {
      for(i=0;i<3;i++){
	comr[i] = 0.0;
	comv[i] = 0.0;
      }
      for(i=0;i<3;i++)
	{
	  for(j=0;j<kdg->nParticles;j++) 
	    {
	      comr[i] += pg[j].r[i]*pg[j].fMass;
	      comv[i] += pg[j].v[i]*pg[j].fMass;
	      if(i==0) totmass += pg[j].fMass;
	    }
	  for(j=0;j<kdd->nParticles;j++) 
	    {
	      comr[i] += pd[j].r[i]*pd[j].fMass;
	      comv[i] += pd[j].v[i]*pd[j].fMass;
	      if(i==0) totmass += pd[j].fMass;
	    }
	  for(j=0;j<kds->nParticles;j++) 
	    {
	      comr[i] += ps[j].r[i]*ps[j].fMass;
	      comv[i] += ps[j].v[i]*ps[j].fMass;
	      if(i==0)totmass += ps[j].fMass;
	    }
	  comr[i] /= totmass;
	  comv[i] /= totmass;
	  fprintf(stderr, "totmass = %f\n", totmass);
	}
    }
	  
	  
  for (i=0;i<kdg->nParticles;i++)
    {
      for (j = 0; j <= 2; j++) {
	pg[i].r[j] -= comr[j];
	pg[i].v[j] -= comv[j];
      }
    }

  for (i=0;i<kdd->nParticles;i++)
    {
      for (j = 0; j <= 2; j++) {
	pd[i].r[j] -= comr[j];
	pd[i].v[j] -= comv[j];
      }
    }

  for (i=0;i<kds->nParticles;i++)
    {
      for (j = 0; j <= 2; j++) {
	ps[i].r[j] -= comr[j];
	ps[i].v[j] -= comv[j];
      }
    }

    
}


void kdSelectInit(KD kd,int d,int k,int l,int r)
{
	PINIT *p,t;
	double v;
	int i,j;

	p = kd->pInit;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


void kdSelectMove(KD kd,int d,int k,int l,int r)
{
	PMOVE *p,t;
	double v;
	int i,j;

	p = kd->pMove;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


void Combine(KDN *p1,KDN *p2,KDN *pOut)
{
	int j;

	/*
	 ** Combine the bounds.
	 */
	for (j=0;j<3;++j) {
		if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
			pOut->bnd.fMin[j] = p2->bnd.fMin[j];
		else
			pOut->bnd.fMin[j] = p1->bnd.fMin[j];
		if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
			pOut->bnd.fMax[j] = p2->bnd.fMax[j];
		else
			pOut->bnd.fMax[j] = p1->bnd.fMax[j];
		}
	}


void UpPassInit(KD kd,int iCell)
{
	KDN *c;
	int l,u,pj,j;

	c = kd->kdNodes;
	if (c[iCell].iDim != -1) {
		l = LOWER(iCell);
		u = UPPER(iCell);
		UpPassInit(kd,l);
		UpPassInit(kd,u);
		Combine(&c[l],&c[u],&c[iCell]);
		}
	else {
		l = c[iCell].pLower;
		u = c[iCell].pUpper;
		for (j=0;j<3;++j) {
			c[iCell].bnd.fMin[j] = kd->pInit[u].r[j];
			c[iCell].bnd.fMax[j] = kd->pInit[u].r[j];
			}
		for (pj=l;pj<u;++pj) {
			for (j=0;j<3;++j) {
				if (kd->pInit[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = kd->pInit[pj].r[j];
				if (kd->pInit[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = kd->pInit[pj].r[j];
				}
			}
		}
	}


void UpPassMove(KD kd,int iCell)
{
	KDN *c;
	int l,u,pj,j;

	c = kd->kdNodes;
	if (c[iCell].iDim != -1) {
		l = LOWER(iCell);
		u = UPPER(iCell);
		UpPassMove(kd,l);
		UpPassMove(kd,u);
		Combine(&c[l],&c[u],&c[iCell]);
		}
	else {
		l = c[iCell].pLower;
		u = c[iCell].pUpper;
		for (j=0;j<3;++j) {
			c[iCell].bnd.fMin[j] = kd->pMove[u].r[j];
			c[iCell].bnd.fMax[j] = kd->pMove[u].r[j];
			}
		for (pj=l;pj<u;++pj) {
			for (j=0;j<3;++j) {
				if (kd->pMove[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = kd->pMove[pj].r[j];
				if (kd->pMove[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = kd->pMove[pj].r[j];
				}
			}
		}
	}


int kdBuildTree(KD kd)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	BND bnd;

	if (kd->bOutDiag) puts(">> kdBuildTree()");
	fflush(stdout);
	if (kd->nInitActive == 0) {
		if (kd->kdNodes) free(kd->kdNodes);
		kd->kdNodes = NULL;
		return(1);
		}
	assert(kd->nInitActive > 0);
	n = kd->nInitActive;
	kd->nLevels = 1;
	l = 1;
	while (n > kd->nBucket) {
		n = n>>1;
		l = l<<1;
		++kd->nLevels;
		}
	kd->nSplit = l;
	kd->nNodes = l<<1;
	if (kd->kdNodes) {
		free(kd->kdNodes);
		kd->kdNodes = NULL;
		}
	kd->kdNodes = (KDN *)malloc(kd->nNodes*sizeof(KDN));
	assert(kd->kdNodes != NULL);
	/*
	 ** Calculate Bounds.
	 */
	for (j=0;j<3;++j) {
		bnd.fMin[j] = kd->pInit[0].r[j];
		bnd.fMax[j] = kd->pInit[0].r[j];
		}
	for (i=1;i<kd->nInitActive;++i) {
		for (j=0;j<3;++j) {
			if (bnd.fMin[j] > kd->pInit[i].r[j]) 
				bnd.fMin[j] = kd->pInit[i].r[j];
			else if (bnd.fMax[j] < kd->pInit[i].r[j])
				bnd.fMax[j] = kd->pInit[i].r[j];
			}
		}
	/*
	 ** Set up ROOT node
	 */
	c = kd->kdNodes;
	c[ROOT].pLower = 0;
	c[ROOT].pUpper = kd->nInitActive-1;
	c[ROOT].bnd = bnd;
	i = ROOT;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;

			m = (c[i].pLower + c[i].pUpper)/2;
			kdSelectInit(kd,d,m,c[i].pLower,c[i].pUpper);

			c[i].fSplit = kd->pInit[m].r[d];
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m-1;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].pLower = m;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			assert(diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;
			SETNEXT(i);
			if (i == ROOT) break;
			}
		}
	UpPassInit(kd,ROOT);
	if (kd->bOutDiag) puts("<< kdBuildTree()");
	fflush(stdout);
	return(1);
	}


int kdBuildMoveTree(KD kd)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	BND bnd;

	if (kd->bOutDiag) puts(">> kdBuildMoveTree()");
	fflush(stdout);
	if (kd->nActive == 0) {
		if (kd->kdNodes) free(kd->kdNodes);
		kd->kdNodes = NULL;
		return(1);
		}
	assert(kd->nActive > 0);
	n = kd->nActive;
	kd->nLevels = 1;
	l = 1;
	while (n > kd->nBucket) {
		n = n>>1;
		l = l<<1;
		++kd->nLevels;
		}
	kd->nSplit = l;
	kd->nNodes = l<<1;
	if (kd->kdNodes) {
		free(kd->kdNodes);
		kd->kdNodes = NULL;
		}
	kd->kdNodes = (KDN *)malloc(kd->nNodes*sizeof(KDN));
	assert(kd->kdNodes != NULL);
	/*
	 ** Calculate Bounds.
	 */
	for (j=0;j<3;++j) {
		bnd.fMin[j] = kd->pMove[0].r[j];
		bnd.fMax[j] = kd->pMove[0].r[j];
		}
	for (i=1;i<kd->nMove;++i) {
		for (j=0;j<3;++j) {
			if (bnd.fMin[j] > kd->pMove[i].r[j]) 
				bnd.fMin[j] = kd->pMove[i].r[j];
			else if (bnd.fMax[j] < kd->pMove[i].r[j])
				bnd.fMax[j] = kd->pMove[i].r[j];
			}
		}
	/*
	 ** Set up ROOT node
	 */
	c = kd->kdNodes;
	c[ROOT].pLower = 0;
	c[ROOT].pUpper = kd->nActive-1;
	c[ROOT].bnd = bnd;
	i = ROOT;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;

			m = (c[i].pLower + c[i].pUpper)/2;
			kdSelectMove(kd,d,m,c[i].pLower,c[i].pUpper);

			c[i].fSplit = kd->pMove[m].r[d];
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m-1;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].pLower = m;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			assert(diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;
			SETNEXT(i);
			if (i == ROOT) break;
			}
		}
	UpPassMove(kd,ROOT);
	if (kd->bOutDiag) puts("<< kdBuildMoveTree()");
	fflush(stdout);
	return(1);
	}



int ScatterCriterion(KD kd, int pi, float radius)
{
  int i;
  float r[3], dr2;
  
  for (i=0;i<3;i++) r[i] = kd->pInit[pi].r[i];

  dr2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
  
  switch (kd->inType)
    {
    case (DARK):
      if(dr2 <= radius*radius) return(1);
      else return(0);
    case (GAS): 
      return(1);
    case (STAR):
      if(kd->pInit[pi].fTimeForm >= 0.0) return(1);
      else return(0);
    }
}



int kdInitResample(KD kd, int nSplitting, int offset, float radius,
		   float Lc, float alpha)
{
  /* 
   * Initialize the child particles based on parent particle properties
   */ 

        int pi,nCnt,j, nPart;
	float fCvg2;
	float xcurr, ycurr, zcurr, bs, currSplit, jx, jy, jz, jtot;
	PINIT *p;
	PMOVE *pm;
	
	if (kd->bOutDiag) puts(">> kdInitResample()");
	fflush(stdout);

	/* 
	** assign pointers
	*/

	p = kd->pInit;
	

	/*
	 ** First count number of particles meeting the criterion.
	 */

	kd->nMove = 0;
	nPart = 0;

	
	if (Lc == 0)
	  {
	    for (pi=0; pi<kd->nParticles;++pi) 
	      {
		xcurr = p[pi].r[0];
		ycurr = p[pi].r[1];
		zcurr = p[pi].r[2];
		
		
		
		if (((xcurr*xcurr + ycurr*ycurr + zcurr*zcurr) < radius*radius) 
		    && p[pi].fTimeForm >= 0.0)
		  nPart++;
	      }
	    fprintf(stderr, "nPart = %d\n", nPart);
	    kd->nMove = nPart*nSplitting;
	  }
	else
	  {
	    for (pi=0; pi<kd->nParticles;++pi)
	      {
		if(p[pi].fTimeForm >= 0.0) 
		  {
		    jx = p[pi].r[1]*p[pi].v[2] - p[pi].r[2]*p[pi].v[1];
		    jy = p[pi].r[2]*p[pi].v[0] - p[pi].r[0]*p[pi].v[2];
		    jz = p[pi].r[0]*p[pi].v[1] - p[pi].r[1]*p[pi].v[0];

		    jtot = sqrt(jx*jx+jy*jy+jz*jz);
		    nPart += (int)ceill((float)nSplitting/(1+jtot/Lc));
		  }
	      }
	    fprintf(stderr, "nPart = %d\n", nPart);
	    kd->nMove = nPart;
	  }
			
	kd->nActive = kd->nMove;
	
	/*
	 ** Allocate moving particles and initialize.
	 */

	kd->pMove = (PMOVE *)malloc(kd->nMove*sizeof(PMOVE));
	assert(kd->pMove != NULL);

	pm = kd->pMove;

	nCnt = 0;
	
	/*
	** initialize particle positions - cut out a sphere of radius ''radius''
	*/

	fprintf(stderr, "offset = %d\n", offset);

	for (pi=0;pi<kd->nInitActive;++pi) 
	  {
	    xcurr = p[pi].r[0];
	    ycurr = p[pi].r[1];
	    zcurr = p[pi].r[2];

	    bs = sqrt(p[pi].fBall2)/2.0;
	    currSplit = 0;

	    /* check if we are splitting by angular momentum or radius */
	    if (Lc > 0) 
	      {
		jx = p[pi].r[1]*p[pi].v[2] - p[pi].r[2]*p[pi].v[1];
		jy = p[pi].r[2]*p[pi].v[0] - p[pi].r[0]*p[pi].v[2];
		jz = p[pi].r[0]*p[pi].v[1] - p[pi].r[1]*p[pi].v[0];
		
		jtot = sqrt(jx*jx+jy*jy+jz*jz);
		currSplit = ceilf((int)nSplitting/(1.0+jtot/Lc));
		if (currSplit < 1) fprintf(stderr,"umm... strange\n");
		//fprintf(stderr, "j = %f\ncurrSplit = %f\n", j, currSplit);
	      }
	    else currSplit = nSplitting;
	    
	    assert(currSplit > 0);

	    if ((((xcurr*xcurr + ycurr*ycurr + zcurr*zcurr) < radius*radius) || Lc > 0) 
		&& p[pi].fTimeForm >= 0.0)
	      {
		for (j=0; j < (int)currSplit; j++) 
		  {
		    pm[nCnt].r[0] = xcurr + (float)rand()/(float)RAND_MAX*bs*2.0 - bs;
		    pm[nCnt].r[1] = ycurr + (float)rand()/(float)RAND_MAX*bs*2.0 - bs;
		    pm[nCnt].r[2] = zcurr + (float)rand()/(float)RAND_MAX*bs*2.0 - bs;
		   

		    pm[nCnt].v[0] = p[pi].v[0];
		    pm[nCnt].v[1] = p[pi].v[1];
		    pm[nCnt].v[2] = p[pi].v[2];
		    
		    pm[nCnt].fMass = p[pi].fMass/(float)currSplit;
		    pm[nCnt].fSoft = p[pi].fSoft/cbrt(currSplit);
		    assert(pm[nCnt].fSoft != 0.0);
		    pm[nCnt].iOrder = nCnt + offset;
		    		    
		    pm[nCnt].iParentOrder = p[pi].iOrder;
		    

		    /* these properties will ideally be given to each split particle based 
		       on the smoothing kernel, but sometimes that fails so we'll start 
		       by assigning them the properties if their parent particle - that can't
		       be too far off anyway
		    */
		    
		    pm[nCnt].fTemp = p[pi].fTemp;
		    pm[nCnt].fDensity = p[pi].fDensity;
		    pm[nCnt].fMetals = p[pi].fMetals/currSplit;
		    pm[nCnt].fMFracOxygen = p[pi].fMFracOxygen;
		    pm[nCnt].fMFracIron = p[pi].fMFracIron;
		    
		    pm[nCnt].CoolParticle.Y_HI = p[pi].CoolParticle.Y_HI;
		    pm[nCnt].CoolParticle.Y_HeI = p[pi].CoolParticle.Y_HeI;
		    pm[nCnt].CoolParticle.Y_HeII = p[pi].CoolParticle.Y_HeII;

		    pm[nCnt].fTimeCoolIsOffUntil = p[pi].fTimeCoolIsOffUntil;
		    pm[nCnt].fTimeForm = p[pi].fTimeForm;
		    pm[nCnt].fMassForm = p[pi].fMassForm/currSplit;
		    
		    nCnt++;
		  }
	
		}
	  }
	if (kd->bOutDiag) puts("<< kdInitResample()");
	fflush(stdout);
	fprintf(stderr, "total number of children = %d\n", nCnt);
	return(kd->nActive);
}




int kdScatterActive(KD kd, float radius)
{
	PINIT *p,t;
	int i,j;
	
	if (kd->bOutDiag) puts(">> kdScatterActive()");
	fflush(stdout);
	p = kd->pInit;
	i = 0;
	j = kd->nParticles-1;
	while (1) {
		while (ScatterCriterion(kd,i,radius))
			if (++i > j) goto done;
		while (!ScatterCriterion(kd,j,radius))
			if (i > --j) goto done;
		t = p[i];
		p[i] = p[j];
		p[j] = t;
		}
 done:
	kd->nInitActive = i;
	fprintf(stderr, "nInitActive = %d\n", kd->nInitActive);
	if (kd->bOutDiag) puts("<< kdScatterActive()");
	fflush(stdout);
	return(i);
	}



int CmpInit(const void *p1,const void *p2)
{
	PINIT *a = (PINIT *)p1;
	PINIT *b = (PINIT *)p2;

	return(a->iOrder - b->iOrder);
	}


int CmpMove(const void *p1,const void *p2)			
{
	PMOVE *a = (PMOVE *)p1;
	PMOVE *b = (PMOVE *)p2;

	return(a->iOrder - b->iOrder);
	}


void Order(KD kd)
{
	if (kd->bOutDiag) puts(">> kdOrder()");
	fflush(stdout);
	if (kd->pInit != NULL) {
		qsort(kd->pInit,kd->nParticles,sizeof(PINIT),CmpInit);
		}
	if (kd->pMove != NULL) {
		qsort(kd->pMove,kd->nMove,sizeof(PMOVE),CmpMove);
		}
	if (kd->bOutDiag) puts("<< kdOrder()");
	fflush(stdout);
	}




void kdOutDensity(KD kd,char *pszFile)
{	
	FILE *fp;
	int i;

	if (pszFile) {
		fp = fopen(pszFile,"w");
		assert(fp != NULL);
		}
	else {
		fp = stdout;
		}
	/*
	 ** Make sure the particles are ordered before outputing the densities.
	 */
	Order(kd);
	fprintf(fp,"%d\n",kd->nParticles);
	for (i=0;i<kd->nParticles;++i) {
		fprintf(fp,"%.10g\n",kd->pInit[i].fDensity);
		}
	fclose(fp);
	}


void kdOutVector(KD kd,char *pszFile)
{
	PINIT *p;
	PMOVE *pm;
	FILE *fp;
	int pi,pmi;
	float hx,hy,hz,dx,dy,dz;

	/*
	 ** Make sure the particles are ordered before outputing the vectors.
	 */
	Order(kd);

	p = kd->pInit;
	pm = kd->pMove;
	hx = 0.5*kd->fPeriod[0];
	hy = 0.5*kd->fPeriod[1];
	hz = 0.5*kd->fPeriod[2];
	fp = fopen(pszFile,"w");
	assert(fp != NULL);
	fprintf(fp,"%d\n",kd->nParticles);
	for (pmi=0,pi=0;pmi<kd->nMove;++pmi,++pi) {
		while (p[pi].iOrder < pm[pmi].iOrder) {
			fprintf(fp,"0\n");
			++pi;
			}
		assert(p[pi].iOrder == pm[pmi].iOrder);
		dx = pm[pmi].r[0] - p[pi].r[0];
		if (dx > hx) dx -= 2*hx;
		if (dx <= -hx) dx += 2*hx;
		fprintf(fp,"%g\n",dx);
		}
	for (;pi<kd->nParticles;++pi) fprintf(fp,"0\n");
	for (pmi=0,pi=0;pmi<kd->nMove;++pmi,++pi) {
		while (p[pi].iOrder < pm[pmi].iOrder) {
			fprintf(fp,"0\n");
			++pi;
			}
		assert(p[pi].iOrder == pm[pmi].iOrder);
		dy = pm[pmi].r[1] - p[pi].r[1];
		if (dy > hy) dy -= 2*hy;
		if (dy <= -hy) dy += 2*hy;
		fprintf(fp,"%g\n",dy);
		}
	for (;pi<kd->nParticles;++pi) fprintf(fp,"0\n");
	for (pmi=0,pi=0;pmi<kd->nMove;++pmi,++pi) {
		while (p[pi].iOrder < pm[pmi].iOrder) {
			fprintf(fp,"0\n");
			++pi;
			}
		assert(p[pi].iOrder == pm[pmi].iOrder);
		dz = pm[pmi].r[2] - p[pi].r[2];
		if (dz > hz) dz -= 2*hz;
		if (dz <= -hz) dz += 2*hz;
		fprintf(fp,"%g\n",dz);
		}
	for (;pi<kd->nParticles;++pi) fprintf(fp,"0\n");
	fclose(fp);
	}




void kdFinish(KD kd)
{
	if (kd->pMove) free(kd->pMove);
	if (kd->pInit) free(kd->pInit);
	if (kd->pGroup) free(kd->pGroup);
	if (kd->kdNodes) free(kd->kdNodes);
	if (kd->piGroup) free(kd->piGroup);
	free(kd);
	}

void kdOutTemperature(KD kd, char *pszFile)
{
  FILE *fp;
  int i;
  PMOVE *pm;


  if (kd->bOutDiag) puts(">> kdOutTemperature()");
  fflush(stdout);
  
  if (pszFile) {
    fp = fopen(pszFile,"w");
    assert(fp != NULL);
  }
  else {
    fp = stdout;
  }

  // Order(kd);

  fprintf(fp, "%d\n", kd->nMove);

  pm = kd->pMove;

  for (i=0;i<kd->nMove;i++) {
    fprintf(fp, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",	  \
	    pm[i].r[0], pm[i].r[1], pm[i].r[2], \
	    pm[i].v[0], pm[i].v[1], pm[i].v[2], \
	    pm[i].fMass, pm[i].fDensity, pm[i].fTemp, pm[i].fMetals);

  }

  fclose(fp);
  
  if (kd->bOutDiag) puts("<< kdOutTemperature()");
  fflush(stdout);
}


void kdSetIord(KD kdg, KD kdd, KD kds, float radius) {
  
  PMOVE *pg, *ps;
  PINIT *pd, *ps_old;

  int i, oldstar=0;

  pg = kdg->pMove;
  pd = kdd->pInit;
  ps = kds->pMove;
  ps_old = kds->pInit;

  for (i=0;i<kdg->nMove;i++)     
    pg[i].iOrder = i;

  for (i=0;i<kdd->nInitActive;i++)
    pd[i].iOrder = i + kdg->nMove;

  for(i=kds->nInitActive;i<kds->nParticles;i++)
    {
      if((ps_old[i].r[0]*ps_old[i].r[0] + 
	  ps_old[i].r[1]*ps_old[i].r[1] + 
	  ps_old[i].r[2]*ps_old[i].r[2]) < radius*radius)
	{
	  ps_old[i].iOrder = i + kdg->nMove + kdd->nInitActive;
	  oldstar++;
	}
    }
  for (i=0;i<kds->nMove;i++)
    ps[i].iOrder = i + kdg->nMove + kdd->nInitActive + oldstar;

}
  
  

void kdWriteTipsyStd(KD kdg, KD kdd, KD kds, char *outFile, float radius, int nSplitting) {
  
  struct gas_particle gas;
  struct dark_particle dark;
  struct star_particle star;
  struct dump header;
  float xcurr, ycurr, zcurr, comr[3], comv[3];
  PINIT *pg, *ps, *pd;

  int i, j, pi, oldgas=0, oldstar=0, olddark=0;
  XDR xdrs;
  FILE *fp;

  if (kdg->bOutDiag) puts(">> kdWriteTipsyStd()");
  fflush(stdout);
  
  if (outFile) {
    fp = fopen(outFile,"w");
    assert(fp != NULL);
  }
  else {
    fp = stdout;
  }
  
  setvbuf(stdin, NULL, _IOFBF, 32*4096);
  setvbuf(stdout, NULL, _IOFBF, 32*4096);
  xdrstdio_create(&xdrs, fp, XDR_ENCODE);

  /* 
  ** Initialize the header and write to the file
  */

  header.time = (double)kdg->fTime;
  header.ndim = 3;

  
  pg = kdg->pInit;
  
  // how many old non-split dark particles
  pd = kdd->pInit;
  /* for (pi=0;pi<kdd->nParticles;pi++)  */
/*     { */
/*       if(((pd[pi].r[0]*pd[pi].r[0] + pd[pi].r[1]*pd[pi].r[1] + pd[pi].r[2]*pd[pi].r[2]) */
/* 	  < radius*radius)) */
/* 	olddark++; */
/*     } */

  olddark = kdd->nInitActive;

  fprintf(stderr, "olddark in writetipsy = %d\n", olddark);

  // how many old non-split star particles
  ps = kds->pInit;
  for (pi=0;pi<kds->nParticles;pi++) 
    {
      if(((ps[pi].r[0]*ps[pi].r[0] + ps[pi].r[1]*ps[pi].r[1] + ps[pi].r[2]*ps[pi].r[2])
	  < radius*radius) && (ps[pi].fTimeForm < 0.0))
	oldstar++;
    }
  
  header.nsph = kdg->nMove;
  header.ndark = olddark;
  header.nstar = kds->nMove + oldstar;
  header.nbodies = header.nsph + header.nstar + header.ndark;

  xdr_header(&xdrs, &header);
  
  /* 
  ** loop over gas particles - only the split particles are used,
  ** so use only the pMove array
  */

  
  for (i=0;i<kdg->nMove;i++) {
    for (j=0;j<3;j++) {
      gas.pos[j] = kdg->pMove[i].r[j];
      gas.vel[j] = kdg->pMove[i].v[j];
    }

    gas.mass = kdg->pMove[i].fMass;
    gas.rho = kdg->pMove[i].fDensity;
    gas.temp = kdg->pMove[i].fTemp;
    assert(kdg->pMove[i].fSoft != 0.0);
    gas.hsmooth = kdg->pMove[i].fSoft;
    
    gas.metals = kdg->pMove[i].fMetals;
    gas.phi = 0.0;

    xdr_gas(&xdrs, &gas);
  }
      
  /* 
  ** loop over dark particles - none of the dark particles are split,
  ** so use only the pInit array
  */

  

  for (i=0;i<kdd->nInitActive;i++) {
    xcurr = pd[i].r[0];
    ycurr = pd[i].r[1];
    zcurr = pd[i].r[2];
    
    if((xcurr*xcurr + ycurr*ycurr + zcurr*zcurr) <= radius*radius)
      {
	for (j=0;j<3;j++) {
	  dark.pos[j] = kdd->pInit[i].r[j];
	  dark.vel[j] = kdd->pInit[i].v[j];
	}

	dark.mass = kdd->pInit[i].fMass;
	dark.eps = kdd->pInit[i].fSoft;
	dark.phi = 0.0;
	
	xdr_dark(&xdrs, &dark);
      }
    else fprintf(stderr, "erm... DARK PARTICLE OUTSIDE RADIUS! \n", (xcurr*xcurr + ycurr*ycurr + zcurr*zcurr));
  }

  /* 
  ** loop over star particles - the old particles aren't relevant for the splitting
  ** so take those from the pInit array, but take the split particles from the pMove array
  */

  
  for (i=kds->nInitActive;i<kds->nParticles;i++)
    {
      xcurr = ps[i].r[0];
      ycurr = ps[i].r[1];
      zcurr = ps[i].r[2];
    
    if((xcurr*xcurr + ycurr*ycurr + zcurr*zcurr) < radius*radius)
      {
	for (j=0;j<3;j++) {
	  star.pos[j] = kds->pInit[i].r[j];
	  star.vel[j] = kds->pInit[i].v[j];
	}

	star.mass = kds->pInit[i].fMass;
	star.tform = kds->pInit[i].fTimeForm;
	star.metals = kds->pInit[i].fMetals;
	star.eps = kds->pInit[i].fSoft;
	star.phi = 0.0;

	xdr_star(&xdrs, &star);
      }
    }

  for (i=0;i<kds->nMove;i++) {
    for (j=0;j<3;j++) {
      star.pos[j] = kds->pMove[i].r[j];
      star.vel[j] = kds->pMove[i].v[j];
    }

    star.mass = kds->pMove[i].fMass;
    star.tform = kds->pMove[i].fTimeForm;
    star.metals = kds->pMove[i].fMetals;
    star.eps = kds->pMove[i].fSoft;
    star.phi = 0.0;

    xdr_star(&xdrs, &star);
  }
 
  xdr_destroy(&xdrs);
     
  if (kdg->bOutDiag) puts("<< kdWriteTipsyStd()");
  fflush(stdout);
  

}


void kdWriteTipsyCheckpoint(KD kdg, KD kdd, KD kds, CHK_HEADER header, 
			    char *outFile, float radius) {

  CHK_PART *cp;
  PMOVE pm;
  PINIT *pd, *ps;
  int i,j, pi, olddark=0, oldstar=0, cp_ind=0;
  FILE *fp, *FDL;
  char c[10000];

  if (kdg->bOutDiag) puts(">> kdWriteTipsyCheckpoint()");
  fflush(stdout);
  
  // how many old non-split dark particles
  pd = kdd->pInit;
  for (pi=0;pi<kdd->nParticles;pi++) 
    {
      if(((pd[pi].r[0]*pd[pi].r[0] + pd[pi].r[1]*pd[pi].r[1] + pd[pi].r[2]*pd[pi].r[2])
	  < radius*radius))
	olddark++;
    }
  
  // how many old non-split star particles
  ps = kds->pInit;
  for (pi=0;pi<kds->nParticles;pi++) 
    {
      if(((ps[pi].r[0]*ps[pi].r[0] + ps[pi].r[1]*ps[pi].r[1] + ps[pi].r[2]*ps[pi].r[2])
	  < radius*radius) && ps[pi].fTimeForm < 0.0)
	oldstar++;
    }
  
  header.number_of_gas_particles = kdg->nMove;
  header.number_of_dark_particles = olddark;
  header.number_of_star_particles = kds->nMove+oldstar;

  fprintf(stderr, "olddark = %d\tnInitDark = %d\n", olddark, kdd->nInitActive);

  header.number_of_particles = 
    header.number_of_gas_particles +
    header.number_of_dark_particles + 
    header.number_of_star_particles;

  header.max_order = kds->pMove[kds->nMove-1].iOrder;
  header.max_order_gas = kdg->pMove[header.number_of_gas_particles-1].iOrder;
  header.max_order_dark = kdd->pInit[header.number_of_dark_particles-1].iOrder;
  
  fprintf(stderr, "max order = %d\nmax gas order = %d\nmax dark order = %d\n", 
	  header.max_order, header.max_order_gas, header.max_order_dark);


  cp = (CHK_PART *)malloc(sizeof(CHK_PART)*header.number_of_particles);

  fprintf(stderr, "header.number_of_particles = %d\n", header.number_of_particles);

  fprintf(stderr, "size of cp = %d, size of chkpart = %d size of chkheader = %d\n", 
	  sizeof(cp), sizeof(CHK_PART), sizeof(CHK_HEADER));
  
  /* 
  ** Write the gas particles
  */

  for (i=0;i<header.number_of_gas_particles;i++)
    {
      if (i<10) fprintf(stderr, "%f\n", kdg->pMove[i].fMass);
      assert(kdg->pMove[i].fMass > 0);
      kdWriteMovedParticle(cp, cp_ind, kdg->pMove[i]);
      cp_ind++;
    }

  /*
  ** Write the dark particles
  */


  for (i=0;i<header.number_of_dark_particles;i++)
    {
      kdWriteInitParticle(cp, cp_ind, kdd->pInit[i]);
      cp_ind++;
    }

  /* 
  ** Write the star particles - need to take some out of the init array
  ** and some out of the move array
  */
  
  for (i=kds->nInitActive;i<kds->nParticles;i++)
    {
      if((ps[i].r[0]*ps[i].r[0] + ps[i].r[1]*ps[i].r[1] + ps[i].r[2]*ps[i].r[2]) 
	 < radius*radius)
      {
	assert(kds->pInit[i].fTimeForm < 0.0);
	kdWriteInitParticle(cp, cp_ind, kds->pInit[i]);
	cp_ind++;
      }
    }

  for (i=0;i<kds->nMove;i++) 
    {
      kdWriteMovedParticle(cp, cp_ind, kds->pMove[i]);
      cp_ind++;
    }
  

  if(getenv("PKDGRAV_CHECKPOINT_FDL") == NULL) {
    fprintf(stderr,"PKDGRAV_CHECKPOINT_FDL environment variable not set\n");
    return;
  }
  else  FDL = fopen(getenv("PKDGRAV_CHECKPOINT_FDL"),"r");
  if(FDL == NULL) {
    fprintf(stderr,"error opening %s\n",getenv("PKDGRAV_CHECKPOINT_FDL"));
    return;
  }
  
  fp = fopen(outFile, "w");
  
  while(!feof(FDL)) {
    fread(&c[0],sizeof(char),1,FDL);
    fprintf(fp,"%c",c[0]);
  }

  fseek(fp,-2,SEEK_CUR);

  fwrite(&header,sizeof(struct chkHeader),1,fp);
  
  fwrite(cp,sizeof(struct chk_particle)*header.number_of_particles,1,fp);
  
  fclose(fp);
  fclose(FDL);

  printf("N: %d Ngas: %d Ndark: %d Nstar: %d\n", 
	 header.number_of_particles,
	 header.number_of_gas_particles,
	 header.number_of_dark_particles,
	 header.number_of_star_particles);
  fflush(stdout);
  

  if (kdg->bOutDiag) puts("<< kdWriteTipsyCheckpoint()");
	fflush(stdout);

  return;
}


void kdWriteMovedParticle(CHK_PART *cp, int i, PMOVE pm)
{
  
  int j;

  cp[i].iOrder = pm.iOrder;
  cp[i].fMass = pm.fMass;
  cp[i].fSoft = pm.fSoft;
  assert(cp[i].fSoft != 0.0);
  for(j=0;j<3;j++) { cp[i].r[j] = pm.r[j]; cp[i].v[j] = pm.v[j]; }
  cp[i].u = pm.u;
  cp[i].fMetals = pm.fMetals;
  cp[i].CoolParticle = pm.CoolParticle;
  cp[i].fTimeCoolIsOffUntil = pm.fTimeCoolIsOffUntil;
  cp[i].fTimeForm = pm.fTimeForm;
  cp[i].fMassForm = pm.fMassForm;
  cp[i].fNSN = pm.fNSN;
  cp[i].fMFracOxygen = pm.fMFracOxygen;
  cp[i].fMFracIron = pm.fMFracIron;
  cp[i].iGasOrder = pm.iGasOrder;
}

void kdWriteInitParticle(CHK_PART *cp, int i, PINIT pi)
{
  
  int j;

  cp[i].iOrder = pi.iOrder;
  cp[i].fMass = pi.fMass;
  cp[i].fSoft = pi.fSoft;
  assert(cp[i].fSoft != 0.0);
  for(j=0;j<3;j++) { cp[i].r[j] = pi.r[j]; cp[i].v[j] = pi.v[j]; }
  cp[i].u = pi.u;
  cp[i].fMetals = pi.fMetals;
  cp[i].CoolParticle = pi.CoolParticle;
  cp[i].fTimeCoolIsOffUntil = pi.fTimeCoolIsOffUntil;
  cp[i].fTimeForm = pi.fTimeForm;
  cp[i].fMassForm = pi.fMassForm;
  cp[i].fNSN = pi.fNSN;
  cp[i].fMFracOxygen = pi.fMFracOxygen;
  cp[i].fMFracIron = pi.fMFracIron;
  cp[i].iGasOrder = pi.iGasOrder;
}

void kdSmbhSoft(KD kd, int nSplitting)
{
  
  int i;
  double max_mass;
  
  // smbh particles are ones with the largest mass
  max_mass = kd->pInit[0].fMass;
  for (i=0;i<kd->nParticles;i++) {
    if (max_mass < kd->pInit[i].fMass) {
      max_mass = kd->pInit[i].fMass;
    }
  }
  
  fprintf(stderr, "MAX MASS = %f\n", max_mass);
  for (i=0;i<kd->nParticles;i++) 
    if(kd->pInit[i].fMass == max_mass) 
      kd->pInit[i].fSoft = kd->pInit[i].fSoft/cbrt((double)nSplitting);
    
}
