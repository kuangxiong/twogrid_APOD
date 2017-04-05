/*************************************************************************
	> File Name: APODfun.h
	> Author:kuangxiong 
	> Mail: kuangxiong@lsec.cc.ac.cn
	> Created Time: 2016年10月21日 星期五 21时26分17秒
 ************************************************************************/
#ifndef  __APODFUN_H_
#define __APODFUN_H_
#endif

#include<stdio.h>
#include"mpi.h"
#include"math.h"
#include"phg.h"

static FLOAT PODT0 = 8.0;
static FLOAT APODT0 = 4.0;
static FLOAT CLT = 2.0;
static INT numt = 30; /* number of coarse grid dof for updata pod basis  */
static FLOAT T = 50;
static FLOAT Tolerr = 0.05;
static INT RefineN = 10;
static FLOAT gamma1 = 0.99;
static FLOAT gamma2 = 0.99;
static FLOAT gamma3 = 0.999;
static FLOAT stime1 = 0.1;
static FLOAT stime = 0.01;
static FLOAT D0 = 0.5;
static FLOAT crtime = 0.0;    /* current time. Modified by program */
static FLOAT ptime = 0.0;     /* previous time. */

static const char *fn="meshdata/testcube6.dat";

static void     /* rhs: f */
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
//   *value = -Cos(y)-Sin(y)*Cos(crtime);
//   *value = -Cos(y)-Sin(y)*Cos(crtime) + 2*exp(-(crtime-100)*(crtime -100)/0.05);
     *value = (1 + 3*crtime)*Cos(x)*Cos(y)*Cos(z);
//	   *value = y + cos(0.1*crtime) + 10*exp(-(crtime -50)*(crtime -50)/0.05);
//     *value = y*Cos(0.5*crtime) + 5*exp(-(crtime-200)*(crtime-200)/0.05);
//       *value = -2*sin(2*crtime)*(x+y+z)*0.5*exp(4*Cos(2*crtime));
//       *value = 0.1*crtime*Cos(x)*Cos(y)*Cos(z);
//         *value = y*0.1*crtime*Cos(x);
//        *value = y*Cos(0.5*crtime) + 5*exp(-(crtime-100)*(crtime-100)/0.05);
//       *value = Sin(crtime)*Sin(y) -Cos(y) - Sin(y)*Cos(crtime);
}
static void     /* rhs: f */
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
//   *value = -Cos(y)-Sin(y)*Cos(crtime);
//     *value = crtime*Cos(x)*Cos(y)*Cos(z);
}

static void 
func_B1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
	*value = Cos(x);
//    *value = 0;
}

static void
func_B2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
	*value = Cos(y) ;
//    *value = 0;
}

static void 
func_B3(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
	*value = Cos(z);
//    *value = 0;
}

static void     /* rhs: f */
func_C1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
//   *value = Sin(y);
    *value = 0;
}
static void     /* rhs: f */
func_C2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
//   *value = Sin(z);
    *value = 0;
}

static void     /* rhs: f */
func_C3(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
//   *value = Sin(x);
    *value = 0;
}

static void
func_u0(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
	*value = 0;
//   *value = -Cos(y)- Sin(y);
}
static FLOAT
eta_t(FLOAT t)
{   FLOAT z;
	z = cos(t);
	return z;
}

static FLOAT 
f_tmp(FLOAT t)
{  FLOAT z;
//  if(t <=500)
		z = cos(t);
//	else 
//		z = cos(t) + sin(t);
    return z;
}

///*********************computing D11*******************/
static FLOAT 
GetD11(DOF *u_h)
{
     ELEMENT *e;
	 DOF *grad_u = phgDofGradient(u_h, NULL, NULL, NULL);

	 FLOAT  norm;
#if USE_MPI
	 FLOAT a, b;
#endif

	 norm = 0;
     ForAllElements(u_h->g, e)
	 {
        norm += phgQuadDofDotDof(e, grad_u, grad_u, QUAD_DEFAULT);
	 }
#if USE_MPI
	 a = norm;
	 MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_SUM, u_h->g->comm);
	 norm = b;
#endif
	 phgDofFree(&grad_u);
	 
	 return norm;
}

static FLOAT
estimate_space_error3(FLOAT time, MAT *G1, MAT *G2, MAT *C, VEC *b, DOF *u_h, DOF *u_p)
{
    MAT *M;
	VEC *bb, *btmp,*btmp1,  *b1;
	MAP *map = b->map;
	GRID *g = u_h->g;
    ELEMENT *e;
    FLOAT eta, h, d;

	bb = phgMapCreateVec(map, 1);
	b1 = phgMapCreateVec(map, 1);
	btmp = phgMapCreateVec(map, 1);
	btmp1 = phgMapCreateVec(map, 1);
	M = phgMapCreateMat(map, map);
	phgMapDofArraysToVec(map, 1, FALSE, &btmp, &u_p, NULL);

	phgMatAXPBY(1.0, G1, 0.0, &M);
	phgMatAXPBY(1.0*Cos(time), G2, 1.0, &M);
    phgMatVec(0, 1.0, C, btmp, 0.0, &bb);
	phgVecAXPBY(1.0, b, 1.0, &bb);
    
	phgMapDofArraysToVec(map, 1, FALSE, &btmp1, &u_h, NULL);
	
	phgMatVec(0, 1.0, M, btmp1, 0.0, &b1);
	phgVecAXPBY(-1.0, b1, 1.0, &bb);

    phgMatDestroy(&M);

	//return phgVecNorm2(b, 0, NULL)/phgVecNorm2(btmp1, 0, NULL);
    return phgVecNorm2(bb, 0, NULL);
}
static FLOAT
estimate_space_error(FLOAT time, MAT *G1, MAT *G2, MAT *C, VEC *g1, VEC *g2, DOF *u_h, DOF *u_p)
{
    MAT *M;
	VEC *bb, *btmp,*btmp1,  *b1;
	MAP *map = g1->map;
	GRID *g = u_h->g;
    ELEMENT *e;
    FLOAT eta, h, d;

	bb = phgMapCreateVec(map, 1);
	b1 = phgMapCreateVec(map, 1);
	btmp = phgMapCreateVec(map, 1);
	btmp1 = phgMapCreateVec(map, 1);
	M = phgMapCreateMat(map, map);
	phgMapDofArraysToVec(map, 1, FALSE, &btmp, &u_p, NULL);

	phgMatAXPBY(1.0, G1, 0.0, &M);
	phgMatAXPBY(1.0*Cos(time), G2, 1.0, &M);
    phgMatVec(0, 1.0, C, btmp, 0.0, &bb);
	phgVecAXPBY(-1.0, g1, 1.0, &bb);
	phgVecAXPBY(-1.0*f_tmp(time), g2, 1.0, &bb);
    
	phgMapDofArraysToVec(map, 1, FALSE, &btmp1, &u_h, NULL);
	
	phgMatVec(0, 1.0, M, btmp1, 0.0, &b1);
	phgVecAXPBY(-1.0, b1, 1.0, &bb);

    phgMatDestroy(&M);

	//return phgVecNorm2(b, 0, NULL)/phgVecNorm2(btmp1, 0, NULL);
    return phgVecNorm2(bb, 0, NULL);
}
static void
phgGetTruncationMat(DOF *u_h, FLOAT *TrM, FLOAT *U, int N, MAP *map, INT M)
{
	INT i, j, k, maxnlocal, descU[9], descM[9], nfdom, ONE = 1;
	FLOAT *adjustU;

    nfdom = DofGetDataCountGlobal(u_h);
	MPI_Reduce(&map->nlocal, &maxnlocal, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxnlocal, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
 //   adjustU = (FLOAT *)malloc((map->nlocal * N) * sizeof(FLOAT));
 //   for(i=0; i< map->nlocal; i++)
//		for(j = 0; j< N; j++)
//			adjustU[i+j * map->nlocal] = U[i + j * maxnlocal];
	descU[0] = 1; descU[1] = 0; descU[2] = M; descU[3]= N; descU[4] = 1;
	descU[5] = 1; descU[6] = 0 ; descU[7] = 0; descU[8] = maxnlocal; 
	descM[0] = 1; descM[1] = 0; descM[2] = M; descM[3]= M; descM[4] = 1;
	descM[5] = 1; descM[6] = 0; descM[7] = 0; descM[8] = maxnlocal; 
    double alpha = 1.0;
	double beta = 0;
	pdgemm_("N", "T", &M, &M ,&N, &alpha, U, &ONE, &ONE, descU, U, &ONE, &ONE, descU, &beta, TrM, &ONE, &ONE, descM);
 //   free(adjustU);
}

static FLOAT
estimate_TGAPOD_error1(DOF *u_h, FLOAT *U, MAP *map, INT N)
{
    DOF *error, *u_p;
	FLOAT  a, x[N];
	GRID *g = map->dofs[0]->g;
	VEC *V, *V1, *V2;
    INT  i, maxnlocal;
	
	MPI_Reduce(&map->nlocal, &maxnlocal, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxnlocal, 1, MPI_INT, 0, MPI_COMM_WORLD);

	V = phgMapCreateVec(map, 1);
    V1 = phgMapCreateVec(map, 1);
    V2 = phgMapCreateVec(map, 1);
    
	phgMapDofArraysToVec(map, 1, FALSE, &V, &u_h, NULL);
	for(i=0; i< N; i++)
	{	
		memcpy(V1->data, &U[i*maxnlocal], map->nlocal*sizeof(FLOAT));
		x[i] = phgVecDot(V1, 0, V, 0, NULL);
	}
	memcpy(V1->data, &U[0], map->nlocal*sizeof(FLOAT));
    phgVecAXPBY(x[0], V1, 0, &V2);
	for(i=1; i< N ; i++)
    {	
		memcpy(V1->data, &U[i*maxnlocal], map->nlocal*sizeof(FLOAT));
        phgVecAXPBY(x[i], V1, 1.0, &V2);
    }
	u_p = phgDofNew(g, DOF_P1, 1, "u_p", DofInterpolation);
    error = phgDofNew(g, DOF_P1, 1, "error", DofInterpolation);
    phgMapVecToDofArrays(map, V2, FALSE, &u_p, NULL);
	phgDofAXPBY(1.0, u_h, 0, &error);
    phgDofAXPBY(-1.0, u_p, 1.0, &error);
    a = phgDofNormL2(error);
	phgDofFree(&error);
	phgDofFree(&u_p);
    return a;
}

static FLOAT
estimate_TGAPOD_error(DOF *u_h, FLOAT *TrM, MAP *map, INT N, INT M)
{
    DOF *error, *u_p;
	FLOAT  a, *adjustU, *x, *y, *tmpy;
	GRID *g = map->dofs[0]->g;
	VEC *V;
    INT  descA[9], descX[9], descY[9], ONE = 1, nfdom, maxnlocal, i;
	
	MPI_Reduce(&map->nlocal, &maxnlocal, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxnlocal, 1, MPI_INT, 0, MPI_COMM_WORLD);
    V = phgMapCreateVec(map, 1);
    u_p = phgDofNew(g, DOF_P1, 1, "u_p", DofInterpolation);
    error = phgDofNew(g, DOF_P1, 1, "error", DofInterpolation);
    x = malloc(maxnlocal*sizeof(double));
	y = malloc(maxnlocal*sizeof(double));
	tmpy = malloc(map->nlocal*sizeof(double));

	phgMapDofArraysToVec(map, 1, FALSE, &V, &u_h, NULL);
	
	for(i=0; i< maxnlocal; i++)
	{
	    if(i< map->nlocal)
			x[i] = V->data[i];
		else 
			x[i] = 0;
	}
	phgMapDofToLocalData(map, 1, &u_h, x);
    descA[0] =1; descA[1] = 0; descA[2] = M; descA[3] = M; descA[4] = 1;
	descA[5]= 1; descA[6] =0 ; descA[7] = 0 ;descA[8] = maxnlocal;
    
	descX[0] =1; descX[1] = 0; descX[2] = M ; descX[3] = 1 ; descX[4] = 1;
	descX[5]= 1; descX[6] =0 ; descX[7] = 0 ;descX[8] = maxnlocal;

    descY[0] =1; descY[1] = 0; descY[2] = M; descY[3] = 1; descY[4] = 1;
	descY[5]= 1; descY[6] = 0; descY[7] =  0;descY[8] = maxnlocal;
    double alpha = 1.0; 
	double beta = 0.0;
	pdgemv_("N", &M, &M, &alpha, TrM, &ONE, &ONE, descA, x, &ONE, &ONE, descX, &ONE, &beta, y, &ONE, &ONE, descY,&ONE);
    for(i=0;i< map->nlocal; i++)
		tmpy[i] = y[i];
	phgMapLocalDataToDof(map, 1, &u_p, tmpy);
	phgDofAXPBY(1.0, u_h, 0, &error);
    phgDofAXPBY(-1.0, u_p, 1.0, &error);
    a = phgDofNormL2(error);
	phgDofFree(&error);
	phgDofFree(&u_p);
    return a;
}

static FLOAT
estimate_space_error1(DOF *u_h, DOF *u_p, DOF *f_h, DOF *error)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *tmp, *jump;
    FLOAT eta, h;
    int i, ONE = 1;
  
    tmp = phgDofGradient(u_h, NULL, NULL, "tmp");
    jump = phgQuadFaceJump(tmp, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    phgDofFree(&tmp);

    eta = 1.0 / stime;
    tmp = phgDofCopy(f_h, NULL, f_h->type, NULL);
    phgDofAXPY(eta, u_p, &tmp);
    eta *= -1.0;
    phgDofAXPY(eta, u_h, &tmp);

    ForAllElements(g, e){
        eta = 0.0;
        for (i = 0; i < NFace; i++) {
            if (e->bound_type[i] & (DIRICHLET | NEUMANN))
                continue;    /* boundary face */
                h = phgGeomGetFaceDiameter(g, e, i);
                eta +=  (*DofFaceData(jump, e->faces[i])) * h;
        }
        h = phgGeomGetDiameter(g, e);
    //    eta += eta * 0.5  + h * h * phgQuadDofDotDof(e, tmp, tmp, QUAD_DEFAULT);
        eta =  h * h * phgQuadDofDotDof(e, tmp, tmp, QUAD_DEFAULT);

        *DofElementData(error, e->index) = eta;
    }
    phgDofFree(&tmp);
    phgDofFree(&jump);
    return phgDofNormInftyVec(error);
}

static FLOAT
estimate_time_error(DOF *u_h, DOF *u_p)
{
    DOF *tmp, *grad;
    FLOAT ss;

    tmp = phgDofCopy(u_h, NULL, DOF_DEFAULT, "tmp2");
    phgDofAXPY( -1.0, u_p, &tmp);
    grad = phgDofGradient(tmp, NULL, NULL, "grad");
    ss = phgDofNormL2(grad);

    phgDofFree(&tmp);
    phgDofFree(&grad);

    return ss * ss * (FLOAT)(1.L/3.L);
}

/* data oscillation */

static void
phgFEMLoopTime3(FLOAT ctime, DOF *u_h, DOF *u_p, MAT *C, MAT *G1, MAT *G2, VEC *bb)
{     
	MAP *map = bb->map;
	VEC *btmp, *btmp1, *b;
	MAT *A;
    SOLVER *solver;
    
    btmp = phgMapCreateVec(map, 1);

    phgMapDofArraysToVec(map, 1, FALSE, &btmp, &u_p, NULL);
    b = phgMapCreateVec(map, 1);
	phgVecDisassemble(b);
	phgMatVec(0, 1.0, C, btmp, 0.0, &b);
	phgVecAXPBY(1.0, bb, 1.0 , &b); 

	A = phgMapCreateMat(map,map);
	phgMatAXPBY(1.0, G1, 0, &A);
	phgMatAXPBY(eta_t(ctime), G2, 1.0, &A);
		
    solver = phgSolverCreate(SOLVER_GMRES, u_h, NULL);
	solver->mat = A;
	solver->rhs = b;
    phgSolverSolve(solver, TRUE, u_h, NULL);
	phgSolverDestroy(&solver);
}
static void
phgFEMLoopTime(FLOAT ctime, DOF *u_h, DOF *u_p, MAT *C, MAT *G1, MAT *G2, VEC *g1, VEC *g2)
{     
	MAP *map = g1->map;
	VEC *btmp, *btmp1, *b;
	MAT *A;
    SOLVER *solver;
    
    btmp = phgMapCreateVec(map, 1);
    btmp1 = phgMapCreateVec(map, 1);

    phgMapDofArraysToVec(map, 1, FALSE, &btmp, &u_p, NULL);
    b = phgMapCreateVec(map, 1);
	phgVecDisassemble(b);
	phgMatVec(0, 1.0, C, btmp, 0.0, &b);
		
	phgVecAXPBY(-1.0, g1, 0.0 , &btmp1); 
	phgVecAXPBY(-1.0*f_tmp(ctime), g2, 1.0, &btmp1);

	phgVecAXPBY(1.0, btmp1, 1.0, &b);

	A = phgMapCreateMat(map,map);
	phgMatAXPBY(1.0, G1, 0, &A);
	phgMatAXPBY(eta_t(ctime), G2, 1.0, &A);
		
    solver = phgSolverCreate(SOLVER_GMRES, u_h, NULL);
	solver->mat = A;
	solver->rhs = b;
    phgSolverSolve(solver, TRUE, u_h, NULL);
	phgSolverDestroy(&solver);

}
static void
phgRefineFEMLoopTime(FLOAT timestp, GRID *g, INT TotalN,  DOF **uhs, INT rank, char *snD1)
{
    INT i, j, k, flag = 0;
	FLOAT err;
	MAT *AA, *GG1, *GG2, *CC;
	VEC *bb, *bbtmp, *bbtmp1;
	MAP *mmap;
	SOLVER *ssolver;
	ELEMENT *ee;
	DOF *uu_h, **uuhs, *ff_h;  /* numerical solution at time t_n */
	DOF *uu_p, *ggrad_u;  /* numerical solution at time t_{n-1} */
	DOF  *BB1, *BB2, *BB3, *CC1, *CC2, *CC3, *eerror;
    
	uuhs = phgAlloc(TotalN * sizeof(*uhs));
	uu_p = phgDofNew(g, DOF_P1, 1, "uu_p", DofInterpolation);
	ff_h = phgDofNew(g, DOF_P1, 1, "ff_h", DofInterpolation);
	phgDofSetDataByFunction(uu_p, func_u0);
	uu_h = phgDofNew(g, DOF_P1, 1, "uu_h", DofInterpolation);
	eerror = phgDofNew(g, DOF_P1, 1, "eerror", DofInterpolation);	
	uuhs[0] = phgDofCopy(uu_p, NULL, DOF_P1, "uu_p1");
	uuhs[0] ->userfunc = DofInterpolation;
	k = 0;
	BB1=phgDofNew(g, DOF_P1, 1, "BB1", func_B1);
	BB2=phgDofNew(g, DOF_P1, 1, "BB2", func_B2);
	BB3=phgDofNew(g, DOF_P1, 1, "BB3", func_B3);
	CC1=phgDofNew(g, DOF_P1, 1, "CC1", func_C1);
	CC2=phgDofNew(g, DOF_P1, 1, "CC2", func_C2);
	CC3=phgDofNew(g, DOF_P1, 1, "CC3", func_C3);

	mmap = phgMapCreate(uu_h, NULL);
	GG1 = phgMapCreateMat(mmap, mmap);
	GG2 = phgMapCreateMat(mmap, mmap);
	CC = phgMapCreateMat(mmap, mmap);

	bbtmp = phgMapCreateVec(mmap, 1);
	phgVecDisassemble(bbtmp);

	build_rhs_Mat(CC, mmap, uu_h);
	build_stiff_Mat1(timestp, D0, GG1, mmap, BB1, BB2, BB3, uu_h);
	build_stiff_Mat2(timestp, GG2, mmap, CC1, CC2, CC3, uu_h);
	crtime = 0;
	while(crtime < T - 1e-8)
	{
		/********************************************************************/	
		ptime = crtime;
		crtime += timestp;
		phgPrintf("\n/********* start new time layer *************/\n");
		phgPrintf("current time layer: [%lf, %lf]\n", (double)ptime, (double)crtime);
		flag++;
		if (flag > 1)
		{   /* update u_p */
			phgDofFree(&uu_p);
			uu_p = phgDofCopy(uu_h, NULL, DOF_P1, "uu_p");
			uu_p->userfunc = DofInterpolation;
		}
		phgMapDofArraysToVec(mmap, 1, FALSE, &bbtmp, &uu_p, NULL);
		bb = phgMapCreateVec(mmap, 1);
		bbtmp1 = phgMapCreateVec(mmap, 1);
		phgVecDisassemble(bb);
		phgVecDisassemble(bbtmp1);

		phgMatVec(0, 1.0, CC, bbtmp, 0.0, &bb);	
		phgDofSetDataByFunction(ff_h, func_f);

		build_rhs_Vec(timestp, bbtmp1, ff_h, mmap, uu_h);
		phgVecAXPBY(1.0, bbtmp1, 1.0, &bb);

		AA = phgMapCreateMat(mmap, mmap);
		phgMatAXPBY(1.0, GG1, 0.0, &AA);
		phgMatAXPBY(eta_t(crtime), GG2, 1.0, &AA);

		ssolver = phgSolverCreate(SOLVER_GMRES, uu_h, NULL);
		ssolver->mat = AA;
		ssolver->rhs = bb;
		phgSolverSolve(ssolver, TRUE, uu_h, NULL);
		phgSolverDestroy(&ssolver);
        k++;
		uuhs[k] = phgDofCopy(uu_h, NULL, DOF_P1, "uu_p1");
		uuhs[k] ->userfunc = DofInterpolation;
		phgPrintf("\ntime step: %lg\n\n", (double)timestp);
	}

	for(i=0; i< TotalN; i++)
	{    
		phgDofAXPBY(1.0, uhs[i], 0, &eerror);
		phgDofAXPBY(-1.0, uuhs[i], 1.0, &eerror);
		err = phgDofNormL2(eerror);
	 //   phgPrintf( "%f\t%f\n", i*timestp, err);
		if(rank==0)
	    	fprintf(snD1, "%f\t%f\n", i*timestp, err);
	}
	for(i=0 ;i< TotalN; i++)
		phgDofFree(uuhs+i);
	phgFree(uuhs);
	phgDofFree(&eerror);
	phgDofFree(&BB1);
	phgDofFree(&BB2);
	phgDofFree(&BB3);
	phgDofFree(&ff_h);

	phgDofFree(&CC1);
	phgDofFree(&CC2);
	phgDofFree(&CC3);
	phgDofFree(&uu_h);
	phgDofFree(&uu_p);
}
static void
phgRefineFEMLoopTime1(GRID *g, INT TotalN,  DOF **uhs, INT rank, char *snD1)
{
    INT i, j, k, flag = 0;
	FLOAT err;
	MAT *AA, *GG1, *GG2, *CC;
	VEC *bb, *bbtmp, *bbtmp1;
	MAP *mmap;
	SOLVER *ssolver;
	ELEMENT *ee;
	DOF *uu_h, **uuhs, *ff_h;  /* numerical solution at time t_n */
	DOF *uu_p, *ggrad_u;  /* numerical solution at time t_{n-1} */
	DOF  *BB1, *BB2, *BB3, *CC1, *CC2, *CC3, *eerror;
    
	uuhs = phgAlloc(TotalN * sizeof(*uhs));
	uu_p = phgDofNew(g, DOF_P1, 1, "uu_p", DofInterpolation);
	ff_h = phgDofNew(g, DOF_P1, 1, "ff_h", DofInterpolation);
	phgDofSetDataByFunction(uu_p, func_u0);
	uu_h = phgDofNew(g, DOF_P1, 1, "uu_h", DofInterpolation);
	eerror = phgDofNew(g, DOF_P1, 1, "eerror", DofInterpolation);	
	uuhs[0] = phgDofCopy(uu_p, NULL, DOF_P1, "uu_p1");
	uuhs[0] ->userfunc = DofInterpolation;
	k = 1;
	BB1=phgDofNew(g, DOF_P1, 1, "BB1", func_B1);
	BB2=phgDofNew(g, DOF_P1, 1, "BB2", func_B2);
	BB3=phgDofNew(g, DOF_P1, 1, "BB3", func_B3);
	CC1=phgDofNew(g, DOF_P1, 1, "CC1", func_C1);
	CC2=phgDofNew(g, DOF_P1, 1, "CC2", func_C2);
	CC3=phgDofNew(g, DOF_P1, 1, "CC3", func_C3);

	mmap = phgMapCreate(uu_h, NULL);
	GG1 = phgMapCreateMat(mmap, mmap);
	GG2 = phgMapCreateMat(mmap, mmap);
	CC = phgMapCreateMat(mmap, mmap);

	bbtmp = phgMapCreateVec(mmap, 1);
	phgVecDisassemble(bbtmp);

	build_rhs_Mat(CC, mmap, uu_h);
	build_stiff_Mat1(D0, GG1, mmap, BB1, BB2, BB3, uu_h);
	build_stiff_Mat2(GG2, mmap, CC1, CC2, CC3, uu_h);
	crtime = 0;
	while(crtime < T - 1e-8)
	{
		/********************************************************************/	
		ptime = crtime;
		crtime += stime;
		phgPrintf("\n/********* start new time layer *************/\n");
		phgPrintf("current time layer: [%lf, %lf]\n", (double)ptime, (double)crtime);
		flag++;
		if (flag > 1)
		{   /* update u_p */
			phgDofFree(&uu_p);
			uu_p = phgDofCopy(uu_h, NULL, DOF_P1, "uu_p");
			uu_p->userfunc = DofInterpolation;
		}
		phgMapDofArraysToVec(mmap, 1, FALSE, &bbtmp, &uu_p, NULL);
		bb = phgMapCreateVec(mmap, 1);
		bbtmp1 = phgMapCreateVec(mmap, 1);
		phgVecDisassemble(bb);
		phgVecDisassemble(bbtmp1);

		phgMatVec(0, 1.0, CC, bbtmp, 0.0, &bb);	
		phgDofSetDataByFunction(ff_h, func_f);

		build_rhs_Vec(bbtmp1, ff_h, mmap, uu_h);
		phgVecAXPBY(1.0, bbtmp1, 1.0, &bb);

		AA = phgMapCreateMat(mmap, mmap);
		phgMatAXPBY(1.0, GG1, 0.0, &AA);
		phgMatAXPBY(eta_t(crtime), GG2, 1.0, &AA);

		ssolver = phgSolverCreate(SOLVER_GMRES, uu_h, NULL);
		ssolver->mat = AA;
		ssolver->rhs = bb;
		phgSolverSolve(ssolver, TRUE, uu_h, NULL);
		phgSolverDestroy(&ssolver);

		uuhs[k] = phgDofCopy(uu_h, NULL, DOF_P1, "uu_p1");
		uuhs[k] ->userfunc = DofInterpolation;
		k++;
		phgPrintf("\ntime step: %lg\n\n", (double)stime);
	}

	for(i=1; i< TotalN; i++)
	{    
		phgDofAXPBY(1.0, uhs[i], 0, &eerror);
		phgDofAXPBY(-1.0, uuhs[i], 1.0, &eerror);
		err = phgDofNormL2(eerror);//phgDofNormL2(uhs[i]);
		if(rank==0)
	    	fprintf(snD1, "%f\t%f\n", i*stime, err);
	}
	for(i=0 ;i< TotalN; i++)
		phgDofFree(uuhs+i);
	phgFree(uuhs);
	phgDofFree(&eerror);
	phgDofFree(&BB1);
	phgDofFree(&BB2);
	phgDofFree(&BB3);
	phgDofFree(&ff_h);

	phgDofFree(&CC1);
	phgDofFree(&CC2);
	phgDofFree(&CC3);
	phgDofFree(&uu_h);
	phgDofFree(&uu_p);
}
