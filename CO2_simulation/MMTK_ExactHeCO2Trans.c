/* C routines for HeCO2FF.py */

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"
//double EnergyEvaluate(double R, double costheta, double *e);
double EnergyEvaluate(double R, double costheta); 
/* This function does the actual energy (and gradient) calculation.
   Everything else is just bookkeeping. */
static void
ExactHeCO2Trans_evaluator(PyFFEnergyTermObject *self,
		   PyFFEvaluatorObject *eval,
		   energy_spec *input,
		   energy_data *energy)
     /* The four parameters are pointers to structures that are
	defined in MMTK/forcefield.h.
	PyFFEnergyTermObject: All data relevant to this particular
                              energy term.
        PyFFEvaluatorObject:  Data referring to the global energy
                              evaluation process, e.g. parallelization
                              options. Not used here.
        energy_spec:          Input parameters for this routine, i.e.
                              atom positions and parallelization parameters.
        energy_data:          Storage for the results (energy terms,
                              gradients, second derivatives).
     */
{
  vector3 *coordinates = (vector3 *)input->coordinates->data;
  vector3 *g;
  int i,j;
  int atom_index1 = (int)self->param[0];  /* atom index */
  int atom_index2 = (int)self->param[1];  /* atom index */
  int atom_index3 = (int)self->param[2];  /* atom index */


/*
  printf("%i %i %i %i %i %i \n",atom_index1,atom_index2,atom_index3,atom_index4,atom_index5,atom_index6);
*/

  //vector12
  double dx = 1.e-9;
  double dx2 = 2.*dx;

  //rab = distance between Ca and Cb
  //  double rab = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  //rab = distance between Ca and Cb
  //  double r0b = sqrt((x2-midx)*(x2-midx)+(y2-midy)*(y2-midy)+(z2-midz)*(z2-midz));
  //r = distance between He and CO2 midpoint
  //  double r = sqrt((x3-midx)*(x3-midx)+(y3-midy)*(y3-midy)+(z3-midz)*(z3-midz));
  // costheta = vector1*vectory2/r*r12
  
  //  double dot = (x1-x2)*(x3-midx)+(y1-y2)*(y3-midy)+(z1-z2)*(z3-midz);
  //  double costheta = dot/(r*rab);
  /* energy_terms is an array because each routine could compute
     several terms that should logically be kept apart. For example,
     a single routine calculates Lennard-Jones and electrostatic interactions
     in a single iteration over the nonbonded list. The separation of
     terms is only done for the benefit of user code (universe.energyTerms())
     returns the value of each term separately), the total energy is
     always the sum of all terms. Here we have only one energy term. */
  double gr[3];
  //double *e;
  //e = (double *) malloc (3*sizeof(double));
  //double *e = (double *) malloc (3*sizeof(double));


  double eplus, eminus;
  double dcostheta = 1.e-5;
  double dr=1.e-5;
  double dr2=2.*dr;
  double dcostheta2 = dcostheta*2.;
  double r_vec[3];

  double rplus;
  double rminus;
  double cosplus;
  double cosminus;
  double drdx[3];
  double dcdx[3];
  double dcostdx[3];
  double bond[3];
  double zunit[3];
  double rab,r,dot,costheta,r3,bondlength;
  double dvdthe;

  double v;
  
  for (i=0;i<3;i++) {
    r_vec[i]=(coordinates[atom_index1][i]-coordinates[atom_index3][i])*10.;
    bond[i]=(coordinates[atom_index2][i]-coordinates[atom_index3][i])*10.;
  }

  r = 0.;
  bondlength=0.;
  for (j=0;j<3;j++) {
    r+=pow(r_vec[j],2.);
    bondlength+=pow(bond[j],2.);
  }
  r = sqrt(r);
  bondlength=sqrt(bondlength);

  costheta=0.;
  for (i=0;i<3;i++){
    zunit[i]=bond[i]/bondlength;
    costheta+=zunit[i]*r_vec[i]/r;
  }
  
  r3 = r*r*r;
  for (j = 0; j < 3; j ++){
    drdx[j] = r_vec[j]/r; 
  }
  //finite dVdtheta
  cosplus = costheta+dcostheta;
  double e = EnergyEvaluate(r, cosplus);
  eplus = e;
  cosminus = costheta-dcostheta;

  e = EnergyEvaluate(r, cosminus);
  eminus = e;
  dvdthe = (eplus-eminus)/dcostheta2;

  dcdx[0] = zunit[0]/r - r_vec[0]*costheta*r/r3;
  dcdx[1] = zunit[1]/r - r_vec[1]*costheta*r/r3;
  dcdx[2] = zunit[2]/r - r_vec[2]*costheta*r/r3;
  
  /*
  dcdx[0] = 1./r - r_vec[0]*r_vec[0]/r3;
  dcdx[1] = -r_vec[0]*r_vec[1]/r3;
  dcdx[2] = -r_vec[0]*r_vec[2]/r3;
  */
  //finit dcosthetadx 
  /*  for (i = 0; i < 3; i ++){
    r_vec[i] += dx;
    r = 0.;
    for (j=0;j<3;j++) {
      r+=pow(r_vec[j]-0.,2.);
    }
    r = sqrt(r); 
    cosplus = r_vec[0]/r;
    r = 0.;
    for (j=0;j<3;j++) {
      r+=pow(r_vec[j]-0.,2.);
    }
    r = sqrt(r); 
    cosminus = r_vec[0]/r;
    dcdx[i] = (cosplus-cosminus)/dx2;
  }
  */
  //  dcdx[3] = -1./r;
  //  double ee[3];
  //  costheta=0.;
  // EnergyEvaluate(r, costheta,e);
  e = EnergyEvaluate(r, costheta);

  // unit convertion to MMTK units
  energy->energy_terms[self->index] = (e*1.196265646e-2);

  /* If only the energy is asked for, stop here. */
  if (energy->gradients == NULL)
    return;

  //  sscanf(line,"%lf ",&e);

  double dvdr;
  rplus = r+dr;
  e = EnergyEvaluate(rplus, costheta);
  eplus = e;
  rminus = r-dr;
  e = EnergyEvaluate(rminus, costheta);
  eminus = e;
  dvdr=(eplus-eminus)/dr2;
  //dvdthe = e[2];

  v = 0;
  for (i=0; i < 3; i ++){
    gr[i]=1.196265646e-1*(drdx[i]*dvdr + dvdthe*dcdx[i]);
    v += gr[i]*r_vec[i]/10.;
  }
  

  energy->energy_terms[self->virial_index] += v;

  /* If only the energy is asked for, stop here. */
  if (energy->gradients == NULL)
    return;
  
  /* Add the gradient contribution to the global gradient array.
     It would be a serious error to use '=' instead of '+=' here,
     in that case all previously calculated forces would be erased. */

  g = (vector3 *)((PyArrayObject*)energy->gradients)->data;

  g[atom_index1][0]+=gr[0];//(double)beads;
  g[atom_index1][1]+=gr[1];//(double)beads;
  g[atom_index1][2]+=gr[2];//(double)beads;
  g[atom_index3][0]-=gr[0];//(double)beads;
  g[atom_index3][1]-=gr[1];//(double)beads;
  g[atom_index3][2]-=gr[2];//(double)beads;
}

/* A utility function that allocates memory for a copy of a string */
static char *
allocstring(char *string)
{
  char *memory = (char *)malloc(strlen(string)+1);
  if (memory != NULL)
    strcpy(memory, string);
  return memory;
}

/* The next function is meant to be called from Python. It creates the
   energy term object at the C level and stores all the parameters in
   there in a form that is convient to access for the C routine above.
   This is the routine that is imported into and called by the Python
   module, HeCO2FF.py. */
static PyObject *
ExactHeCO2TransTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self;
  int atom_index1; 
  int atom_index2;
  int atom_index3; 
  //int beads;

  /* Create a new energy term object and return if the creation fails. */
  self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  /* Convert the parameters to C data types. */
  if (!PyArg_ParseTuple(args, "O!iii",
			&PyUniverseSpec_Type, &self->universe_spec,
			&atom_index1,&atom_index2,&atom_index3))
    return NULL;
  /* We keep a reference to the universe_spec in the newly created
     energy term object, so we have to increase the reference count. */
  Py_INCREF(self->universe_spec);
  /* A pointer to the evaluation routine. */
  self->eval_func = ExactHeCO2Trans_evaluator;
  /* The name of the energy term object. */
  self->evaluator_name = "ExactHeCO2Trans";
  /* The names of the individual energy terms - just one here. */
  self->term_names[0] = allocstring("ExactHeCO2Trans");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  /* self->param is a storage area for parameters. Note that there
     are only 40 slots (double) there, if you need more space, you can use
     self->data, an array for up to 40 Python object pointers. */
  self->param[0] = (double) atom_index1;
  self->param[1] = (double) atom_index2;
  self->param[2] = (double) atom_index3;
  //self->param[1] = (double) beads;
  /* Return the energy term object. */
  return (PyObject *)self;
}

/* This is a list of all Python-callable functions defined in this
   module. Each list entry consists of the name of the function object
   in the module, the C routine that implements it, and a "1" signalling
   new-style parameter passing conventions (only veterans care about the
   alternatives). The list is terminated by a NULL entry. */
static PyMethodDef functions[] = {
  {"ExactHeCO2TransTerm", ExactHeCO2TransTerm, 1},
  {NULL, NULL}		/* sentinel */
};


/* The initialization function for the module. This is the only function
   that must be publicly visible, everything else should be declared
   static to prevent name clashes with other modules. The name of this
   function must be "init" followed by the module name. */
DL_EXPORT(void)
initMMTK_ExactHeCO2Trans(void)
{
  PyObject *m;

  /* Create the module and add the functions. */
  m = Py_InitModule("MMTK_ExactHeCO2Trans", functions);

  /* Import the array module. */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules. */
  import_MMTK_universe();
  import_MMTK_forcefield();

  /* Check for errors. */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_heco2");
}

double EnergyEvaluate(double R, double costheta){
  double energy;
  double LegendrePoly[9];
  double dLegendrePoly[9];
  double De[9];
  double Re[7];
  double betaco[8];
  double PD[4]; //partial derivatives
  double CmVAL[2];
  double dCmVAL[2];
  double C6bar_0 = 72457.62;
  double C8bar_0 = 4.65*C6bar_0;
  double C6bar_2 = 19136.06;
  double C8bar_2 = 4.85*C6bar_0;
  int i; //iteretion
  double r[9];
  int MMLR[2];

  double p[17];
  double dp[17];
  p[0] = 1.;
  p[1] = costheta;
  dp[0] = 0.;
  dp[1] = 1.;

  for (i = 2; i < 17; i ++){
    p[i] = ((2*i-1)*costheta*p[i-1]-(i-1)*p[i-2])/i;
    //    printf("p[%i] = %lf\n", i, p[i]);
  }
  LegendrePoly[0] = p[0];
  for (i = 1; i < 9; i ++){
    LegendrePoly[i] = p[2*i];
  }
  //  printf("Legendre[2] = %lf\n", LegendrePoly[1]);//-0.125
  for (i = 2; i < 17; i ++){
    dp[i] = ((2*i-1)*(p[i-1]+costheta*dp[i-1])-(i-1)*dp[i-2])/i;
    //    printf("dp[%i] = %lf\n", i, dp[i]);
  }
  dLegendrePoly[0] = dp[0];
  for (i = 1; i < 9; i ++){
    dLegendrePoly[i] = dp[2*i];
    //    printf("dLegendrePoly[%i] = %lf\n", i, dLegendrePoly[i]);
  }

  r[0] = R;
  for (i = 1; i < 9; i++){
    r[i] = r[i-1]*R;
  }
 
  De[0] = 32.039; 
  De[1] = -14.787;
  De[2] = 14.616;
  De[3] = -8.049;
  De[4] = 4.315;
  De[5] = -2.135;
  De[6] = 0.975;
  De[7] = -0.439;
  De[8] = 0.184;
  Re[0] = 3.619;
  Re[1] = 0.83644;
  Re[2] = -0.26198;
  Re[3] = 0.09367;
  Re[4] = -0.03118;
  Re[5] = 0.00830;
  Re[6] = -0.00142;
  betaco[0] = 0.0304;
  betaco[1] = 1.0010;
  betaco[2] = -0.0172;
  betaco[3] = 0.6961;
  betaco[4] = 0.180;
  betaco[5] = -0.235;
  betaco[6] = -0.203;
  betaco[7] = 0.15;

  CmVAL[0] = C6bar_0*LegendrePoly[0]+C6bar_2*LegendrePoly[1];
  CmVAL[1] = C8bar_0*LegendrePoly[0]+C8bar_2*LegendrePoly[1];
  
  dCmVAL[0] = C6bar_2 * dLegendrePoly[1];
  dCmVAL[1] = C8bar_2 * dLegendrePoly[1];

  MMLR[0] = 6;
  MMLR[1] = 8;

  double RE = 0.;
  double dREdtheta = 0.;
  for (i = 0; i < 7; i ++){
    RE += Re[i]*LegendrePoly[i];
    dREdtheta += Re[i]*dLegendrePoly[i];
  }
  double uLR = 0.;
  double duLRdr = 0.;
  double duLRdtheta = 0.;
  double uLRE = 0.;
  double duLRE = 0.;
  double dCM = 0.;
  double AA = 0.;
  double BB = 0.;
  for (i = 0; i < 2; i ++){
    AA = CmVAL[i]/r[MMLR[i]-1];
    uLR += AA;
    duLRdr = duLRdr - MMLR[i]*AA/R;
    //    duLRdtheta += dCmVAL[i]/r[MMLR[i]-1];
    BB = CmVAL[i]/pow(RE,MMLR[i]);
    uLRE += BB;
    duLRE = duLRE - MMLR[i]*BB/RE;
  }

    //  duLRE = duLRE*dREdtheta + dCM;
  //dU/dCM
  double U = uLR/uLRE;
  double dUdCM[2];
  for(i = 0 ;i < 2; i ++){
    dUdCM[i] = (1/r[MMLR[i]-1]-U/pow(RE,MMLR[i]))/uLRE;
  }
  double y_3; //when p = 3 since p > (8-6)
  y_3 = (r[2]-pow(RE,3))/(r[2]+pow(RE,3));
  /*  double dy_3dr;
  dy_3dr = 3. * r[1] / (r[2] + pow(RE, 3))
    - 3. * (r[2] - pow(RE, 3))
    * pow(r[2] + pow(RE,3), 2) * r[1];
  double dy_3dtheta;
  dy_3dtheta = -dREdtheta/(r[2]+pow(RE,3))+(r[2]-pow(RE,3))
    *(-(r[2]+pow(RE,3))*(r[2]+pow(RE,3))*dREdtheta);
    */
  double DE = 0.;
  double dDEdr = 0.;
  double dDEdtheta = 0.;
  for (i = 0; i < 9; i ++){
    DE += De[i]*LegendrePoly[i];
    dDEdtheta += De[i]*dLegendrePoly[i];
  }
  double DE_2 = DE*2.;
  
  double beta[4];
  double db[4];
  beta[0] = 0.;
  db[0] = 0.;
  //WE CHANGED THIS TO i<3 TO REMOVE THE LAST BETACO FROM THE TOTAL
  //THEN CORRECTED THE BETA[1] COEFFICIENT TO INCLUDE THIS TERM 
  //THE DERIVATIVE WAS ALSO CORRECTED
  for (i = 0; i < 3; i ++){
    beta[0] += betaco[i]*LegendrePoly[i];
    db[0] += betaco[i]*dLegendrePoly[i];
  }
  beta[1]= betaco[3]*LegendrePoly[0]+betaco[4]*LegendrePoly[1];
  //db[1] = 0.;
  db[1] = betaco[4]*dLegendrePoly[1];
  beta[2] = betaco[5]*LegendrePoly[0]+betaco[6]*LegendrePoly[1];
  db[2] = betaco[6]*dLegendrePoly[1];
  beta[3] = betaco[7]*LegendrePoly[0];
  db[3] = 0.;
  double betainf;
  betainf = log(2.*DE/uLRE);
  double dbetainf = dDEdtheta-duLRE/uLRE;

  double sumbeta = 0.;
  double dsumdr = 0.;
  double dsumdtheta = 0.;
  double yPOW = 1. - y_3;
  double sum;
  double dsum = 0.;
  sum = beta[0]*yPOW;
  for (i = 1; i < 4; i++){
    dsum = dsum + beta[i]*i*yPOW;
    yPOW = yPOW*y_3;
    sum += yPOW*beta[i];
    //    sumbeta += beta[i]*pow(y_3,i);
    //    dsumdr += beta[i]*i*pow(y_3,i-1)*dy_3dr;
    //    dsumdtheta += db[i]*pow(y_3,i)+beta[i]*(i*pow(y_3,i-1)*dy_3dtheta);
  }
  double xp = sum + betainf*y_3;
  //printf("sumbeta=%lf\n",xp);
  double xpw = exp(-xp*y_3)*uLR/uLRE;
  //  printf("E2=%lf\n",xpw);
  double yc = DE*(1.-xpw)*(1.-xpw)-DE; //potential energy
  /*
  double der = 2.*DE*(1.-xpw)*xpw;
  yPOW = der*y_3*(1.-y_3);
  //dV/dbetai
  for (i = 0; i < 4; i ++){
    PD[i] = yPOW;
    yPOW = yPOW*y_3;
  }
  //dVdbetainf
  //  double dVdbinf=der*(-y_3)*y_3;
  //dV/De
  double dVdDE;
  dVdDE = (1.-xpw)*(1.-xpw)+(1.-xpw)*y_3*y_3*uLRE-1.;
  //dype/dre
  double dype = -0.5*(3./RE)*(1.-y_3*y_3);
  dsum = betainf - sum/(1.-y_3) + dsum;
  //dV/dre
  double dVdre = der*(dype * (xp+y_3*dsum)
		      +(1./uLRE+y_3*y_3*uLRE/DE_2)*duLRE);
  
  double dVdr = der*(-dype * (xp+y_3*dsum)
		      -(1./uLR)*duLRdr);
  
  //printf("dVdr = %lf\n", dVdr);
  double dVdtheta1 = 0.;
  for (i = 0; i < 4; i++){
    dVdtheta1 += PD[i]*db[i];
  }
  
  double dVdtheta2 = 0.;
  for (i = 0; i < 2; i++){
    dVdtheta2 += der*dUdCM[i]*dCmVAL[i];
  }
 
  double dVdtheta3 = dVdDE*dDEdtheta;

  double dVdtheta4 = dVdre*dREdtheta;

  // double dVdtheta5 = dVdbinf*dbetainf;
  double dVdtheta = dVdtheta1+dVdtheta2+dVdtheta3+dVdtheta4;
  /* 
  printf("theta1=%lf\n", dVdtheta1);
  printf("theta2=%lf\n", dVdtheta2);
  printf("theta3=%lf\n", dDEdtheta);
  printf("theta4=%lf\n", dVdtheta4);

  printf("theta=%lf\n", dVdtheta);
  printf("CmVAL[0] = %lf\n",CmVAL[0]);
  printf("CmVAL[1] = %lf\n",CmVAL[1]);
  printf("Re = %lf\n", RE);
  printf("De = %lf\n", DE);
  printf("der = %lf\n", der);
  printf("uLR = %lf\n", uLR);
  printf("uLRE = %lf\n", uLRE);* / */
    
  energy= yc;
  return energy;
}
