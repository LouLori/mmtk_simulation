/* C routines for HeHeFF.py */

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"
//double Aziz(double r);
double betafit(double r);
double dVdr(double r);
double num_dVdr(double r);
/* This function does the actual energy (and gradient) calculation.
   Everything else is just bookkeeping. */
static void
hehe_evaluator(PyFFEnergyTermObject *self,
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
  int i;
  int atom_index1 = (int)self->param[0];  /* atom index */
  int atom_index2 = (int)self->param[1];  /* atom index */
  //int beads=(int)self->param[2]; /*number of beads */

/*
  printf("%i %i %i %i %i %i \n",atom_index1,atom_index2,atom_index3,atom_index4,atom_index5,atom_index6);
*/

  double x1 = coordinates[atom_index1][0];
  double y1 = coordinates[atom_index1][1];
  double z1 = coordinates[atom_index1][2];
  
  double x2 = coordinates[atom_index2][0];
  double y2 = coordinates[atom_index2][1];
  double z2 = coordinates[atom_index2][2];

  /* energy_terms is an array because each routine could compute
     several terms that should logically be kept apart. For example,
     a single routine calculates Lennard-Jones and electrostatic interactions
     in a single iteration over the nonbonded list. The separation of
     terms is only done for the benefit of user code (universe.energyTerms())
     returns the value of each term separately), the total energy is
     always the sum of all terms. Here we have only one energy term. */
 
  // unit convertion to MMTK units
  //double kB=3.1668288610848352283e-6; // Hartree/Kelvin
  //double eps = 10.2*kB; //units in Hartree
  //double sigma = 0.228*10.; //units in angstrum
  //double sigma12 = pow(sigma, 12);
  //double sigma6 = pow(sigma,6);
  double e;
  double gr[3];

  double r=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));

  //  e=Aziz(r*10.);
  e = betafit(r*10.);
  //  double dist1 = pow((sigma/r),12);
  //  double dist2 = pow((sigma/r),6);
  //  e = 4*eps*(dist1-dist2);
  //  double dr=1.e-4;
  // double dr2=2.*dr;
  //  double dvdr=(Aziz((r+dr)*10.)-Aziz((r-dr)*10.))*2625.42025277/dr2;
  // double dvdr=(betafit((r+dr))-betafit((r-dr)))/dr2; 
  //   double dvdr = (4*eps*((-12)*sigma12*pow(r,-13) - (-6)*sigma6*pow(r,-7)))*2625.42025277;
   //double dvdr = dVdr(r*10.);
   double dvdr = num_dVdr(r*10.);
  // components
   gr[0]=1.196265646e-1*dvdr*(x1-x2)/r;
   gr[1]=1.196265646e-1*dvdr*(y1-y2)/r;
   gr[2]=1.196265646e-1*dvdr*(z1-z2)/r;
   //gr[3]=1.196265646e-1*dvdr*(x2-x1)/r;
   //gr[4]=1.196265646e-1*dvdr*(y2-y1)/r;
   //gr[5]=1.196265646e-1*dvdr*(z2-z1)/r;
   
  energy->energy_terms[self->index] = (e*1.196265646e-2);
  //energy->energy_terms[self->index] = (e*1.196265646e-2/(double)beads);

  /* If only the energy is asked for, stop here. */
  if (energy->gradients == NULL)
    return;
  
  /* Add the gradient contribution to the global gradient array.
     It would be a serious error to use '=' instead of '+=' here,
     in that case all previously calculated forces would be erased. */

  g = (vector3 *)((PyArrayObject*)energy->gradients)->data;

  // g[atom_index1][0]+=gr[0]/(double)beads;
  

  g[atom_index1][0]+=gr[0];//(double)beads;
  g[atom_index1][1]+=gr[1];//(double)beads;
  g[atom_index1][2]+=gr[2];//(double)beads;
  g[atom_index2][0]-=gr[0];//(double)beads;
  g[atom_index2][1]-=gr[1];//(double)beads;
  g[atom_index2][2]-=gr[2];//(double)beads;


  /*
  g[atom_index1][0]+=gr[0];
  g[atom_index1][1]+=gr[1];
  g[atom_index1][2]+=gr[2];

  g[atom_index2][0]+=gr[3];
  g[atom_index2][1]+=gr[4];
  g[atom_index2][2]+=gr[5];
  */

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
   module, HeHeFF.py. */
static PyObject *
HeHeTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self;
  int atom_index1;
  int atom_index2;
  //int beads;

  /* Create a new energy term object and return if the creation fails. */
  self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  /* Convert the parameters to C data types. */
  if (!PyArg_ParseTuple(args, "O!ii",
			&PyUniverseSpec_Type, &self->universe_spec,
			&atom_index1, &atom_index2))
    //			&atom_index1, &atom_index2,&beads))
    return NULL;
  /* We keep a reference to the universe_spec in the newly created
     energy term object, so we have to increase the reference count. */
  Py_INCREF(self->universe_spec);
  /* A pointer to the evaluation routine. */
  self->eval_func = hehe_evaluator;
  /* The name of the energy term object. */
  self->evaluator_name = "hehe";
  /* The names of the individual energy terms - just one here. */
  self->term_names[0] = allocstring("hehe");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  /* self->param is a storage area for parameters. Note that there
     are only 40 slots (double) there, if you need more space, you can use
     self->data, an array for up to 40 Python object pointers. */
  self->param[0] = (double) atom_index1;
  self->param[1] = (double) atom_index2;
  //self->param[2] = (double) beads;

  /* Return the energy term object. */
  return (PyObject *)self;
}

/* This is a list of all Python-callable functions defined in this
   module. Each list entry consists of the name of the function object
   in the module, the C routine that implements it, and a "1" signalling
   new-style parameter passing conventions (only veterans care about the
   alternatives). The list is terminated by a NULL entry. */
static PyMethodDef functions[] = {
  {"HeHeTerm", HeHeTerm, 1},
  {NULL, NULL}		/* sentinel */
};


/* The initialization function for the module. This is the only function
   that must be publicly visible, everything else should be declared
   static to prevent name clashes with other modules. The name of this
   function must be "init" followed by the module name. */
DL_EXPORT(void)
initMMTK_hehe(void)
{
  PyObject *m;

  /* Create the module and add the functions. */
  m = Py_InitModule("MMTK_hehe", functions);

  /* Import the array module. */
#ifdef import_array
  import_array();
#endif

  /* Import MMTK modules. */
  import_MMTK_universe();
  import_MMTK_forcefield();

  /* Check for errors. */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_hehe");
}
  
double betafit(double r){
  double betacoeff[9];
  betacoeff[0]= -2.472186e-1;
  betacoeff[1]=  2.109e-01;   
  betacoeff[2]= -3.444e-01; 
  betacoeff[3]=  3.36e-01; 
  betacoeff[4]=  2.e-02; 
  betacoeff[5]=  1.360443; 
  betacoeff[6]=  7.3e-01; 
  betacoeff[7] = -1.1;
  betacoeff[8]= -2.2063; 

  double betainf = 0.031316;
  double c6 = 7.041084322e3;  //cm-1
  double c8 = 1.904226812e4;  //cm-1
  double c10 = 6.934738071e4;  //cm-1
  double epsilon = 7.61481; //cm-1
  double Re = 2.96830;  //angstrom
  double De = 7.614810;  //cm-1
  int i;
  double sum = 0;
  double y_5; //p=5
  double y_3;  //q = 3
  y_5 = (pow(r,5)-pow(Re,5))/(pow(r,5)+pow(Re,5));
  y_3 = (pow(r,3)-pow(Re,3))/(pow(r,3)+pow(Re,3));
  for (i = 0; i < 9; i ++){
    sum = sum + betacoeff[i]*pow(y_3,i);
  }
  double betaR;
  betaR = betainf*y_5 + (1-y_5)*sum;
  //potential function
  double uLR;
  uLR = c6/pow(r,6)+c8/pow(r,8)+c10/pow(r,10);
  double uLRe;
  uLRe = c6/pow(Re,6)+c8/pow(Re,8)+c10/pow(Re,10);
  double VMLR;
  VMLR = De*pow((1-(uLR/uLRe)*exp(-betaR*y_5)),2)-De;
    //VMLR = De - (2*De*exp(-betainf)/uLRe)*uLR;
  //  VMLR = De - uLR;
  return VMLR;
}

/*
Function to numerically calculate dVdr(double r) by finite difference
to produce the derivitive of the potential energy
*/

double num_dVdr(double r){
  double div;
  double dr=1.e-5;
  double dr2=2.*dr;

  double eplus = betafit(r+dr);
  double eminus = betafit(r-dr);
  div=(eplus-eminus)/dr2;

  //free(eplus);
  //free(eminus);

  return div;
}

/*
Function dVdr(double r) consume different distances of any two He atoms
to produce the derivitive of potential energy
 */
double dVdr(double r){
  double betacoeff[9];
  betacoeff[0]= -2.472186e-1;
  betacoeff[1]=  2.109e-01;   
  betacoeff[2]= -3.444e-01; 
  betacoeff[3]=  3.36e-01; 
  betacoeff[4]=  2.e-02; 
  betacoeff[5]=  1.360443; 
  betacoeff[6]=  0.73; 
  betacoeff[7] = -1.1;
  betacoeff[8]= -2.2063; 

  double betainf = 0.031316;
  double c6 = 7.041084322e3;  //cm-1
  double c8 = 1.904226812e4;  //cm-1
  double c10 = 6.934738071e4;  //cm-1
  double epsilon = 7.61481; //cm-1
  double Req = 2.96830;  //angstrom
  double Req2 =Req*Req;
  double Req3 =Req2*Req;
  double Req4 =Req2*Req2;
  double Req5 =Req2*Req3;
  double Req6 =Req2*Req4;
  double Req8= Req4*Req4;
  double Req10 =Req5*Req5;

  double r2=r*r;
  double r3=r2*r;
  double r4=r2*r2;
  double r5=r3*r2;
  double r6=r3*r3;
  double r7=r3*r4;
  double r8=r4*r4;
  double r9=r4*r5;
  double r10=r5*r5;
  double r11=r5*r6;

  double De = 7.614810;  //cm-1
  int i;
  //term seperation
  /*the potential equation is Morse/Long-Range (MLR)
    Vmlr(r) = De{1-(uLR/uLRe)*exp(-betaR*yp(r))}^2
    yp(r) = y_5 = (r^5-re^5)/(r^5+re^5)
   */ 
  double sum = 0;
  double y_5; //p=5
  double y_3;  //q = 3
  double betaR;
  double U;
  double W;
  double uLR;
  double uLRe;
  double E;
  double dsum = 0;
  double dy_5;
  double dy_3;
  double dbetaR;
  double duLR;
  double dU;
  double dW;
  double dE;

  y_5 = (r5-Req5)/(r5+Req5);
  y_3 = (r3-Req3)/(r3+Req3);
  for (i = 0; i < 9; i ++){
    sum = sum + betacoeff[i]*pow(y_3,i);
  }
  betaR = betainf*y_5 + (1-y_5)*sum; 
  uLR = c6/r6+c8/r8+c10/r10;
  uLRe = c6/Req6+c8/Req8+c10/Req10;
  U = uLR/uLRe;
  E = exp(-betaR*y_5);
  W = U*E;
  dy_5 = 5. * r4 / (r4 + Req4)
    - 5. * (r5 - Req5) * pow(r5 + Req5, -2) * r4;
  dy_3 = 3. * r * r / (r3 + Req3)
    - 3. * (r3 - Req3) * pow(r3 + Req3, -2) * r2;
  for (i = 0;i < 9; i ++){
    dsum = dsum + betacoeff[i]*i*pow(y_3,i-1)*dy_3;
  }
  
  dbetaR = betainf*dy_5 + (-dy_5)*sum + (1-y_5)*dsum;
  //potential function
 
  duLR = (-6.*c6/r7)+(-8.*c8/r9)+(-10.*c10/r11);
  dE = -E*(dbetaR*y_5+betaR*dy_5);
  dU = duLR/uLRe;
  dW = dU*E+U*dE;
  double div;
  div = 2.*De*(1.-W)*(-dW);
  return div;
}

