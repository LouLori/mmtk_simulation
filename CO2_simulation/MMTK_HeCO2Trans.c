/* C routines for HeCO2FF.py */

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/forcefield_private.h"
#include <math.h>
#include <stdio.h>
//double EnergyEvaluate(double R, double costheta, double *e);
//double EnergyEvaluate(double R, double costheta); 
/* This function does the actual energy (and gradient) calculation.
   Everything else is just bookkeeping. */
static void
HeCO2Trans_evaluator(PyFFEnergyTermObject *self,
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

  PyArrayObject *pot_array = (PyArrayObject *)self->data[0];
  double *pot = (double *)pot_array->data;  /* atomic charges */
  
  PyArrayObject *dvdr_array = (PyArrayObject *)self->data[1];
  double *dvdr = (double *)dvdr_array->data;  /* atomic charges */
  
  PyArrayObject *dvdct_array = (PyArrayObject *)self->data[2];
  double *dvdct = (double *)dvdct_array->data;  /* atomic charges */
  

  //printf("%i %i %i \n",atom_index1,atom_index2,atom_index3);
  
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

  double e;
  double gr[3];
  
  double r_vec[3];
  double bond[3];
  double zunit[3];
  
  double drdx[3];
  double dcdx[3];

  double r,costheta,r3,bondlength;
  double v;

  double rmin=1.5; //Ang
  int nr=23501;
  double dr=0.001; //Ang

  double ctmin=-1.0;
  int nct=601;
  double dct=2.0/600.0;

  int rindex, ctindex, arrayindex;
  int ihighr, ihighct;
  double vadj, dvdradj, dvdctadj;
  
  //Convert lengths to Angstrom
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

  //printf("%f %f %f \n",r,costheta,drdx[0]);
  
  rindex=(int)round((r-rmin)/dr);
  ctindex=(int)round((costheta-ctmin)/dct);
  arrayindex=ctindex+nct*rindex;

  if (rindex>nr-1)
    rindex=nr-1;
  if (rindex<0)
    rindex=0;
  if (ctindex>nct-1)
    ctindex=nct-1;
  if (ctindex<0)
    ctindex=0;
  
  //printf("%i %i %i \n",rindex,ctindex,arrayindex);
  
  ihighr=arrayindex+nct;
  ihighct=arrayindex+1;
  
  vadj=pot[arrayindex];
  dvdradj=dvdr[arrayindex];
  dvdctadj=dvdct[arrayindex];


  if (ihighr < nct*nr && ihighct < (rindex+1)*nct)
    {
      vadj=vadj+((pot[ihighr]-pot[arrayindex])/dr)*
	(r-(rindex*dr+rmin));
      vadj=vadj+((pot[ihighct]-pot[arrayindex])/dct)*(costheta-(ctindex*dct+ctmin));
      
      dvdradj=dvdradj+((dvdr[ihighr]-dvdr[arrayindex])/dr)*
	(r-(rindex*dr+rmin));
      dvdradj=dvdradj+((dvdr[ihighct]-dvdr[arrayindex])/dct)*(costheta-(ctindex*dct+ctmin));
      
      dvdctadj=dvdctadj+((dvdct[ihighr]-dvdct[arrayindex])/dr)*
	(r-(rindex*dr+rmin));
      dvdctadj=dvdctadj+((dvdct[ihighct]-dvdct[arrayindex])/dct)*(costheta-(ctindex*dct+ctmin));
    }

  e=vadj;

  //printf("%f \n",e);
  
  dcdx[0] = zunit[0]/r - r_vec[0]*costheta*r/r3;
  dcdx[1] = zunit[1]/r - r_vec[1]*costheta*r/r3;
  dcdx[2] = zunit[2]/r - r_vec[2]*costheta*r/r3;
  
  // unit convertion to MMTK units
  energy->energy_terms[self->index] = e;

  /* If only the energy is asked for, stop here. */
  if (energy->gradients == NULL)
    return;

  //  sscanf(line,"%lf ",&e);

  v = 0;
  for (i=0; i < 3; i ++){
    gr[i]=(drdx[i]*dvdradj + dvdctadj*dcdx[i]*10.0);  //(1/Ang -> 1/nm) in dc/dx
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
HeCO2TransTerm(PyObject *dummy, PyObject *args)
{
  PyFFEnergyTermObject *self;
  PyArrayObject *pot;
  PyArrayObject *dvdr;
  PyArrayObject *dvdct;
  int atom_index1; 
  int atom_index2;
  int atom_index3; 

  /* Create a new energy term object and return if the creation fails. */
  self = PyFFEnergyTerm_New();
  if (self == NULL)
    return NULL;
  /* Convert the parameters to C data types. */
  if (!PyArg_ParseTuple(args, "O!O!O!O!iii",
			&PyUniverseSpec_Type, &self->universe_spec,
                        &PyArray_Type, &pot,
                        &PyArray_Type, &dvdr,
                        &PyArray_Type, &dvdct,
			&atom_index1,&atom_index2,&atom_index3))
    return NULL;
  /* We keep a reference to the universe_spec in the newly created
     energy term object, so we have to increase the reference count. */
  Py_INCREF(self->universe_spec);
  /* A pointer to the evaluation routine. */
  self->eval_func = HeCO2Trans_evaluator;
  /* The name of the energy term object. */
  self->evaluator_name = "HeCO2Trans";
  /* The names of the individual energy terms - just one here. */
  self->term_names[0] = allocstring("HeCO2Trans");
  if (self->term_names[0] == NULL)
    return PyErr_NoMemory();
  self->nterms = 1;
  /* self->param is a storage area for parameters. Note that there
     are only 40 slots (double) there, if you need more space, you can use
     self->data, an array for up to 40 Python object pointers. */
  self->param[0] = (double) atom_index1;
  self->param[1] = (double) atom_index2;
  self->param[2] = (double) atom_index3;

  self->data[0] = (PyObject *)pot;
  self->data[1] = (PyObject *)dvdr;
  self->data[2] = (PyObject *)dvdct;


  Py_INCREF(pot);
  Py_INCREF(dvdr);
  Py_INCREF(dvdct);
  
  /* Return the energy term object. */
  return (PyObject *)self;
}

/* This is a list of all Python-callable functions defined in this
   module. Each list entry consists of the name of the function object
   in the module, the C routine that implements it, and a "1" signalling
   new-style parameter passing conventions (only veterans care about the
   alternatives). The list is terminated by a NULL entry. */
static PyMethodDef functions[] = {
  {"HeCO2TransTerm", HeCO2TransTerm, 1},
  {NULL, NULL}		/* sentinel */
};


/* The initialization function for the module. This is the only function
   that must be publicly visible, everything else should be declared
   static to prevent name clashes with other modules. The name of this
   function must be "init" followed by the module name. */
DL_EXPORT(void)
initMMTK_HeCO2Trans(void)
{
  PyObject *m;

  /* Create the module and add the functions. */
  m = Py_InitModule("MMTK_HeCO2Trans", functions);

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

