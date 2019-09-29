# Quartic potential with respect to a fixed point in space
from MMTK.ForceFields.ForceField import ForceField
from MMTK_hehe import HeHeTerm

class HeHeForceField(ForceField):

    """HeHe potential with respect to a fixed point in space

    Constructor: HeHeForceField(|atom|, |center|,
                                             |force_constant|)

    Arguments:

    |atom| -- an atom object or an integer atom index, specifying the
              atom on which the force field acts

    """
    #def __init__(self, atom1, atom2, beads):

    def __init__(self, atom1, atom2):
        self.id1, self.id2 = self.getAtomParameterIndices([atom1, atom2])
        self.arguments = (self.id1, self.id2)
        #self.arguments = (self.id1, self.id2, beads)
        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, 'hehe')
        # Store the parameters for later use.
        #self.beads=beads

    # The following method is called by the energy evaluation engine
    # to inquire if this force field term has all the parameters it
    # requires. This is necessary for interdependent force field
    # terms. In our case, we just say "yes" immediately.
    def ready(self, global_data):
	return True

    def supportsPathIntegrals(self):
        return True

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        if subset1 is not None or subset2 is not None:
            raise ValueError, "sorry, no subsets here"
        # Here we pass all the parameters as "simple" data types to
        # the C code that handles energy calculations.
#        return [paraH2Term(universe._spec,
#                                       self.atom1,self.atom2)]

        f, offsets = self.beadOffsetsAndFactor([self.id1, self.id2], global_data)
        return [HeHeTerm(universe._spec,
                            self.id1 + o1,
                            self.id2 + o2) for o1,o2 in offsets]

       # return [HeHeTerm(universe._spec,
       #                  atom1.index, atom2.index,self.beads)]

    # This method returns the string that is inserted into the universe
    # descriptions in trajectories. It is the class name followed by
    # the arguments, just what it takes to re-create an equivalent object.
    def description(self):
        return "HeHeFF.HeHeForceField" + `self.arguments`
