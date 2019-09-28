# Quartic potential with respect to a fixed point in space
from MMTK.ForceFields.ForceField import ForceField
from MMTK_ExactHeCO2Trans import ExactHeCO2TransTerm

class ExactHeCO2TransForceField(ForceField):

    """HeCO2 potential with respect to an x-axis aligned point CO2 free to translate

    Constructor: HeCO2TransForceField(|co2 atom|, |he atom|, |center|,
                                             |force_constant|)

    Arguments:

    |atom| -- an atom object or an integer atom index, specifying the
              atom on which the force field acts

    """
    def __init__(self, atom1, atom2,atom3):
        self.index1, self.index2, self.index3 = self.getAtomParameterIndices([atom1, atom2, atom3])
        self.arguments = (self.index1, self.index2, self.index3)

        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, 'ExactHeCO2Trans')

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

        f, offsets = self.beadOffsetsAndFactor([self.index1,self.index2,self.index3], global_data)
        return [ExactHeCO2TransTerm(universe._spec, self.index1+o1, self.index2+o2,self.index3+o3) for o1, o2, o3 in offsets]


    # This method returns the string that is inserted into the universe
    # descriptions in trajectories. It is the class name followed by
    # the arguments, just what it takes to re-create an equivalent object.
    def description(self):
        return "ExactHeCO2TransFF.ExactHeCO2TransForceField" + `self.arguments`
