import numpy as np
import parmed as pmd

class molecule(object):
    def __init__(self, top_file, traj_file=None):
        self._top_file = top_file
        self._traj_file = traj_file

    def where(self, atoms):
        if not type(atoms) == list:
            atoms = list(atoms)

        top = pmd.load_file(self._top_file).to_dataframe()['name'].values
        return [np.where((top == atom) == True)[0] for atom in atoms]
