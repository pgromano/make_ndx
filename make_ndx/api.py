import numpy as np
import parmed as pmd

class Molecule(object):
    def __init__(self, topfile):
        '''
        Identifies the file format of the specified topology file and returns
        its parsed contents.

        Parameters
        ----------
        topfile : str
            The name of the topology file to try to parse. If the filename starts with
            http:// or https:// or ftp://, it is treated like a URL and the file will be
            loaded directly from its remote location on the web
        '''
        # Create topology dataframe
        topology = pmd.load_file(topfile, structure=True).to_dataframe()

        # Create chainid key
        topology['chainid'] = topology['chain']

        # Create chainid values
        chains = np.unique(topology['chain'].as_matrix())
        for i,chain in enumerate(chains):
            topology['chainid'] = topology['chainid'].replace(chain, i)

        self._name = topology[['chain', 'resname', 'name']].as_matrix()
        self._index = topology[['chainid', 'resid', 'number']].as_matrix()
        self._xyz = topology[['xx', 'xy', 'xz']].as_matrix()

    def chains(self, chainid=None):
        ''' List of chain types.
        chainid: int (Default=None)
            Index value to get chain type by chain sequence. Default is None which
            returns all chain types.
        '''
        return np.squeeze(np.unique(self._name[:,0])[chainid])

    @property
    def n_chains(self):
        'Number of chains in topology'
        return len(self.chains())

    def residues(self, chainid=None):
        ''' List of residue types.
        chainid: int (Default=None)
            Index value to get residues from specific chain type. Default is
            None which returns all residues from all chain types.
        '''
        res_idx = list(np.concatenate([[0], np.where(self._index[:-1,1] != self._index[1:,1])[0]+1]))
        chain_per_res = self._name[res_idx,0]
        residues = self._name[res_idx,1]

        if chainid == None:
            return np.squeeze(residues)
        else:
            return np.squeeze(residues[chain_per_res == self.chains(chainid)])

    def n_residues(self, chainid=None):
        '''Number of residues in topology
        chainid: int (Default=None)
            Index value to get number of residues within specific chain type.
            Default is None which returns the total number of residues across
            all chains.
        '''
        return len(self.residues(chainid))

    def atoms(self, chainid=None, resid=None):
        ''' List of residue types.
        chainid: int (Default=None)
            Index value to get atoms from specific chain type. Default is
            None which returns all residues from all chain types.
        resid: int (Default=None)
            Index value to get atoms from specific residue atom. Default is None
            which returns all atoms from all residues. Note if chainid is
            supplied, the residue index begins specific to chain type.
        '''
        if chainid is None:
            if resid is None:
                return self._name[:,2]
            else:
                return self._name[self._index[:,1] == resid, 2]
        else:
            if resid is None:
                return self._name[self._index[:,0] == chainid, 2]
            else:
                return self._name[np.logical_and(self._index[:,0] == chainid, self._index[:,1] == resid), 2]

    @property
    def n_atoms(self):
        'Number of atoms in topology'
        return self._name.shape[0]

    def where(self, atoms, chainid=None, resid=None):
        ''' Find index for all supplied atom types
        Parameters:
        -----------
        atoms: str, list
            A string or list of strings of the atom types to be index from
            topology file.
        chainid: int (Default=None)
            Index value to get atoms from specific chain type. Default is
            None which returns all residues from all chain types.
        resid: int (Default=None)
            Index value to get atoms from specific residue atom. Default is None
            which returns all atoms from all residues. Note if chainid is
            supplied, the residue index begins specific to chain type.
        Returns:
        --------
        index: int, list
            Returns index or list of indices where queried atom types are found.
        '''
        if type(atoms) == str:
            atoms = list(atoms)

        if chainid is None:
            if resid is None:
                return np.unique(np.concatenate([np.where(self._name[:,2] == atom)[0] for atom in atoms]))
            else:
                idx = np.where(self._index[:,1] == resid)[0]
                return np.unique(np.concatenate([idx[np.where(self._name[idx,2] == atom)[0]] for atom in atoms]))
        else:
            if resid is None:
                idx = np.where(self._index[:,0] == chainid)[0]
                return np.unique(np.concatenate([idx[np.where(self._name[idx,2] == atom)[0]] for atom in atoms]))
            else:
                if chainid != 0:
                    resid += np.sum([self.n_residues(i) for i in range(chainid)])
                idx = np.where(self._index[:,1] == resid)[0]
                return np.unique(np.concatenate([idx[np.where(self._name[idx,2] == atom)[0]] for atom in atoms]))

    @property
    def positions(self):
        '''Return the atomic cartesian coordinates array'''
        return self._xyz
