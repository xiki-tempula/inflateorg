"""
inflateorg.py
Inflate or shrink the membrane to resolve clash between membrane and protein.

Handles the primary functions
"""

import os
import shutil
from subprocess import call
from pkg_resources import resource_filename

import MDAnalysis as mda
import numpy as np

os.environ["GMX_MAXBACKUP"] = "-1"
gromacs = '/usr/local/gromacs/2018.8/bin/gmx'
mdp = resource_filename(__name__, 'data/minim.mdp')
grompp = '{gromacs} grompp -f minim.mdp -c {gro} -p {topol} -o em.tpr -maxwarn {maxwarn} -r {gro}'
mdrun = '{gromacs} mdrun -deffnm em'
trjconv_pbc = 'echo 0 | {gromacs} trjconv -f em.gro -s em.tpr -o em.gro -pbc mol'
trjconv_check = 'echo 0 | {gromacs} trjconv -f pre_inflation.gro -s em.tpr -o em.gro -pbc mol'

class InflateORG():
    def __init__(self, start_file='start.gro', topol='topol.top', center='protein', mobile='not protein', sep=None,
                 scaling_factor = 0.95, dim = [1,1,0], cutoff=1, maxwarn=0):
        '''

        :param start_file: The coordinate file for the InflateAny program.
        :param topol: The topology file compatible with the gromacs program
        :param center: The center of the inflation which is not modified
        :param mobile: The peripheral which will be expanded and shirked
        :param sep: Define how to separate the peripheral.
        :param scaling_factor: The factor of inflation at each
        :param dim: The dimension of the scaling on x, y and z axis. default is 1,1,1.
        :param cutoff: Cutoff distance where two particles are considered as separate.
        '''
        self.start_file = start_file
        self.topol = topol
        self.center = center
        self.mobile = mobile
        self.sep = sep
        self.scaling_factor = scaling_factor
        self.dim = np.array(dim)
        self.cutoff =cutoff
        self.maxwarn = maxwarn

        self.sanity_check()
        self.inflate_system()
        self.shrink_system()

    def sanity_check(self):
        '''
        Check if the input is correct.
        '''
        u = mda.Universe(self.start_file)

        u.select_atoms('({}) or ({})'.format(self.center, self.mobile)).write('pre_inflation.gro')
        try:
            call(grompp.format(gromacs=gromacs, gro='pre_inflation.gro', topol=self.topol, maxwarn=self.maxwarn), shell=True)
            call(trjconv_check.format(gromacs=gromacs), shell=True)
            shutil.move('em.gro', 'pre_inflation.gro')
        except:
            print('Make sure the mdp file (minim.mdp), the topology file ({}) and the input coordinate file ({}) is correct.'.format(
                self.topol,  self.start_file
            ))

        os.mkdir('InflateAny')
        os.remove('em.tpr')
        shutil.move('pre_inflation.gro', 'InflateAny/pre_inflation.gro')
        shutil.copy('minim.mdp', 'InflateAny/minim.mdp')
        os.chdir('InflateAny')
        self.topol = '../' + self.topol


    def separate_molecule(self, selection=None, sep=None):
        '''
        Separate the selection into defferent entities.
        :param u: input universe
        :return: A list of the atom groups which will be scaled.
        '''

        if self.sep is None:
            # By default the separtion is based on residue id
            mobile_parts = u.select_atoms(self.mobile).residues
            return [residue.atoms for residue in mobile_parts]
        else:
            return [u.select_atoms(part) for part in self.sep]

    def inflate(self, u, scaling_factor):
        '''

        :param u: the input MDAnalysis Universe to which scaling will be applied
        '''
        # dimensions = u.dimensions
        # center_of_geometry = u.select_atoms(self.center).center_of_geometry()
        # u.atoms.translate(-center_of_geometry)
        #
        # for part in self.separate_mobile(u):
        #     vector = part.center_of_geometry()
        #     part.translate((vector * (scaling_factor - 1))*self.dim)
        # u.atoms.translate(center_of_geometry * scaling_factor)
        # dimensions[:3] = dimensions[:3] * scaling_factor
        # u.dimensions = dimensions
        # return u

    def inflate_system(self):
        '''
        Inflate the system.
        :return:
        '''
        u = mda.Universe('pre_inflation.gro')
        repeat = True
        count = 0
        print('Start inflating the {}'.format(self.mobile))
        while repeat:
            count += 1
            u = self.inflate(u, 1 / self.scaling_factor)
            check = u.select_atoms('{} and around {} ({})'.format(self.mobile, self.cutoff, self.center))
            print('Interation {}:'.format(count))
            print('Atoms with {}A of {}:'.format(self.cutoff, self.center))
            print(check)
            if len(check) == 0:
                'No atom in contact with {}.'.format(self.center)
                repeat = False
        print('Begin the shrinking process.')
        u.atoms.write('inflated.gro')
        call(grompp.format(gromacs=gromacs, gro='inflated.gro', topol=self.topol, maxwarn=self.maxwarn), shell=True)
        self.mdrun(gromacs=gromacs)
        call(trjconv_pbc.format(gromacs=gromacs), shell=True)
        shutil.copy('em.gro', 'inflated_em.gro')
        self.count = count

    def shrink_system(self):
        for i in range(self.count):
            print('Interation {}:'.format(i))
            u = mda.Universe('em.gro')
            u = self.inflate(u, self.scaling_factor)
            u.atoms.write('shrinked_{}.gro'.format(i))
            call(grompp.format(gromacs=gromacs, gro='shrinked_{}.gro'.format(i), topol=self.topol, maxwarn=self.maxwarn), shell=True)
            self.mdrun(gromacs=gromacs)
            call(trjconv_pbc.format(gromacs=gromacs), shell=True)
        shutil.copy2('em.gro', 'equilibrated.gro')
        shutil.copy2('em.gro', '../equilibrated.gro')
        os.chdir('../')

    def mdrun(self, gromacs=gromacs, additional=''):
        # Try to get around the situation where opencl won't start
        repeat = True
        while repeat:
            returncode = call(mdrun.format(gromacs=gromacs) + ' ' + additional, shell=True)
            if returncode == 0:
                repeat = False
