import os
from subprocess import call
from pkg_resources import resource_filename

os.environ["GMX_MAXBACKUP"] = "-1"
class Gromacs():
    def __init__(self, topol, gmx=None):
        if gmx is None:
            self.gmx = gmx
        else:
            pass
        self.topol = topol

        self.grompp = '{gromacs} grompp -f minim.mdp -c {gro} -p {topol} -o ' \
                  'em.tpr -maxwarn {maxwarn} -r {gro}'

    def autoimage(self, inpcrd, outcrd):
        trjconv_pbc = 'echo 0 | {gromacs} trjconv -f {inpcrd} -s em.tpr ' \
                      '-o {outcrd} -pbc mol'
        call(trjconv_pbc.format(gromacs=self.gmx, inpcrd=inpcrd, outcrd=outcrd),
             shell=True)

    def mdrun(self, inpcrd, outcrd):
        mdrun = '{gromacs} mdrun -deffnm em'
        call(mdrun.format(gromacs=self.gmx), shell=True)

    def grompp(self, inpcrd):
        grompp = '{gromacs} grompp -f {minim} -c {gro} -p {topol} -o ' \
                 'em.tpr -maxwarn -1 -r {gro}'
        call(grompp.format(gromacs=self.gmx, gro=inpcrd, topol=self.topol,
                           minim=resource_filename(__name__, 'data/minim.mdp')),
             shell=True)