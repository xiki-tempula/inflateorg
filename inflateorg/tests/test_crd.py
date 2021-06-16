import pytest
from pkg_resources import resource_filename

import MDAnalysis as mda
from inflateorg.crd import MemProt

class TestSegement():
    @staticmethod
    @pytest.fixture(scope='class')
    def memprot():
        '''Simulate the case where the lipid has three residues and is
        classified into segements'''
        u = mda.Universe(resource_filename(__name__, 'test_data/segement.gro'))
        core_segment = u.add_Segment(segid='A')
        u.select_atoms('protein').residues.segments = core_segment
        core_segment = u.add_Segment(segid='B')
        u.select_atoms('resname SOL NA CL').residues.segments = core_segment
        for i in range(186):
            core_segment = u.add_Segment(segid='P{}'.format(i))
            u.segments[0].atoms[:118].residues.segments = core_segment

        return MemProt(u, reference='protein', mobile='resname MY PC',
                 headgroup='name P*', threshold=5)

    def test_cal_density(self, memprot):
        assert [59.9, 58.7] == pytest.approx(memprot.cal_density(), 0.1)

    def test_separate_molecule_segment(self, memprot):
        assert len(memprot.mobile_list) == 186

    def test_separate_molecule_sep(self, memprot):
        assert len(memprot.separate_molecule(['bynum 3622-3631',
                                              'bynum 3632-3641',]
                                             )) == 2

    def test_process_protein(self, memprot):
        assert 17.2 == pytest.approx(memprot.protein_radius, 0.1)
        assert [41.2, 38.8, 52.8] == pytest.approx(memprot.center, 0.1)

    def test_find_closest_lipid_radius(self, memprot):
        assert 8.2 == pytest.approx(memprot.find_closest_lipid_radius(), 0.1)

    def test_translate_lipid(self, memprot):
        PO4 = memprot.u.select_atoms('name P*')[0]
        memprot.translate_lipid(
            memprot.protein_radius - \
            memprot.find_closest_lipid_radius())
        assert [24.6, 22.6, 72.0] == pytest.approx(PO4.position, 0.1)

    def test_scale_lipid(self, memprot):
        PO4 = memprot.u.select_atoms('name P*')[0]
        memprot.scale_lipid(1.05)
        assert [24.6, 22.6, 72.0] == pytest.approx(PO4.position, 0.1)