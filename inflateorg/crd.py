import numpy as np

class MemProt():
    def __init__(self, u, reference='protein', mobile='resname POPC',
                 headgroup='name P*', threshold=5, sep=None):
        self.u = u
        self.reference = reference
        self.mobile = mobile
        self.threshold = threshold
        self.headgroup = headgroup

        self.mobile_list = self.separate_molecule(sep)
        self.center, self.protein_radius = self.process_protein()

    def separate_molecule(self, sep=None):
        select = self.u.select_atoms(self.mobile)
        if sep is None:
            # If topology is used as input, segment is used to separate the
            # molecules.
            if len(select.segments) > 1:
                return [segment.atoms for segment in select.segments]
            else:
                # If segment is not present, residues are used instead.
                return [residue.atoms for residue in select.residues]
        else:
            return [select.select_atoms(part) for part in sep]

    def process_protein(self):
        """Find the radius of the protein at the plane of lipid head group.

        Returns
        ----------
        radius : float
            The radius of the circle that describes the protein.
        """
        # Find the atoms that are touching the membrane headgroup
        boundry = self.u.select_atoms('{} and around {} ({} and {})'.format(
            self.reference, self.threshold, self.headgroup, self.mobile))
        positions = boundry.positions
        min_x = min(positions[:, 0])
        max_x = max(positions[:, 0])
        min_y = min(positions[:, 1])
        max_y = max(positions[:, 1])
        min_z = min(positions[:, 2])
        max_z = max(positions[:, 2])

        center = [(min_x + max_x) / 2,
                  (min_y + max_y) / 2,
                  (min_z + max_z) / 2]

        distance = np.sqrt(np.sum((positions[:, :2] - center[:2])**2, axis=1))

        return np.array(center), np.mean(distance)

    def find_closest_lipid_radius(self):
        """Find the minimum distance of lipid to the protein center.

        Returns
        ----------
        radius : float
            The radius of the circle that describes the protein.
        """
        lipid = self.u.select_atoms('({} and {}) and around {} {}'.format(
            self.headgroup, self.mobile, self.threshold, self.reference))
        positions = lipid.positions[:, :2]
        center = self.center[:2]
        distance = np.sqrt(np.sum((positions-center)**2, axis=1))
        return min(distance)

    def translate_lipid(self, distance):
        dimensions = self.u.dimensions
        center_of_geometry = self.center
        self.u.atoms.translate(-center_of_geometry)

        for part in self.mobile_list:
            vector = part.center_of_geometry()
            translate = vector[:2] / np.sqrt(np.sum(vector[:2]**2))*distance
            part.translate([translate[0], translate[1], 0])

        dimensions[0] = dimensions[0] + 2 * distance
        dimensions[1] = dimensions[1] + 2 * distance
        self.u.atoms.translate(dimensions[:3] / 2)

        self.u.dimensions = dimensions

    def scale_lipid(self, scaling_factor):
        dimensions = self.u.dimensions
        center_of_geometry = self.center
        self.u.atoms.translate(-center_of_geometry)

        for part in self.mobile_list:
            vector = part.center_of_geometry()
            part.translate((vector * (scaling_factor - 1)))

        self.u.atoms.translate(center_of_geometry * scaling_factor)
        dimensions[:3] = dimensions[:3] * scaling_factor

        self.u.dimensions = dimensions

    def cal_density(self):
        headgroup = self.u.select_atoms('{} and {}'.format(self.headgroup,
                                                           self.mobile))
        mean = np.mean(headgroup.positions[:, -1])
        upper = headgroup.select_atoms('prop z >= {}'.format(mean))
        lower = headgroup.select_atoms('prop z < {}'.format(mean))

        area = self.u.dimensions[0] * self.u.dimensions[1] - \
               np.pi * self.protein_radius ** 2

        return area / len(upper), area / len(lower)