import copy
from pathlib import Path
import collections
import itertools
import simtk.openmm as mm
import simtk.openmm.app as app
from simtk.unit import Quantity

import potentials

COMMENT_CHAR = ';'


NB_FUNCTION_LENNARD_JONES = '1'
NB_FUNCTION_BUCKINGHAM = '2'

AtomType = collections.namedtuple('AtomType', ('mass', 'charge', 'ptype', 'v', 'w'))
NonBondLJ = collections.namedtuple('NonBondLJ', ('v', 'w'))
NonBondBuckingham = collections.namedtuple('NonBondBuckingham', ('a', 'b', 'c6'))

Atom = collections.namedtuple('Atom', (
    'index', 'type_', 'residue_number', 'residue_name', 'atom_name',
    'charge_group', 'charge', 'mass',
))
Interaction = collections.namedtuple('Interaction', ('atoms', 'parameters'))


def _interaction_parser(number_of_atoms, interaction_name):
    def wrapped(self, line):
        splitted = line.split()
        atoms = tuple(int(field) for field in splitted[:number_of_atoms])
        parameters = tuple(splitted[number_of_atoms:])
        interaction = Interaction(atoms=atoms, parameters=parameters)
        self.current_molecule.add_interaction(interaction_name, interaction)
    return wrapped


class SectionParser(type):
    """
    Metaclass (!) that populates the `METH_DICT` attribute of new classes. The
    contents of `METH_DICT` are set by reading the `_section_names` attribute
    of all its attributes. You can conveniently set `_section_names` attributes
    using the :meth:`section_parser` decorator.
    """
    def __new__(mcs, name, bases, attrs, **kwargs):
        obj = super().__new__(mcs, name, bases, attrs, **kwargs)
        if not hasattr(obj, 'METH_DICT'):
            obj.METH_DICT = {}

        for attribute_name in dir(obj):
            attribute = getattr(obj, attribute_name)
            try:
                section_name = attribute._section_name
            except AttributeError:
                pass
            else:
                obj.METH_DICT[section_name] = attribute
        return obj

    @staticmethod
    def section_parser(name):
        """
        Parameters
        ----------
        name: str
        """
        def wrapper(method):
            method._section_name = name
            return method
        return wrapper


class TopParser(metaclass=SectionParser):
    bonds = SectionParser.section_parser('bonds')(_interaction_parser(2, 'bonds'))
    pairs = SectionParser.section_parser('pairs')(_interaction_parser(2, 'pairs'))
    angles = SectionParser.section_parser('angles')(_interaction_parser(3, 'angles'))
    dihedrals = SectionParser.section_parser('dihedrals')(_interaction_parser(4, 'dihedrals'))
    constraints = SectionParser.section_parser('constraints')(_interaction_parser(2, 'constraints'))
    settles = SectionParser.section_parser('settles')(_interaction_parser(1, 'settles'))


    def __init__(self, input_path, include_dir, defines):
        self._reset()
        self.parse(input_path, include_dir, defines)

    def parse(self, input_path, include_dir, defines):
        self._reset()
        unknown_sections = set()
        current_section = None
        line_iterator = TopIterator(input_path, include_dir, defines)
        for line in line_iterator:
            section = get_section_or_none(line)
            if section is not None:
                current_section = section
            elif current_section in self.METH_DICT:
                self.METH_DICT[current_section](self, line)
            else:
                # TODO: replace by an error.
                if current_section not in unknown_sections:
                    unknown_sections.add(current_section)
                    print(f'Unknown section [{current_section}].')

    def _reset(self):
        self.gmx_topology = GmxTopology()
        self.current_molecule = None

    @SectionParser.section_parser('defaults')
    def defaults(self, line):
        if self.gmx_topology.has_defaults_defined:
            raise IOError('Defaults are provided more than once.')
        splitted = line.split()
        if len(splitted) < 2:
            raise IOError(
                'The [defaults] section must define at least the non-bonded '
                'function type and the combination rule; not enough values '
                'are being provided.'
            )
        if len(splitted) > 5:
            raise IOError('Too many values provided for the [defaults] section.')
        if splitted[0] in '12':
            self.gmx_topology.non_bonded_function_type = splitted[0]
        else:
            raise IOError(
                f'The value "{splitted[0]}" in the [defaults] section is not '
                'valid for the non-bonded function type (first value of the '
                'line).'
            )
        if splitted[1] in '123':
            self.gmx_topology.combination_rule = splitted[1]
        else:
            raise IOError(
                f'The value "{splitted[1]}" in the [defaults] section is not '
                'valid for the combination rule (second value of the line).'
            )
        if len(splitted) >= 3:
            generate_pair_conversion = {'yes': True, 'no': False}
            raw_generate_pair = splitted[2].lower()
            try:
                self.gmx_topology.generate_pairs = (
                    generate_pair_conversion[raw_generate_pair]
                )
            except KeyError:
                raise IOError(
                    f'The value "{splitted[2]}" in the [defaults] is not '
                    'valid when defining if 1-4 pairs should be generated '
                    '(third value of the line).'
                )
        if len(splitted) >= 4:
            self.gmx_topology.fudge_lennard_jones = float(splitted[3])
        if len(splitted) == 5:
            self.gmx_topology.fudge_coulomb = float(splitted[4])

    @SectionParser.section_parser('atomtypes')
    def atom_types(self, line):
        splitted = line.split()
        if len(splitted) != 6:
            raise ValueError(
                'A line in the [atomtypes] section should have 6 field; '
                f'{len(splitted)} found instead.'
            )
        name = splitted[0]
        self.gmx_topology.atom_types[name] = AtomType(
            mass=float(splitted[1]),
            charge=float(splitted[2]),
            ptype=splitted[3],
            v=float(splitted[4]),
            w=float(splitted[4]),
        )

    @SectionParser.section_parser('nonbond_params')
    def nonbond_params(self, line):
        splitted = line.split()
        if splitted[2] == NB_FUNCTION_LENNARD_JONES:
            self.nonbond_params_lj(line)
        elif splitted[2] == NB_FUNCTION_BUCKINGHAM:
            self.nonbond_params_buckingham(line)
        else:
            raise IOError(f'Unexpected function {splitted[2]} in [nonbond_params].')

    def nonbond_params_lj(self, line):
        splitted = line.split()
        self.gmx_topology.non_bonded_params[splitted[0]][splitted[1]] = NonBondLJ(
            v=float(splitted[3]), w=float(splitted[4]),
        )

    def nonbond_params_buckingham(self, line):
        splitted = line.split()
        self.gmx_topology.non_bonded_params[splitted[0]][splitted[1]] = NonBondBuckingham(
            a=float(splitted[3]), b=float(splitted[4]), c6=float(splitted[5]),
        )

    @SectionParser.section_parser('moleculetype')
    def moleculetype(self, line):
        name, exclusion_distance_str = line.split()
        exclusion_distance = int(exclusion_distance_str)
        # TODO: Warn if a molecule name is already in use.
        self.current_molecule = self.gmx_topology.add_molecule_type(name, exclusion_distance)

    @SectionParser.section_parser('atoms')
    def atoms(self, line):
        splitted = line.split()
        
        charge = 0
        mass = None
        try:
            charge = float(splitted[6])
        except IndexError:
            pass
        try:
            mass = float(splitted[7])
        except IndexError:
            pass

        atom = Atom(
            index=int(splitted[0]),
            type_=splitted[1],
            residue_number=int(splitted[2]),
            residue_name=splitted[3],
            atom_name=splitted[4],
            charge_group=int(splitted[5]),
            charge=charge,
            mass=mass,
        )
        self.current_molecule.add_atom(atom)

    @SectionParser.section_parser('system')
    def system(self, line):
        self.gmx_topology.name = line.strip()

    @SectionParser.section_parser('molecules')
    def molecules(self, line):
        name, number_str = line.split()
        number = int(number_str)
        self.gmx_topology.molecules.append((name, number))


class MoleculeType:
    def __init__(self, exclusion_distance):
        self.exclusion_distance = exclusion_distance
        self.atoms = {}
        self.interactions = collections.defaultdict(list)

    def add_atom(self, atom):
        self.atoms[atom.index] = atom

    def add_interaction(self, interaction_name, interaction):
        self.interactions[interaction_name].append(interaction)


class GmxTopology:
    def __init__(self):
        self.name = ''
        self.non_bonded_function_type = None
        self.combination_rule = None
        self.generate_pairs = False
        self.fudge_lennard_jones = 1
        self.fudge_coulomb = 1

        self.molecule_types = {}
        self.molecules = []
        self.atom_types = {}
        self.non_bonded_params = collections.defaultdict(dict)

    @property
    def has_defaults_defined(self):
        return self.non_bonded_function_type is not None

    def add_molecule_type(self, name, exclusion_distance):
        molecule_type = MoleculeType(exclusion_distance)
        self.molecule_types[name] = molecule_type
        return molecule_type

    def create_omm_topology(self):
        topology = app.Topology()
        for molecule_name, number_of_molecule_copies in self.molecules:
            for _ in range(number_of_molecule_copies):
                self.add_moltype_to_omm_topology(
                    self.molecule_types[molecule_name], topology,
                )
        return topology

    def setup_system_non_bonded(self, system: mm.System, cutoff: Quantity):
        if self.non_bonded_function_type != NB_FUNCTION_LENNARD_JONES:
            # TODO: Implement Buckingham potential
            raise NotImplementedError(
                f'Non bonded function {self.non_bonded_function_type} '
                'is not implemented.'
            )

        if self.combination_rule in '1':
            potential = potentials.create_c6c12(cutoff=cutoff)
        elif self.combination_rule in '2':
            potential = potentials.create_epsilon_sigma(cutoff=cutoff)
        else:
            # TODO: implement combination rule 3
            raise NotImplementedError(
                f'Combination rule {self.combination_rule} is not implemented.'
            )

        for molecule_name, number_of_molecule_copies in self.molecules:
            molecule = self.molecule_types[molecule_name]
            for _ in range(number_of_molecule_copies):
                for atom in molecule.atoms.values():
                    ptype = self.atom_types[atom.type_]
                    potential.addParticle([ptype.v, ptype.w])

        system.addForce(potential)
        return potential

    def createSystem(self, cutoff: Quantity):
        topology = self.create_omm_topology()
        system = mm.System()
        self.setup_system_non_bonded(system, cutoff=cutoff)
        return system

    @staticmethod
    def add_moltype_to_omm_topology(molecule_type: MoleculeType, topology: app.Topology):
        gmx_to_omm_index = {}
        atom_index = topology.getNumAtoms() - 1
        chain = topology.addChain()
        residue_iterator = itertools.groupby(
            molecule_type.atoms.items(),
            key=lambda x: (x[1].residue_number, x[1].residue_name),
        )
        for (residue_index, residue_name), residue_atoms in residue_iterator:
            residue = topology.addResidue(name=residue_name, chain=chain)
            for atom_index, (_, atom) in enumerate(residue_atoms, start=atom_index + 1):
                topology.addAtom(atom_index, atom.type_, residue)
                gmx_to_omm_index[atom.index] = atom_index
        for interaction in molecule_type.interactions['bonds']:
            topology.addBond(
                gmx_to_omm_index[interaction.atoms[0]],
                gmx_to_omm_index[interaction.atoms[1]],
            )


class TopIterator:
    """
    Iterate over the lines of a TOP/ITP file, and its included files.
    """
    def __init__(self, input_path, include_dir, defines):
        self._file_stack = []
        self._line_stack = []
        self._include_dir = include_dir
        self._defines = copy.copy(defines)
        self._starting_file = input_path

    def _iter_from(self, input_path):
        with open(input_path) as infile:
            self._file_stack.append(input_path)
            self._line_stack.append(0)
            for line in clean_lines(infile):
                self._line_stack[-1] += 1
                if line.startswith('#include'):
                    yield from self._include(line)
                elif line.startswith('#if'):
                    self._do_if_statement(line, infile)
                elif line.startswith('#endif'):
                    pass
                elif line.startswith('#define'):
                    self._do_define_statement(line)
                elif line:
                    yield line
                # "else" would be an empty line, we want to ignore those.
            self._file_stack.pop()
            self._line_stack.pop()

    def __iter__(self):
        yield from self._iter_from(self._starting_file)

    @property
    def current_file(self):
        return self._file_stack[-1]

    @property
    def current_line(self):
        return self._line_stack[-1]

    def _do_if_statement(self, line, infile):
        variable = line.split(maxsplit=1)[-1]
        if line.startswith('#ifdef'):
            if variable not in self._defines:
                self._skip_lines_until_endif(infile)
        elif line.startswith('#ifndef'):
            if variable in self._defines:
                self._skip_lines_until_endif(infile)

    def _skip_lines_until_endif(self, infile):
        for line in clean_lines(infile):
            self._line_stack[-1] += 1
            if line.startswith('#endif'):
                break
        else:  # no break
            raise IOError(f'Missing "#endif" in {self.current_file}.')

    def _include(self, line):
        include_path = line.split(maxsplit=1)[-1][1:-1]
        include_path = find_path(include_path, self.current_file, self._include_dir)
        yield from self._iter_from(include_path)

    def _do_define_statement(self, line):
        _, variable, *value = line.split()
        self._defines[variable] = value


def uncomment(line):
    comment_start = line.find(COMMENT_CHAR)
    if comment_start < 0:
        return line
    return line[:comment_start]


def clean_lines(lines):
    for line in lines:
        yield uncomment(line).strip()


def find_path(path, current_path, include_dir):
    current_path = Path(current_path)
    parent = current_path.parent
    path_to_include = parent / path
    if path_to_include.exists():
        return str(path_to_include)
    path_to_include = Path(include_dir) / path
    if path_to_include.exists():
        return str(path_to_include)
    raise FileNotFoundError(f'Could not find a file matching "{path}".')


def get_section_or_none(line):
    stripped = line.strip()
    if stripped and stripped[0] == '[' and stripped[-1] == ']':
        return stripped[1:-1].strip()
    return None
