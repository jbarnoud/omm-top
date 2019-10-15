import simtk.openmm as mm
from simtk.unit import Quantity


def create_c6c12(cutoff: Quantity) -> mm.CustomNonbondedForce:
    potential = mm.CustomNonbondedForce(
        'C12_ij/(r6^2) - C6_ij/r6; '
        'r6 = r^6; '
        'C6_ij = sqrt(C61*C62); '
        'C12_ij = sqrt(C121*C122)'
    )
    potential.addPerParticleParameter('C6')
    potential.addPerParticleParameter('C12')
    potential.setCutoffDistance(cutoff)
    potential.setNonbondedMethod(potential.CutoffPeriodic)
    return potential


def create_epsilon_sigma(cutoff: Quantity) -> mm.CustomNonbondedForce:
    potential = mm.CustomNonbondedForce(
        '4 * epsilon_ij * (frac6^2 - frac6); '
        'frac6 = (sigma_ij/r)^6; '
        'epsilon_ij = sqrt(epsilon1*epsilon2); '
        'sigma_ij = 0.5 * (sigma1 + sigma2)'
    )
    potential.addPerParticleParameter('sigma')
    potential.addPerParticleParameter('epsilon')
    potential.setCutoffDistance(cutoff)
    potential.setNonbondedMethod(potential.CutoffPeriodic)
    return potential


