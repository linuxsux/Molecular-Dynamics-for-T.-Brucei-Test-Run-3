from openmm.app import PDBFile
from openmm.app import *
from openmm import *
from openmm.app import CharmmPsfFile
from openmm.unit import *
from sys import stdout
from copy import deepcopy
import openmm.app as app
import os

def add_backbone_posres(system, positions, atoms, restraint_force):
  force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
  force_amount = restraint_force * kilocalories_per_mole/angstroms**2
  force.addGlobalParameter("k", force_amount)
  force.addPerParticleParameter("x0")
  force.addPerParticleParameter("y0")
  force.addPerParticleParameter("z0")
  for i, (atom_crd, atom) in enumerate(zip(positions, atoms)):
    if atom.name in  ('CA', 'C', 'N', 'P', 'S', 'O'):
      force.addParticle(i, atom_crd.value_in_unit(nanometers))
  posres_sys = deepcopy(system)
  posres_sys.addForce(force)
  return posres_sys

def add_backbone_posres1(system, positions, atoms, restraint_force):
  force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
  force_amount = restraint_force * kilocalories_per_mole/angstroms**2
  force.addGlobalParameter("k", force_amount)
  force.addPerParticleParameter("x0")
  force.addPerParticleParameter("y0")
  force.addPerParticleParameter("z0")
  for i, (atom_crd, atom) in enumerate(zip(positions, atoms)):
    if atom.name in  ('CA', 'C', 'N'):
      force.addParticle(i, atom_crd.value_in_unit(nanometers))
  posres_sys1 = deepcopy(system)
  posres_sys1.addForce(force)
  return posres_sys1

# Setup
pdb = PDBFile('centered_output.pdb')

# Popular protein force field plus respective water model force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
modeller=Modeller(pdb.topology, pdb.positions)

# Ensuring entire protein is protonated
modeller.addHydrogens(forcefield)
modeller.addExtraParticles(forcefield)

# Adds water explicitly
# Provides .15M of NaCl 
modeller.addSolvent(forcefield, ionicStrength=0.15*molar)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)
posres_sys = add_backbone_posres(system, modeller.positions, modeller.topology.atoms(), 100)
integrator = LangevinMiddleIntegrator(10*kelvin, 1/picosecond, 2*femtoseconds)
integrator.setConstraintTolerance(0.0001)
simulation = Simulation(modeller.topology, system, integrator)
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True, density=True))

# Minimize
# Run for 4 ps
print('Minimizing...')
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()
simulation.step(2000)

# Begin equilibration step 1: HP 50 kcalmol-1k-2, 100ps at 10 K, restraining heavy atoms incl. oxygen
# Lasts 100 ps
print('Step 1 equilibration commencing...')
simulation.context.setPositions(modeller.positions)
posres_sys = add_backbone_posres(system, modeller.positions, modeller.topology.atoms(), 50)
integrator = LangevinMiddleIntegrator(10*kelvin, 1/picosecond, 2*femtoseconds)
integrator.setConstraintTolerance(0.0001)
simulation.step(50000)

# Begin equilibration step 2: HP 50 kcalmol-1k-2, 100ps at 10 K, restraining heavy atoms excl. oxygen
# Lasts 100 ps
print('Step 2 equilibration commencing...')
simulation.context.setPositions(modeller.positions)
posres_sys1 = add_backbone_posres1(system, modeller.positions, modeller.topology.atoms(), 50)
integrator = LangevinMiddleIntegrator(10*kelvin, 1/picosecond, 2*femtoseconds)
integrator.setConstraintTolerance(0.0001)
simulation.step(50000)

# Begin equilibration step 3: HP 5 kcalmol-1k-2, 100ps at 10 K, restraining heavy atoms excl. oxygen
# Lasts 100 ps
print('Step 3 equilibration commencing...')
simulation.context.setPositions(modeller.positions)
posres_sys1 = add_backbone_posres1(system, modeller.positions, modeller.topology.atoms(), 5)
integrator = LangevinMiddleIntegrator(10*kelvin, 1/picosecond, 2*femtoseconds)
integrator.setConstraintTolerance(0.0001)
simulation.step(50000)

# Begin equilibration step 4: HP 0 kcalmol-1k-2, 100ps at 10 K, restraining heavy atoms excl. oxygen
# Lasts 100 ps
print('Step 4 equilibration commencing...')
simulation.context.setPositions(modeller.positions)
posres_sys1 = add_backbone_posres1(system, modeller.positions, modeller.topology.atoms(), 0)
integrator = LangevinMiddleIntegrator(10*kelvin, 1/picosecond, 2*femtoseconds)
integrator.setConstraintTolerance(0.0001)
simulation.step(50000)

# Begin equilibration step 5: Increasing temperature to 310K
# Lasts ~2000ps total
print('Step 5 equilibration commencing...')
simulation.context.setVelocitiesToTemperature(10*kelvin)

print('Warming up the system...')
T = 310
temperature = T
integrator.setTemperature(temperature)
simulation.step(1000000)
print('Temperature of 310K reached')

# Begin production and assign barostat
print('Begin production steps...')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
system.addForce(MonteCarloBarostat(1.01325*bar, 310*kelvin))
integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 2*femtoseconds)
integrator.setConstraintTolerance(0.0001)
simulation.context.setPositions(modeller.positions)
simulation.context.setParameter('k', 0)

# Added reporters
# Includes XTC trajectory file to be used to find RMSD/RMSF vs Time
# PDB produces 30 frames, 0.5 every ns
# Production runs for 60 ns
simulation.reporters.append(XTCReporter('Cterm0623xtc.xtc', 1000))
simulation.reporters.append(StateDataReporter('Cterm0607csv.csv', 1000, time=True, temperature=True, kineticEnergy=True, potentialEnergy=True))
simulation.reporters.append(PDBReporter('Cterm0623pdb.pdb', 250000))
simulation.step(30000000)