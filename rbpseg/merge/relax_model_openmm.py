# import packages
import numpy as np
import sys
import os
from openmm import unit
from openmm.app.internal.pdbstructure import PdbStructure
import openmm.app as omm_app
import openmm as omm
import openmm.unit as unt
from openmm.app.forcefield import PME, HBonds
#import pdbfixer
#import io

def relax_amber14_score(pdb_path, relaxed_pdb_path, minimize_energy=True,local_min=False, iter=0, tol=0.1,pltfm='CUDA'):

    fixer = pdbfixer.PDBFixer(pdb_path)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.removeHeterogens()
    #print(fixer.missingTerminals)

    omm_app.pdbfile.PDBFile.writeFile(fixer.topology, fixer.positions, open(relaxed_pdb_path, 'w'))

    pdb = omm_app.Modeller(fixer.topology, fixer.positions)
    forcefield = omm_app.ForceField('amber14-all.xml')
    pdb.addHydrogens(forcefield)
    
    platform = omm.Platform.getPlatformByName(pltfm)  # OpenCL works out the box on a Mac platform = omm.Platform.getPlatformByName('CPU')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=omm_app.CutoffNonPeriodic)
    integrator = omm.LangevinIntegrator(298*unt.kelvin, 1.0/unt.picoseconds, 2.0*unt.femtoseconds)
    simulation = omm_app.Simulation(pdb.topology, system, integrator, platform=platform)
    simulation.context.setPositions(pdb.positions)
    simulation.reporters.append(omm_app.PDBReporter('output.pdb', 1000))
    simulation.reporters.append(
        omm_app.StateDataReporter(
            sys.stdout, 1000,
            step=True, potentialEnergy=True, temperature=True
        )
    )

    if minimize_energy: 
        if local_min: 
            omm.LocalEnergyMinimizer.minimize(simulation,maxIterations=iter,tolerance=tol)
        else: 
            simulation.minimizeEnergy(maxIterations=iter,tolerance=tol)
            
    current_state = simulation.context.getState(getEnergy=True,getPositions=True)
    
    top = simulation.topology
    pos=current_state.getPositions()
    if minimize_energy:
        omm_app.pdbfile.PDBFile.writeFile(top, pos, open(relaxed_pdb_path, 'w'))

    pot_en = current_state.getPotentialEnergy()
    kin_en = current_state.getKineticEnergy()
    
    return [pot_en.value_in_unit(pot_en.unit),kin_en.value_in_unit(kin_en.unit)]

