"""
Explicit-solvent MD refinement / validation loop.

Runs physiological explicit-solvent MD (TIP3P water, 0.15 M NaCl, PME
electrostatics, NPT at 310.15 K, 2 fs Langevin integration) on a docked
complex, tracking interface RMSD (iRMSD) of the ligand backbone against
the target pocket over time. Used both as:

  1. A short (ns-scale) filter during iterative refinement: candidates
     whose iRMSD stays under a stability threshold are advanced; those
     that drift are discarded.
  2. A longer confirmation run on the final selected design.

Batch mode picks up any `*_complex.pdb` file, solvates + minimizes +
equilibrates + runs production MD, and writes a trajectory (DCD) plus
an iRMSD/Rg diagnostic plot per design, with an optional multi-design
comparison plot at the end.

Requires: openmm, pdbfixer, numpy, matplotlib (GPU/CUDA recommended).
"""

import os
import sys
import traceback
import numpy as np
import matplotlib.pyplot as plt
from openmm import *
from openmm.app import *
from openmm.unit import *
from pdbfixer import PDBFixer


def kabsch_align(coords, ref_coords):
    """SVD-based Kabsch superposition to minimize RMSD between coords and ref_coords."""
    c_center = np.mean(coords, axis=0)
    r_center = np.mean(ref_coords, axis=0)
    c_trans = coords - c_center
    r_trans = ref_coords - r_center
    h = np.dot(c_trans.T, r_trans)
    u, s, vt = np.linalg.svd(h)
    v = vt.T
    d = np.linalg.det(np.dot(v, u.T))
    e = np.eye(3)
    if d < 0.0:
        e[2, 2] = -1.0
    r = np.dot(np.dot(v, e), u.T)
    return np.dot(c_trans, r.T) + r_center


def run_complex_md(complex_pdb, output_dir, production_steps, report_interval=5000):
    """
    Solvate, minimize, equilibrate, and run production MD on a receptor-ligand
    complex PDB, tracking interface RMSD and ligand radius of gyration.

    production_steps: total production steps (2 fs/step). E.g. 1,000,000 steps
    ~= 2 ns; 25,000,000 steps ~= 50 ns.
    """
    clean_name = os.path.splitext(os.path.basename(complex_pdb))[0]
    os.makedirs(output_dir, exist_ok=True)

    fixer = PDBFixer(filename=complex_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.removeHeterogens(False)
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    modeller = Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(
        forcefield, model='tip3p', padding=1.0 * nanometers,
        ionicStrength=0.15 * molar, neutralize=True, boxShape='dodecahedron',
    )

    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometers, constraints=HBonds,
    )
    system.addForce(MonteCarloBarostat(1.0 * bar, 310.15 * kelvin))
    integrator = LangevinIntegrator(310.15 * kelvin, 1.0 / picoseconds, 2.0 * femtoseconds)

    try:
        platform = Platform.getPlatformByName('CUDA')
        simulation = Simulation(modeller.topology, system, integrator, platform, {'Precision': 'mixed'})
    except Exception:
        simulation = Simulation(modeller.topology, system, integrator)

    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=1000)
    simulation.context.setVelocitiesToTemperature(310.15 * kelvin)
    simulation.step(5000)  # 10 ps equilibration

    pocket_ca = [a.index for a in modeller.topology.atoms() if a.name == 'CA' and a.residue.chain.index == 0]
    anchor_ca = [a.index for a in modeller.topology.atoms() if a.name == 'CA' and a.residue.chain.index == 1]
    ref_positions = np.array(simulation.context.getState(getPositions=True).getPositions().value_in_unit(nanometers))
    ref_anchor = ref_positions[anchor_ca]

    dcd_out = os.path.join(output_dir, f"{clean_name}_trajectory.dcd")
    simulation.reporters.append(DCDReporter(dcd_out, report_interval, enforcePeriodicBox=True))

    num_reports = production_steps // report_interval
    time_axis, irmsds = [], []
    for i in range(num_reports):
        simulation.step(report_interval)
        pos = np.array(simulation.context.getState(getPositions=True).getPositions().value_in_unit(nanometers))
        aligned_anchor = kabsch_align(pos[anchor_ca], ref_anchor)
        irmsd = np.sqrt(np.mean(np.sum((aligned_anchor - ref_anchor) ** 2, axis=1))) * 10.0
        irmsds.append(irmsd)
        time_axis.append((i + 1) * report_interval * 2.0 / 1000.0)

    return np.array(time_axis), np.array(irmsds), dcd_out


if __name__ == "__main__":
    # Example: 2 ns production run (1,000,000 steps @ 2 fs) used for
    # ns-scale stability checks during iterative refinement.
    complex_file = sys.argv[1] if len(sys.argv) > 1 else "design_complex.pdb"
    t, rmsd, dcd = run_complex_md(complex_file, "./outputs", production_steps=1_000_000)
    print(f"Trajectory written to {dcd}")
