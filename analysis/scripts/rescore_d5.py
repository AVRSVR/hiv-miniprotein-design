"""
Score the D5 redesign variants from trajectories already on disk.

The first redesign pass ran all the MD correctly but threw every result away:
on Windows the DCDReporter's file handle was still open when mdtraj tried to
read the trajectory, so each variant failed at the scoring step with
"premature end of file". The trajectories themselves are intact, so scoring is
re-run here rather than repeating ~1 h of MD per variant.
"""
import os, sys, json, glob, warnings
import numpy as np
warnings.filterwarnings("ignore")
import mdtraj as mdt
from openmm.app import PDBFile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from redesign_d5 import score, OUT, VARIANTS

out = []
for name, muts in VARIANTS:
    dcd = os.path.join(OUT, f"{name}_traj.dcd")
    top = os.path.join(OUT, f"{name}_top.pdb")
    if not (os.path.exists(dcd) and os.path.exists(top)):
        print(f"[skip] {name}: no trajectory", flush=True)
        continue
    try:
        pdb = PDBFile(top)
        tr = mdt.load(dcd, top=top)
        keep = tr.xyz[len(tr) // 2:]
        dg, sem, dep = score(name, pdb.topology, keep)
    except Exception as e:
        print(f"[fail] {name}: {e}", flush=True)
        continue
    r = {"variant": name, "mutations": muts,
         "dG_kJ_mol": round(dg, 1), "dG_kcal_mol": round(dg / 4.184, 1),
         "sem_kcal_mol": round(sem / 4.184, 1),
         "liability_dependence_kcal": dep,
         "liability_total_kcal": round(sum(dep.values()), 2),
         "n_frames": len(keep)}
    out.append(r)
    print(f"[ok] {name:<12} dG {r['dG_kcal_mol']:>8.1f} +/- {r['sem_kcal_mol']:.1f} kcal/mol | "
          f"liability {r['liability_total_kcal']:>6.2f} | {dep}", flush=True)
    with open(os.path.join(OUT, "redesign_d5.json"), "w") as f:
        json.dump(out, f, indent=2)
print("RESCORE DONE", flush=True)
