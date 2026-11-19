"""
Write every design's binder into the 1GC1 reference frame (gp120 core + CD4 D1D2),
so all designs, CD4, and gp120 share one coordinate system for figures and for the
co-binding analysis.

Outputs into 04_figures/scene/:
  gp120_ref.pdb   gp120 core from 1GC1 (chain G)
  cd4_ref.pdb     CD4 D1D2 from 1GC1 (chain C)
  binder_<tag>.pdb  each design, superposed via the CD4 epitope
"""
import os, glob, warnings
import numpy as np
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, PDBIO, Select

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SDIR = os.path.join(PROJECT, "paper_analysis", "01_structures")
SCENE = os.path.join(PROJECT, "paper_analysis", "04_figures", "scene")
os.makedirs(SCENE, exist_ok=True)

# superposition quality gate (A); above this the fit is unusable
MAX_FIT_RMSD = 10.0

parser = PDBParser(QUIET=True)
io = PDBIO()


class ChainSel(Select):
    def __init__(self, cid): self.cid = cid
    def accept_chain(self, c): return c.id == self.cid
    def accept_residue(self, r): return r.id[0] == " "


def ca_map(chain):
    return {r.id[1]: r["CA"].coord for r in chain if r.id[0] == " " and "CA" in r}


def superpose(mob, ref):
    mc, rc = mob.mean(0), ref.mean(0)
    m, r = mob - mc, ref - rc
    V, S, Wt = np.linalg.svd(m.T @ r)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    return R, rc - (mc @ R)


ref_st = parser.get_structure("ref", os.path.join(SDIR, "1GC1.pdb"))
ref_model = next(iter(ref_st))
ref_gp120, ref_cd4 = ref_model["G"], ref_model["C"]
ref_ca = ca_map(ref_gp120)
cd4_xyz = np.array([a.coord for a in ref_cd4.get_atoms() if a.element != "H"])

cd4_epitope = set()
for r in ref_gp120:
    if r.id[0] != " ":
        continue
    ra = np.array([a.coord for a in r if a.element != "H"])
    if ra.size and np.linalg.norm(ra[:, None, :] - cd4_xyz[None, :, :], axis=-1).min() < 5.0:
        cd4_epitope.add(r.id[1])

io.set_structure(ref_st)
io.save(os.path.join(SCENE, "gp120_ref.pdb"), ChainSel("G"))
io.save(os.path.join(SCENE, "cd4_ref.pdb"), ChainSel("C"))
print(f"wrote gp120_ref.pdb, cd4_ref.pdb (CD4 epitope {len(cd4_epitope)} residues)")

for f in sorted(glob.glob(os.path.join(SDIR, "*_split.pdb"))):
    tag = os.path.basename(f).replace("_split.pdb", "")
    st = parser.get_structure(tag, f)
    m = next(iter(st))
    if "A" not in m or "B" not in m:
        continue
    gp = ca_map(m["A"])
    common = sorted(set(gp) & set(ref_ca))
    local = [i for i in common if i in cd4_epitope]
    use = local if len(local) >= 15 else common
    if len(use) < 15:
        print(f"[skip] {tag}: too few shared residues")
        continue
    mob = np.array([gp[i] for i in use])
    rf = np.array([ref_ca[i] for i in use])
    R, t = superpose(mob, rf)
    rms = float(np.sqrt((((mob @ R) + t - rf) ** 2).sum(1).mean()))
    if rms > MAX_FIT_RMSD:
        print(f"[skip] {tag}: fit RMSD {rms:.1f} A exceeds {MAX_FIT_RMSD} A")
        continue

    for a in m["B"].get_atoms():
        a.set_coord((a.coord @ R) + t)
    io.set_structure(st)
    io.save(os.path.join(SCENE, f"binder_{tag}.pdb"), ChainSel("B"))
    print(f"[ok] binder_{tag}.pdb (fit {rms:.2f} A)")
