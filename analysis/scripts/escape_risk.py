"""
Escape-risk map: which gp120 positions can the virus cheaply mutate to break a binder?

The escape argument in this paper is structural, so it should be tested
structurally. Two quantities decide whether a contact is an escape liability:

  dependence  how much binding is lost when that gp120 sidechain is removed,
              from the same decoupling alanine scan applied to the receptor side
  constraint  how conserved the position is across HIV-1 Env sequences

A contact the binder leans on heavily at a position the virus varies freely is an
escape route. A contact it leans on at an invariant position is escape-resistant,
because the substitution that would break the binder also costs the virus.

  escape_risk = (ddG normalised to the design's maximum) x (1 - conservation)

High-risk positions are the ones a redesign should move off; escape-resistant
contacts are what a durable binder should be built on.

Usage: python escape_risk.py <tag> [n_frames]
"""
import os, sys, json, warnings
from collections import Counter
import numpy as np
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mdtraj as md
from Bio import SeqIO, Align
from alanine_scan import (build_wt, make_sim, Decoupler, energy, BACKBONE,
                          SKIP_RES, IFACE_CUT, LABEL)
from openmm.app import Modeller

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
SPLIT = os.path.join(PROJECT, "paper_analysis", "01_structures")
MD = os.path.join(PROJECT, "paper_analysis", "03_md")
OUT = os.path.join(PROJECT, "paper_analysis", "05_tables")
# Superseded: the 47-sequence reviewed set below is no longer the source of
# conservation. conservation_core reads hiv_env_full.fasta (4,497 full-length
# UniProt records) and caches the per-position result. Kept only so the HXB2
# reference sequence that follows stays next to its provenance.
_UNUSED_FASTA = os.path.join(SPLIT, "hiv_env.fasta")

HXB2_GP120 = (
    "TEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQ"
    "MHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFF"
    "YKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTH"
    "GIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGK"
    "IGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTW"
    "FNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIF"
    "RPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKR"
)
HXB2_START = 31


def conservation_by_position():
    """Per-HXB2-position conservation from the shared cache (conservation_core).

    Previously this re-implemented the alignment locally, duplicating
    conservation.py. Both now read one cache built from the full-length UniProt
    Env set, clustered at 98% identity so an over-sequenced lineage cannot
    dominate a column.
    """
    import conservation_core
    cons, meta = conservation_core.load()
    print("conservation: %d sequences, %d clusters, %d positions"
          % (meta["sequences_used"], meta["clusters_at_98pct_identity"],
             meta["positions_scored"]), flush=True)
    return cons


def main(tag, n_frames=10):
    cons = conservation_by_position()
    print(f"conservation computed for {len(cons)} gp120 positions", flush=True)

    split_pdb = os.path.join(SPLIT, f"{tag}_split.pdb")
    top, pos = build_wt(split_pdb)
    ref = os.path.join(MD, f"{tag}_mmgbsa_top.pdb")
    traj = md.load(os.path.join(MD, f"{tag}_traj.dcd"), top=ref)

    atoms = list(top.atoms())
    rec_i = np.array([i for i, a in enumerate(atoms) if a.residue.chain.id == "A"])
    lig_i = np.array([i for i, a in enumerate(atoms) if a.residue.chain.id == "B"])
    idx = np.linspace(traj.n_frames // 2, traj.n_frames - 1, n_frames).astype(int)
    X = traj.xyz[idx]

    mod_l = Modeller(top, pos)
    mod_l.delete([r for r in top.residues() if r.chain.id != "B"])
    mod_r = Modeller(top, pos)
    mod_r.delete([r for r in top.residues() if r.chain.id != "A"])

    sim_c = make_sim(top, pos)
    sim_l = make_sim(mod_l.topology, mod_l.positions)
    sim_r = make_sim(mod_r.topology, mod_r.positions)
    rec_pos_in_rec = {g: k for k, g in enumerate(rec_i)}

    e_lig = np.array([energy(sim_l, x[lig_i]) for x in X])
    wt = np.array([energy(sim_c, x) for x in X]) - \
         np.array([energy(sim_r, x[rec_i]) for x in X]) - e_lig
    print(f"{tag}: wild-type dG = {wt.mean():.1f} kJ/mol", flush=True)

    # gp120 interface residues
    x0 = X[len(X) // 2]
    heavy_l = [i for i in lig_i if atoms[i].element.symbol != "H"]
    iface = set()
    for i in rec_i:
        if atoms[i].element.symbol == "H":
            continue
        if np.linalg.norm(x0[heavy_l] - x0[i], axis=1).min() < IFACE_CUT:
            iface.add(atoms[i].residue.index)
    targets = [r for r in sorted(iface)
               if next(a for a in atoms if a.residue.index == r).residue.name not in SKIP_RES]
    print(f"  {len(iface)} gp120 interface residues, {len(targets)} scannable", flush=True)

    dec_c = Decoupler(sim_c)
    dec_r = Decoupler(sim_r)

    rows = []
    for n, rid in enumerate(targets):
        res = next(a.residue for a in atoms if a.residue.index == rid)
        side_c = [a.index for a in res.atoms() if a.name not in BACKBONE]
        if not side_c:
            continue
        side_r = [rec_pos_in_rec[i] for i in side_c if i in rec_pos_in_rec]

        dec_c.zero(side_c)
        dec_r.zero(side_r)
        mut = np.array([energy(sim_c, x) for x in X]) - \
              np.array([energy(sim_r, x[rec_i]) for x in X]) - e_lig
        dec_c.restore()
        dec_r.restore()

        ddg = float((mut - wt).mean())
        c = cons.get(int(res.id))
        rows.append({"resname": res.name, "hxb2_pos": int(res.id),
                     "ddG_kcal_mol": round(ddg / 4.184, 2),
                     "conservation": round(c, 3) if c is not None else None})
        print(f"  [{n+1}/{len(targets)}] {res.name}{res.id}: "
              f"{ddg/4.184:+.2f} kcal/mol, cons "
              f"{'n/a' if c is None else f'{100*c:.0f}%'}", flush=True)

    # escape risk: heavy dependence on a variable position
    scored = [r for r in rows if r["conservation"] is not None and r["ddG_kcal_mol"] > 0]
    mx = max((r["ddG_kcal_mol"] for r in scored), default=1.0) or 1.0
    for r in rows:
        if r["conservation"] is None or r["ddG_kcal_mol"] <= 0:
            r["escape_risk"] = None
            continue
        r["escape_risk"] = round((r["ddG_kcal_mol"] / mx) * (1.0 - r["conservation"]), 3)

    ranked = sorted([r for r in rows if r["escape_risk"] is not None],
                    key=lambda z: -z["escape_risk"])
    resistant = [r for r in rows
                 if r["conservation"] is not None and r["conservation"] >= 0.95
                 and r["ddG_kcal_mol"] >= 1.0]

    out = {"design": tag, "label": LABEL.get(tag, tag), "n_frames": len(X),
           "wt_dG_kJ_mol": round(float(wt.mean()), 1),
           "contacts": rows,
           "top_escape_risk": ranked[:8],
           "escape_resistant_contacts": sorted(resistant,
                                               key=lambda z: -z["ddG_kcal_mol"])}
    path = os.path.join(OUT, f"escape_risk_{tag}.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)

    print(f"\n{tag} top escape-risk positions:")
    for r in ranked[:6]:
        print(f"  {r['resname']}{r['hxb2_pos']}: ddG {r['ddG_kcal_mol']:+.2f}, "
              f"cons {100*r['conservation']:.0f}%, risk {r['escape_risk']:.2f}")
    print(f"\n{tag} escape-resistant contacts (>=95% conserved, >=1 kcal/mol): "
          f"{len(resistant)}")
    for r in resistant[:6]:
        print(f"  {r['resname']}{r['hxb2_pos']}: ddG {r['ddG_kcal_mol']:+.2f}, "
              f"cons {100*r['conservation']:.0f}%")
    print(f"wrote {path}")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 10)
