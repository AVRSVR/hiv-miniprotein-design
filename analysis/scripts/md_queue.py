"""
Sequential MD queue for the design panel.

Runs as a second stream alongside the in-flight S310 job: these systems are small
(~8k atoms in implicit solvent) and leave the GPU underutilised, so two concurrent
contexts raise aggregate throughput rather than halving it.

docked_complex is included here even though it was excluded from the CD4-competition
analysis - that exclusion was a superposition failure against 1GC1, which does not
affect whether its own interface is stable.

Usage: python md_queue.py [tag:ns ...]
"""
import os, sys, time, subprocess

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
PY = r"C:\Users\anant\miniconda3\python.exe"
RUN = os.path.join(PROJECT, "paper_analysis", "scripts", "run_md.py")
SPLIT = os.path.join(PROJECT, "paper_analysis", "01_structures")
MD = os.path.join(PROJECT, "paper_analysis", "03_md")

DEFAULT_JOBS = [
    ("Chai1_Best", 3.0),     # D2 de novo bispecific
    ("Latest_Best", 3.0),    # D5 de novo, coreceptor-site blocker
    ("model_4_chai", 2.0),   # D4 de novo
    ("Cyclic_Ring", 2.0),    # D3 CD4-mimetic control
    ("docked_complex", 2.0), # unassigned epitope; stability only
]


def already_complete(tag):
    """True only if a finished run exists - a checkpoint does not count."""
    p = os.path.join(MD, f"{tag}_md_metrics.json")
    if not os.path.exists(p):
        return False
    try:
        import json
        return bool(json.load(open(p)).get("complete", True))
    except Exception:
        return False


if __name__ == "__main__":
    jobs = DEFAULT_JOBS
    if len(sys.argv) > 1:
        jobs = []
        for a in sys.argv[1:]:
            t, _, n = a.partition(":")
            jobs.append((t, float(n or 2.0)))

    for tag, ns in jobs:
        src = os.path.join(SPLIT, f"{tag}_split.pdb")
        if not os.path.exists(src):
            print(f"skip {tag}: no split PDB", flush=True)
            continue
        if already_complete(tag):
            print(f"skip {tag}: already complete", flush=True)
            continue
        log = os.path.join(MD, f"{tag}_run.log")
        print(f"=== {tag}: {ns} ns ===  {time.strftime('%H:%M:%S')}", flush=True)
        with open(log, "w") as fh:
            rc = subprocess.call([PY, "-u", RUN, src, tag, str(ns)],
                                 stdout=fh, stderr=subprocess.STDOUT)
        print(f"{tag} finished rc={rc}  {time.strftime('%H:%M:%S')}", flush=True)
    print("QUEUE COMPLETE", flush=True)
