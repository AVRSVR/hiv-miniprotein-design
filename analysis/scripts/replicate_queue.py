"""
Independent replicate trajectories for the two designs the paper's claims rest on.

Each run re-randomises initial velocities, so replicates are statistically
independent samples of the same ensemble. Three replicates per design lets the
stability metrics be reported as mean +/- SD across trajectories instead of as a
single number with unquantified run-to-run variance.
"""
import os, subprocess, time

PROJECT = r"C:\Users\anant\Desktop\hiv_miniprotein_design"
PY = r"C:\Users\anant\miniconda3\python.exe"
RUN = os.path.join(PROJECT, "paper_analysis", "scripts", "run_md.py")
SPLIT = os.path.join(PROJECT, "paper_analysis", "01_structures")
MD = os.path.join(PROJECT, "paper_analysis", "03_md")

NS = 1.5
JOBS = [
    ("Latest_Best", "Latest_Best_rep2"),   # D5 - lead design
    ("Latest_Best", "Latest_Best_rep3"),
    ("S310", "S310_rep2"),                 # D1 - CD4-mimetic control
    ("S310", "S310_rep3"),
    ("Chai1_Best", "Chai1_Best_rep2"),     # D2
]

if __name__ == "__main__":
    for src, tag in JOBS:
        out = os.path.join(MD, f"{tag}_md_metrics.json")
        if os.path.exists(out):
            print(f"skip {tag}", flush=True)
            continue
        print(f"=== {tag}: {NS} ns ===  {time.strftime('%H:%M:%S')}", flush=True)
        with open(os.path.join(MD, f"{tag}_run.log"), "w") as fh:
            rc = subprocess.call(
                [PY, "-u", RUN, os.path.join(SPLIT, f"{src}_split.pdb"), tag, str(NS)],
                stdout=fh, stderr=subprocess.STDOUT)
        print(f"{tag} finished rc={rc}  {time.strftime('%H:%M:%S')}", flush=True)
    print("REPLICATES COMPLETE", flush=True)
