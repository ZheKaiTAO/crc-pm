"""Run from any working directory; every writable artifact stays in approved roots."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/crc-pm-bulk-mpl")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import argparse
from pipeline import Pipeline

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["all", "prepare", "differential", "validate", "survival", "enrichment", "figures", "report"], default="all")
    args = parser.parse_args()
    pipeline = Pipeline()
    stages = ["prepare", "differential", "validate", "survival", "enrichment", "report"] if args.stage == "all" else [args.stage]
    for stage in stages:
        print(f"Starting {stage}", flush=True)
        getattr(pipeline, stage)()
        print(f"Finished {stage}", flush=True)
