"""manifest.py -- SHA-256 manifest of every regenerated JSON (for reproducibility checks)."""
import glob, hashlib, json, os
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SKIP = ("PROVENANCE",)
out = {}
for pat in ("project/sc_params/results/*.json", "project/ultrasonic/results/*.json",
            "project/ultrasonic_U/results/*.json"):
    for f in sorted(glob.glob(os.path.join(ROOT, pat))):
        if any(s in os.path.basename(f) for s in SKIP):
            continue
        out[os.path.relpath(f, ROOT)] = hashlib.sha256(open(f, "rb").read()).hexdigest()
os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
with open(os.path.join(ROOT, "results", "REGEN_MANIFEST.json"), "w") as fh:
    json.dump(dict(n_files=len(out), sha256=out), fh, indent=1)
print(f"manifest: {len(out)} JSON files hashed -> results/REGEN_MANIFEST.json")
