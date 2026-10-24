"""refresh_inputs.py -- copy the regenerated stage-2 exports into the paper's vendored inputs.

Single direction of data flow: ultrasonic/results (regenerated) -> paper_correction/vendor/results.
Every copied file is verified by SHA-256 after copying (recorded in vendor/results/INPUTS_SHA256.json).
"""
import hashlib, json, os, shutil
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "ultrasonic", "results")
DST = os.path.join(HERE, "vendor", "results")
FILES = ["groups.json", "scans.json", "results.json", "vso_scan.json"]
os.makedirs(DST, exist_ok=True)
rec = {}
for f in FILES:
    s = os.path.join(SRC, f)
    d = os.path.join(DST, f)
    shutil.copyfile(s, d)
    h_s = hashlib.sha256(open(s, "rb").read()).hexdigest()
    h_d = hashlib.sha256(open(d, "rb").read()).hexdigest()
    assert h_s == h_d, f"copy verification failed for {f}"
    rec[f] = h_d
json.dump(rec, open(os.path.join(DST, "INPUTS_SHA256.json"), "w"), indent=1)
print("refreshed", len(FILES), "inputs (sha256 verified)")
