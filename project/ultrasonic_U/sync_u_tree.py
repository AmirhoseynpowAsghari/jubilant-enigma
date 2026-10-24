"""sync_u_tree.py -- regenerate ultrasonic_U/*_u.py from ultrasonic/*.py (module renames only).

Keeps the U-family stage-2 copy physically identical to the audited main tree.  Run:
    python3 sync_u_tree.py     (rewrites the *_u.py files; then run `make test`)
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
MAIN = os.path.join(HERE, "..", "ultrasonic")
PAIRS = {"absorption.py": "absorption_u.py", "dos_ua.py": "dos_ua_u.py", "kgrid_ua.py": "kgrid_ua_u.py",
         "parameters_ua.py": "parameters_ua_u.py", "physics_ua.py": "physics_ua_u.py",
         "scanning_ua.py": "scanning_ua_u.py", "io_utils_ua.py": "io_utils_ua_u.py",
         "sc_loader.py": "sc_loader_u.py"}
# NOTE: main_ua_u.py and plotting_ua_u.py are U-FAMILY ENTRY POINTS, not copies of the main
# entry points; they are NOT synced (restored from git HEAD after an accidental overwrite).
REN = [("from parameters_ua import", "from parameters_ua_u import"), ("from dos_ua import", "from dos_ua_u import"),
       ("from kgrid_ua import", "from kgrid_ua_u import"), ("from absorption import", "from absorption_u import"),
       ("from physics_ua import", "from physics_ua_u import"), ("from scanning_ua import", "from scanning_ua_u import"),
       ("import io_utils_ua as io", "import io_utils_ua_u as io"), ("import sc_loader\n", "import sc_loader_u as sc_loader\n"),
       ("import sc_loader, io_utils_ua as io", "import sc_loader_u as sc_loader, io_utils_ua_u as io"),
       ("from sc_loader import", "from sc_loader_u import"), ("import scanning_ua as su", "import scanning_ua_u as su")]
HDR = ("# AUTO-SYNCED from project/ultrasonic/{src} by sync_u_tree.py.  Only module names are renamed\n"
       "# (*_u); the physics is identical to the main tree.  Do not edit by hand.\n")
if __name__ == "__main__":
    for src, dst in PAIRS.items():
        s = open(os.path.join(MAIN, src)).read()
        for a, b in REN:
            s = s.replace(a, b)
        s = s.replace("from __future__ import annotations", HDR.format(src=src) + "from __future__ import annotations", 1)
        open(os.path.join(HERE, dst), "w").write(s)
        print("synced", src, "->", dst)
