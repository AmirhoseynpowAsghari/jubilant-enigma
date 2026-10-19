# minimal reproducibility entry points (relative paths only)
PY ?= python3
test:
	cd project/ultrasonic && $(PY) validate_ua.py | tail -4
	cd project/ultrasonic && $(PY) case2_checks.py | tail -3
	cd project/ultrasonic && $(PY) levy1963_check.py | tail -3
	cd project/ultrasonic && $(PY) dualspace_check.py | tail -3
master:
	cd project/ultrasonic && $(PY) main_ua.py master
figures:
	cd project/ultrasonic && $(PY) plotting_ua.py
