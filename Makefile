# Reproducibility entry points.  Every recipe runs under bash with pipefail and errexit, so a
# failing check or script can NOT be masked by a pipe (audit F-05).  No `|| true` anywhere.
SHELL := /bin/bash
.SHELLFLAGS := -o pipefail -c
PY ?= python3

.PHONY: test test-kernels test-fixes master figures

test: test-kernels test-fixes
	cd project/ultrasonic && $(PY) validate_ua.py
	cd project/ultrasonic && $(PY) case2_checks.py > /dev/null
	cd project/ultrasonic && $(PY) levy1963_check.py > /dev/null
	cd project/ultrasonic && $(PY) dualspace_check.py > /dev/null
	@echo "MAKE TEST: all checks passed"

test-kernels:
	cd project/ultrasonic && $(PY) kernel_bdg_check.py

test-fixes:
	cd project/tests && $(PY) test_fixes.py

master:
	cd project/ultrasonic && $(PY) main_ua.py master

figures:
	cd project/ultrasonic && $(PY) plotting_ua.py
