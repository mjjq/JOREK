#/usr/bin/env python
import glob
from FRUIT import *

# Automatically detect files in tests/fruit and create a single driver
test_modules = glob.glob('particles/tests/*.f90')
driver = "particles/tests/particle_test_driver.f90"
if driver in test_modules: test_modules.remove(driver)

suite = test_suite(test_modules)
suite.build_run(driver, build_command="make particle_test_driver")
suite.summary()
