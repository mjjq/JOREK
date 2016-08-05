#/usr/bin/env python
import glob
from FRUIT import *

# Automatically detect files in tests/fruit and create a single driver
test_modules = glob.glob('tests/fruit/*.f90')
driver = "tests/fruit_test_driver.f90"

suite = test_suite(test_modules)
suite.build_run(driver, build_command="make fruit_test_driver")
suite.summary()
