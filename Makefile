# ---------------------------------------------------------------------
#  move-motor
#  Copyright (C) William Dickson, 2008.
#  
#  wbd@caltech.edu
#  www.willdickson.com
#
#  Released under the LGPL Licence, Version 3
#  
#  This file is part of move-motor.
#
#  move-motor is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#    
#  move-motor is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with move-motor.  If not, see
#  <http://www.gnu.org/licenses/>.
#
# ----------------------------------------------------------------------
#  Makefile
#
#  Purpose: to build move-motor library and install python ctypes 
#  interface.
#
#  Author: Will Dickson
# ---------------------------------------------------------------------- 
PWD=$(shell pwd)
LIB_DIR=$(PWD)/lib
PYLIB_DIR=$(PWD)/python

default: libmove-motor

libmove-motor:
	echo $(LIB_DIR)
	$(MAKE) -C $(LIB_DIR)

.PHONY: clean install uninstall

clean:
	$(MAKE) clean -C $(LIB_DIR)
	-rm *~

install:
	$(MAKE) install -C $(LIB_DIR)
	cd $(PYLIB_DIR); python $(PYLIB_DIR)/setup.py develop; cd $(PWD)

uninstall:
	$(MAKE) uninstall -C $(LIB_DIR)
	cd $(PYLIB_DIR); python $(PYLIB_DIR)/setup.py develop --uninstall; cd $(PWD)
