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
#  Purpose: to build move-motor library
#
#  Author: Will Dickson
# ---------------------------------------------------------------------- 
LIBRARY_NAME=move-motor
MINOR_NUM=1
RELEASE_NUM=0.1
INSTALL_DIR=/usr/local/lib

SRC=${LIBRARY_NAME}.c
OBJ=${LIBRARY_NAME}.o
TARGET=lib${LIBRARY_NAME}.so.${MINOR_NUM}

CC= $(shell rtai-config --cc)
LXRT_CFLAGS= $(shell rtai-config --lxrt-cflags)
LXRT_LDFLAGS= $(shell rtai-config --lxrt-ldflags)
INCLUDES= -I 
CFLAGS= -fPIC -c

$(TARGET): $(OBJ)
	${CC} -shared -W1,-soname,$@ -o $@.$(RELEASE_NUM) $? -lm 

$(OBJ): $(SRC)	
	$(CC) $(CFLAGS) $(LXRT_CFLAGS) $? -o $@ 

test: $(LIBRARY_NAME).c test-$(LIBRARY_NAME).c
	$(CC) $(LXRT_CFLAGS) $? -o test-$(LIBRARY_NAME) $(LXRT_LDFLAGS)

.PHONY: clean install uninstall

clean: 
	if [ -f $(TARGET).$(RELEASE_NUM) ]; then rm $(TARGET).$(RELEASE_NUM); fi
	if [ -f test-$(LIBRARY_NAME) ]; then rm test-$(LIBRARY_NAME); fi
	-rm *.o
	-rm *~

install:
	cp $(TARGET).$(RELEASE_NUM) $(INSTALL_DIR)
	if [ ! -h $(INSTALL_DIR)/$(TARGET) ]; then ln -s ${INSTALL_DIR}/$(TARGET).$(RELEASE_NUM) ${INSTALL_DIR}/$(TARGET); fi
	ldconfig

uninstall:
	rm ${INSTALL_DIR}/$(TARGET)
	rm ${INSTALL_DIR}/$(TARGET).$(RELEASE_NUM)
	ldconfig	