PWD=$(shell pwd)
LIB_DIR=$(PWD)/lib
PYLIB_DIR=$(PWD)/python/pylibmove-motor

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
