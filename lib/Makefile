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
	$(CC) $(LXRT_CFLAGS) $? -o test-$(LIBRARY_NAME)

.PHONY: clean install uninstall

clean: 
	if [ -f $(TARGET).$(RELEASE_NUM) ]; then rm $(TARGET).$(RELEASE_NUM); fi
	if [ -f test-$(LIBRARY_NAME) ]; then rm test-$(LIBRARY_NAME); fi
	if [ -f *.o ]; then rm *.o; fi
	if [ -f *~ ]; then rm *~; fi

install:
	cp $(TARGET).$(RELEASE_NUM) $(INSTALL_DIR)
	if [ ! -h $(INSTALL_DIR)/$(TARGET) ]; then ln -s ${INSTALL_DIR}/$(TARGET).$(RELEASE_NUM) ${INSTALL_DIR}/$(TARGET); fi
	ldconfig

uninstall:
	rm ${INSTALL_DIR}/$(TARGET)
	rm ${INSTALL_DIR}/$(TARGET).$(RELEASE_NUM)
	ldconfig	