# Makefile for all programs under this folder
CP = /bin/cp -f 

APPS = removeUnnecessaryGap

.PHONY: all clean $(APPS)

all: $(APPS)

$(APPS):
	make -f $@.makefile

clean:
	$(foreach var,$(APPS), make clean -f  $(var).makefile;)
