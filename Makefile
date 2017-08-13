# Makefile for all programs under this folder
CP = /bin/cp -f 

all:
	make -f removeUnnecessaryGap.makefile
clean:
	make clean -f removeUnnecessaryGap.makefile
