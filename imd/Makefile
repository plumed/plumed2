-include ../Makefile.conf

all: Makefile.conf libimd.$(SOEXT)

Makefile.conf: libimd.$(SOEXT)
	echo IMDLIB=-L$(PWD) -limd     >  Makefile.conf
	echo IMDDEF=-D__PLUMED_HAS_IMD >> Makefile.conf


libimd.$(SOEXT): imd.o vmdsock.o
	$(LDSO) -o $@ $^ -lnsl

vmdsock.c imd.c vmdsock.h imd.h: imdapi.tar.gz
	tar xzOf imdapi.tar.gz imdapi/$@ > $@

vmdsock.o imd.o: vmdsock.h imd.h

imdapi.tar.gz:
	wget http://www.ks.uiuc.edu/Research/vmd/imd/code/imdapi.tar.gz


clean:
	rm -fr imdapi.tar.gz *.c *.h *.$(SOEXT) Makefile.conf *.o




