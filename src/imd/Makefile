# include the machine dependent configuration
ifneq ($(MAKECMDGOALS),clean)
  -include ../../Makefile.conf
endif

all: libimd.$(SOEXT)

libimd.$(SOEXT): imd.o vmdsock.o ActionIMD.o
	$(LDSO) -o $@ $^ # -lnsl

vmdsock.c imd.c vmdsock.h imd.h: imdapi.tar.gz
	tar xzOf imdapi.tar.gz imdapi/$@ > $@

vmdsock.o imd.o: vmdsock.h imd.h

imdapi.tar.gz:
	wget http://www.ks.uiuc.edu/Research/vmd/imd/code/imdapi.tar.gz


clean:
	rm -fr imdapi.tar.gz *.c *.h *.$(SOEXT) *.o

obj: all

links:
	# do nothing


