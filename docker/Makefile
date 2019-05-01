
.PHONY: all plumed2.tgz clean

all: plumed2.tgz
	docker build -t plumed .

plumed2.tgz:
	cd ../ ; git archive HEAD -o "$(CURDIR)"/plumed2.tgz --prefix plumed2/ ; cd -

clean:
	rm -f plumed2.tgz
