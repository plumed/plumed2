
.PHONY: ubuntu plumed2.tgz clean

ubuntu: plumed2.tgz
	docker build -t plumed .

centos7: plumed2.tgz
	docker build -t plumed -f centos7 .

plumed2.tgz:
	cd ../ ; git archive HEAD -o "$(CURDIR)"/plumed2.tgz --prefix plumed2/ ; cd -

clean:
	rm -f plumed2.tgz
