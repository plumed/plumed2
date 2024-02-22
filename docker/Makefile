
.PHONY: ubuntu plumed2.tgz clean centos7 fedora34

ubuntu: plumed2.tgz
	docker build --shm-size 1024M -t plumed .

centos7: plumed2.tgz
	docker build --shm-size 1024M -t plumed -f centos7 .

fedora34: plumed2.tgz
	docker build --shm-size 1024M -t plumed -f fedora34 .

plumed2.tgz:
	cd ../ ; git archive HEAD -o "$(CURDIR)"/plumed2.tgz --prefix plumed2/ ; cd -

clean:
	rm -f plumed2.tgz
