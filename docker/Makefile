
.PHONY: ubuntu plumed2.tgz clean fedora42 fedora42-pycv rocky8 rocky8-pycv 

ubuntu: plumed2.tgz
	docker build --shm-size 1024M -t plumed --progress=plain .

fedora42: plumed2.tgz
	docker build --shm-size 1024M -t plumed:fedora42 --progress=plain -f fedora42 .
fedora42-pycv: plumed2.tgz
	docker build --shm-size 1024M -t plumed:fedora42-pycv --progress=plain -f fedora42-pycv .


rocky8: plumed2.tgz
	docker build --shm-size 1024M -t plumed:rocky8 --progress=plain -f rocky8 .

rocky8-pycv: plumed2.tgz
	docker build --shm-size 1024M -t plumed:rocky8-pycv --progress=plain -f rocky8-pycv .

plumed2.tgz:
	cd ../ ; git archive HEAD -o "$(CURDIR)"/plumed2.tgz --prefix plumed2/ ; cd -

clean:
	rm -f plumed2.tgz
