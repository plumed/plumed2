
.PHONY: ubuntu plumed2.tgz clean fedora44 fedora44-pycv rocky8 rocky8-pycv rocky9 rocky9-pycv 

ubuntu: plumed2.tgz
	docker build --shm-size 1024M -t plumed --progress=plain .

fedora44: plumed2.tgz
	docker build --shm-size 1024M -t plumed:fedora44 --progress=plain -f fedora44 .

fedora44-pycv: plumed2.tgz
	docker build --shm-size 1024M -t plumed:fedora44-pycv --progress=plain -f fedora44-pycv .

rocky8: plumed2.tgz
	docker build --shm-size 1024M -t plumed:rocky8 --progress=plain -f rocky8 .

rocky8-pycv: plumed2.tgz
	docker build --shm-size 1024M -t plumed:rocky8-pycv --progress=plain -f rocky8-pycv .

plumed2.tgz:
	cd ../ ; git archive HEAD -o "$(CURDIR)"/plumed2.tgz --prefix plumed2/ ; cd -

clean:
	rm -f plumed2.tgz

rocky9: plumed2.tgz
	docker build --shm-size 1024M -t plumed:rocky9 --progress=plain -f rocky9 .

rocky9-pycv: plumed2.tgz
	docker build --shm-size 1024M -t plumed:rocky9-pycv --progress=plain -f rocky9-pycv .
