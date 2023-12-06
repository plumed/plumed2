
.PHONY: ubuntu plumed2.tgz clean fedora38 fedora38-pycv rocky8 rocky8-pycv 

ubuntu: plumed2.tgz
	docker build -t plumed --progress=plain .

fedora38: plumed2.tgz
	docker build -t plumed:fedora38 --progress=plain -f fedora38 .
fedora38-pycv: plumed2.tgz
	docker build -t plumed:fedora38-pycv --progress=plain -f fedora38-pycv .


rocky8: plumed2.tgz
	docker build -t plumed:rocky8 --progress=plain -f rocky8 .

rocky8-pycv: plumed2.tgz
	docker build -t plumed:rocky8-pycv --progress=plain -f rocky8-pycv .

plumed2.tgz:
	cd ../ ; git archive HEAD -o "$(CURDIR)"/plumed2.tgz --prefix plumed2/ ; cd -

clean:
	rm -f plumed2.tgz
