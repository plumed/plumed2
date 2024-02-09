#!/bin/bash
#small utilty that recurses the src folder and creates a standard CMakeListst.txt
#for modules where is not present
#thinked to be launched in repodir/src

createCMakeLists (){
    dir=$1
    if test -f "$dir/module.type" ; then
	if test -f $dir/CMakeLists.txt ;then
        echo "$dir has the CMakeLists.txt"
	# if grep -q "automatically generated CMakeLists.txt, if it does not work" $dir/CMakeLists.txt; then 
        #this update non modified CMakeLists.txt, decomment if needed
        # rm -v $dir/CMakeLists.txt 
	# fi
    else
        echo "$dir"
	fi 

	if test ! -f $dir/CMakeLists.txt
	then
	    (   
		cd $dir ||exit
		{
		    echo "message(WARNING \"${dir} has an automatically generated CMakeLists.txt, if it does not work modify it and remove this warning\")"
		    echo "#the variable module_name is set up as a sugar to reduce \"copy-paste\" errors"
		    echo "set (module_name \"${dir}\")"
		    echo "#Note that the macros here require this directory added as a subdir of plumed/src"

		    if [[ $(wc -l < Makefile) -gt 4 ]]; then
			echo "message (FATAL_ERROR \"\${module_name} has a non standard Makefile (more than 4 lines) you need to modify the CMakeLists.txt!\")"
		    fi

		    case "$(cat "module.type")" in
			(always)      echo "set(module_\${module_name} ON CACHE INTERNAL \"always active module \${module_name}\")";;
			(default-on)  echo "option(module_\${module_name} \"activate module \${module_name}\" ON)" ;;
			(default-off) echo "option(module_\${module_name} \"activate module \${module_name}\" OFF)" ;;
		    esac
		    echo "ADDMODULETOKERNEL(\${module_name}" 
		    ls -1 *.cpp
		    echo ")"
		    if grep -q USE Makefile; then
			echo "ADDMODULENEEDS(\${module_name}" 
			t=$(awk '/USE=/{print }' < Makefile)
			echo -e "\t${t#USE*=}"
			echo ")"
			echo "ADDMODULEDEPENDENCIES(\${module_name}" 
			t=$(awk '/USE=/{print }' < Makefile)
			echo -e "\t${t#USE*=}"
			echo ")"
		    fi
		} > CMakeLists.txt
		if ! grep -q "!/CMakeLists.txt" .gitignore && [[ -f .gitignore ]]; then
		    echo  "!/CMakeLists.txt" >>.gitignore
		fi
	    )
	fi
    fi
}


for dir in */; do
    createCMakeLists "${dir///}"
done
