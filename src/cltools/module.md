This module contains the set of simple command line tools that you can use with PLUMED. You can use
one of these tools even when you do not have a molecular dynamics code patched with PLUMED.  

To use one of these tools you issue a command something like:

````
plumed <toolname> <list of input flags for that tool>
````

There are tools that allow you to postprocess trajectories, benchmark your PLUMED calculations and 
to run molecular dynamics on a system of Lennard Jones particles. 

# Using bash autocompletion

When possible, PLUMED tries to install bash autocompletion so that
you can just use the `<TAB>` key to complete
plumed commands (e.g. `plumed dr<TAB>`) or even options (e.g. `plumed driver --i<TAB>`).
In case this does not work, you might have to add the following lines to your .bashrc file:

````
_plumed() { eval "$(plumed --no-mpi completion 2>/dev/null)";}
complete -F _plumed -o default plumed
````

## Effect

When typing on the the shell you should observe the following behavior.

````
> plumed <TAB>
````

will autocomplete with the names of the available PLUMED commands (e.g. `driver`, `help`, etc).

````
> plumed -<TAB>
````

will autocomplete with the available PLUMED options (e.g. `--no-mpi`, etc).

PLUMED also knows which are the options available for each command
(e.g. `plumed driver --natoms`). So, the following

````
> plumed driver -<TAB>
````

(notice the `-`) will autocomplete to the options of `plumed driver`. On the contrary

````
> plumed driver --ixtc <TAB>
````

(notice the there is no `-` before `<TAB>`) will autocomplete to the files in the current directory.

Also notice that every time you use the `<TAB>` key to autocomplete the command `plumed` will be invoked.
This should allow the correct commands and options to be reported depending on the exact `plumed` command
in the current execution path. For instance, if you have multiple PLUMED versions installed with
env modules, you should be able to see the commands available in the currently loaded version.
Clearly, this feature will only be available if `plumed` can run on this machine (that is: will not work
if you are cross compiling). This is not a problem since you are not expecting to run the `plumed` command
in this specific case.

## Technicalities

At configure time if the variable `BASH_COMPLETION_DIR` is defined it will be used to decide where PLUMED
autocompletion should be installed. Otherwise, configure will look for the presence of the `bash-completion` package
and, in case it is installed on the same prefix as PLUMED, also PLUMED autocompletion will be installed.
Finally, if none of these two conditions are satisfied, autocompletion will not be enabled. You will
have to change your bashrc file once adding the following lines:

````
_plumed() { eval "$(plumed --no-mpi completion 2>/dev/null)";}
complete -F _plumed -o default plumed
````

The command `plumed completion` just writes on its standard output the body of a bash function that is
then used by bash to construct the autocompletion.
The `--no-mpi` flag makes it more likely that the command can be executed correctly e.g. when you are on the login node of a cluster and
PLUMED was compiled with MPI but the login node does not support MPI. In other cases, it is harmless.
The `-o default` options will make sure that if `plumed --no-mpi completion` returns an error the default bash completion
will be used. This is what will happen if you load an older PLUMED version for which the `completion` command is not available yet.
In future PLUMED versions the `plumed completion` command might return more sophisticated functions. You should
be able to benefit of these features without ever changing your bash configuration file again.

## Multiple versions and suffixes

If you have multiple versions of PLUMED installed in separate env modules there is nothing more to do.
However, if you have have multiple versions of PLUMED installed with different suffixes you should
consistently add more lines to your profile file. For instance, if you installed two executables named
`plumed` and `plumed_mpi` your configuration file should look like:

````
_plumed() { eval "$(plumed --no-mpi completion 2>/dev/null)";}
complete -F _plumed -o default plumed
_plumed_mpi() { eval "$(plumed_mpi --no-mpi completion 2>/dev/null)";}
complete -F _plumed_mpi -o default plumed_mpi
````
