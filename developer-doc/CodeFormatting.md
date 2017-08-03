\page CodeFormatting How to format code properly

Since version 2.3.2, we format code using <a href="http://astyle.sourceforge.net/"> astyle </a>.
As a convention, we use `astyle` version 3.00, with the options that are
reported in the file `.astyle.options` located in the root directory of PLUMED.
You might want to automatize the application of `astyle` using those options.
As of now, you can use the following command from root PLUMED directory:
\verbatim
> make astyle
> git commit
\endverbatim
Notice that this command will both apply `astyle` to all the relevant files as well as
add them to the next git commit. After having inspected the changes, you can commit them.
Also notice that running this command from PLUMED root directory will also
compile the `astyle` version that is distributed with PLUMED. We decided to distribute
`astyle` within the PLUMED repository to make sure that everyone is using exactly the same version.
In addition, you can run `make astyle` directly within a module directory
so as to only reformat that specific module.

Additional care must be used while merging branches. In this case, you should
make sure that both branches are formatted with `astyle` before merging them.
The procedure discussed below should be made once for each not-yet-formatted branch that you are maintaining.

Let's say that you are working on a branch `feature` that is not yet formatted, and you
would like to merge changes from branch `v2.3` that is already formatted.
You should use the following list of commands
\verbatim
# Bring in all the changes in v2.3 up to astyle formatting (excluded):
> git merge astyle-v2.3~1
# Notice that this will include the astyle scripts, but not
# the big commit formatting the whole code.

# Mark the following commit (astyle-v2.3) as merged, even though
# it is completely ignored:
> git merge -s ours astyle-v2.3
# This is necessary since this commit is too big to be really merged.

# Indeed, instead of merging it, we apply directly astyle to the
# current branch:
> make astyle

# Now the two branches are both formatted and can be merged one into the other.
# Merge all the newer commits on v2.3
> git merge v2.3
\endverbatim

Notice that here `astyle-v2.3` is a tag that refers to the commit where we introduced
formatting in v2.3 branch. 
In a similar way you can bring in any changes from another branch, just replace
`v2.3` with the proper name (e.g. `master`).
After merging, your branch will be also formatted correctly.
Notice that you cannot work it in the opposite direction (that is, unformat an already
formatted branch). Finally, consider that rebasing might be more complicated.
When all branches will be formatted this will not be an issue anymore.

