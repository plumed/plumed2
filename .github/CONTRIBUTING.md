
You are welcome to contribute to making PLUMED a better tool.

If you want to submit a small change in the code, such as fixing a bug
or adding a small feature, the best way to do it is to open a pull request.
If it is a bug fix, please open a merge request on the oldest maitained branch
where the bug appears. If you are submitting a new feature, please open a merge request on the
latest developer branch (master). Please ensure that your change does not break any of the
existing regression tests on Travis. In addition, add regression tests for any new features that you have added.  If you do not do 
this it is very likely that someone else will break the code for your new feature it in the future.
We (the developers) need to be certain that we will be able to maintain
your code in the future.  Consequently, please  be ready to answer to specific questions about your changes.
Finally, we are very happy to accept contributions to the documentation.

Notice that when you open a pull request you
*implictly agree to transfer the copyright of your code to the PLUMED developers*.
We understand that you might think this unfair.  However, we want to be 100% sure that in the
future we can make drastic changes to the code, including changes to the  license and that we will not have to 
contact all the developers that contributed a small number of lines when doing so.

If you want to contribute some large change, notice that
we prefer not to include your code directly in PLUMED. Large external contributions are difficult for us to maintain.
To be more specific, a contribution can be considered small (and maintainable) only if one
of the current PLUMED developers is able to read every single line of it and understand what it does.
We have, however, put considerable effort into making the code modular.  Consequently, adding a new collective variable is often a single matter of including one additional
file. We would therefore suggest that, one option you might consider is putting your extra files on the web in some way.
In other words, you can just distribute a single cpp file together with instructions on how to include this with PLUMED.
Alternatively, and this is the option we would recommend, you can distribution a separate fork of the full PLUMED with your added code inside.  Creating such forks is relatively straightforward if you use github.
In addition, if you use this second object people will be able to download a fully validated version of plumed with your extra feature,
and will, furthermore, be able to merge the latest changes of the code.  On a final note, we are more than happy to put a link to your webpage/repository/whatever in the PLUMED manual so as to make it is visible to our community of users.
