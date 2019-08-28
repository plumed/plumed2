
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
*implictly agree to transfer the copyright of your code to the PLUMED developers* (or to the authors
of the code that you are modifying).
We understand that you might think this unfair.  However, we want to be 100% sure that in the
future we can make drastic changes to the code, including changes to the  license and that we will not have to 
contact all the developers that contributed a small number of lines when doing so.

If you want to contribute some large change, please consider adding a new module.
Documentation about adding new modules is still limited, but you can get inspiration
from the existing ones. This procedure will allow you to keep the ownership on your code.
On the other hand, we expect that you will maintain it in the future.
In order to incorporate a new module into the main repository, we ask contributors to declare that
the module is available with an open source license.

Finally, notice that you can always share modified versions of PLUMED with your changes.
We are happy if you want to host on github a fork of PLUMED with additional features.
