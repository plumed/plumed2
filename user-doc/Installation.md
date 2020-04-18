\page Installation Installation

This page describes the various ways that you can install PLUMED. Start by selecting the codes that you would like to use PLUMED with and the method that you would like 
to employ in building PLUMED from the drop-down menu below.

@htmlonly
<script>
function getMDCodes() {
 var mdcodedata = ""; 
@MDCODEFUNC@
 var mdcodeother_cbox=document.getElementById("mdcodeother_cbox");
 if( mdcodeother_cbox.checked ) {
   var mdcodeother_data=document.getElementById("mdcodeother");
   mdcodedata += mdcodeother_data.innerHTML;
 }
 return mdcodedata;
}

var current_instructions="unset";

function updateCheck( name ) {
 var mybox = document.getElementById( name + "_cbox");
 if ( mybox.checked ) { mybox.checked=true; }
 else { mybox.checked=false; }
 showInstructions( current_instructions );
}

function showInstructions( name ) {
 current_instructions=name;
 var mydiv = document.getElementById("installdiv");
 if( name=="local" ) {
   var mydata1 = document.getElementById("local-1");
   var mydata2 = document.getElementById("config-opts");
   var mydata3 = document.getElementById("compiling");
   var mydata4 = document.getElementById("installing");
   var mydata5 = document.getElementById("testing");
   var mydata6 = document.getElementById("modules-2");
   mydiv.innerHTML = mydata1.innerHTML + mydata2.innerHTML + mydata3.innerHTML + mydata4.innerHTML + mydata5.innerHTML + mydata6.innerHTML + getMDCodes();
 } else if( name=="cluster" ) { 
   var mydata1 = document.getElementById("cluster-1");
   var mydata2 = document.getElementById("config-opts");
   var mydata3 = document.getElementById("compiling");
   var mydata4 = document.getElementById("installing");
   var mydata5 = document.getElementById("testing");
   var mydata6 = document.getElementById("modules-1"); 
   var mydata7 = document.getElementById("modules-3");
   mydiv.innerHTML = mydata1.innerHTML + mydata2.innerHTML + mydata3.innerHTML + mydata4.innerHTML + mydata5.innerHTML + mydata6.innerHTML + mydata7.innerHTML + getMDCodes(); 
 } else if( name=="multiple" ) {
   var mydata1 = document.getElementById("multiple-1");
   var mydata2 = document.getElementById("config-opts");
   var mydata3 = document.getElementById("compiling");
   var mydata4 = document.getElementById("installing");
   var mydata5 = document.getElementById("testing");
   var mydata6 = document.getElementById("modules-1"); 
   var mydata7 = document.getElementById("modules-3"); 
   mydiv.innerHTML = mydata1.innerHTML + mydata2.innerHTML + mydata3.innerHTML + mydata4.innerHTML + mydata5.innerHTML + mydata6.innerHTML + mydata7.innerHTML + getMDCodes(); 
 } else if( name=="developer" ) {
   var mydata1 = document.getElementById("developer-1");
   var mydata2 = document.getElementById("compiling");
   var mydata4 = document.getElementById("--standalone-executable-content");
   var mydata3 = document.getElementById("developer-2");
   mydiv.innerHTML = mydata1.innerHTML + mydata2.innerHTML + mydata4.innerHTML + mydata3.innerHTML + getMDCodes();
 } else if( name=="cross" ) {
   var mydata1 = document.getElementById("cross-1");
   var mydata2 = document.getElementById("config-opts");
   var mydata3 = document.getElementById("compiling");
   var mydata4 = document.getElementById("installing");
   var mydata5 = document.getElementById("cross-testing");
   mydiv.innerHTML = mydata1.innerHTML + mydata2.innerHTML + mydata3.innerHTML + mydata4.innerHTML + mydata5.innerHTML + getMDCodes(); 
 } else if( name=="python" ) {
   var mydata1 = document.getElementById("python-1");
   var mydata2 = document.getElementById("config-opts");
   var mydata3 = document.getElementById("compiling");
   var mydata4 = document.getElementById("installing");
   var mydata5 = document.getElementById("testing");
   mydiv.innerHTML = mydata1.innerHTML + mydata2.innerHTML + mydata3.innerHTML + mydata4.innerHTML + mydata5.innerHTML + getMDCodes(); 
 } else if( name=="unset" ) {
   mydiv.innerHTML = getMDCodes();
 } else {
   showData( name, "installdiv") + getMDCodes();
 }
}
function showData( name, indiv ) {
 var mydiv = document.getElementById(indiv);
 var mydata = document.getElementById(name);
 mydiv.innerHTML = mydata.innerHTML;
}
function swapConfigure(name) {
 var btn = document.getElementById(name + "_button");
 var mydiv = document.getElementById("conf_" + name); 
 if( btn.textContent=="show defaults" ) { 
   btn.textContent = "hide defaults";
   var dataField = document.getElementById(name + "_long");
   mydiv.innerHTML = dataField.innerHTML;
 } else if( btn.textContent=="hide defaults" ) {
   btn.textContent = "show defaults"; 
   var dataField = document.getElementById(name + "_short");
   mydiv.innerHTML = dataField.innerHTML;
 }
}
function openModal( name ) {
 var mymodal = document.getElementById( name );
 mymodal.style.display = "block";  
}
// When the user clicks anywhere outside of the modal, close it
window.onclick = function(event) {
 var cxxmodal = document.getElementById("CXX")
 if (event.target == cxxmodal) { cxxmodal.style.display = "none"; }
 var ccmodal = document.getElementById("CC");
 if (event.target == ccmodal ) { ccmodal.style.display = "none"; }
 var cxxfmodal = document.getElementById("CXXFLAGS");
 if (event.target == cxxfmodal ) { cxxfmodal.style.display = "none"; }
 var ldfmodal = document.getElementById("LDFLAGS");
 if (event.target == ldfmodal ) { ldfmodal.style.display = "none"; }
 var cppfmodal = document.getElementById("CPPFLAGS");
 if (event.target == cppfmodal) { cppfmodal.style.display = "none"; }
 var libsmodal = document.getElementById("LIBS");
 if (event.target == libsmodal) { libsmodal.style.display = "none"; }
 var slibsmodal = document.getElementById("STATIC_LIBS");
 if (event.target == slibsmodal) { slibsmodal.style.display = "none"; } 
 var dislibfmodal = document.getElementById("--disable-libsearch");
 if (event.target == dislibfmodal ) { dislibfmodal.style.display = "none"; }
 var enmmodal = document.getElementById("--enable-modules");
 if (event.target == enmmodal ) { enmmodal.style.display = "none"; }
 var asmjitmodal = document.getElementById("--enable-asmjit");
 if (event.target == asmjitmodal) { asmjitmodal.style.display = "none"; }
 var prefixmodal = document.getElementById("--prefix");
 if (event.target == prefixmodal) { prefixmodal.style.display = "none"; }
 var debugmodal = document.getElementById("--enable-debug");
 if (event.target == debugmodal) { debugmodal.style.display = "none"; }
 var mpimodal = document.getElementById("--enable-mpi");
 if (event.target == mpimodal) { mpimodal.style.display = "none"; }
 var lapackmodal = document.getElementById("--enable-external-lapack");
 if (event.target == lapackmodal) { lapackmodal.style.display = "none"; }
 var modfilemodal = document.getElementById("--enable-molfile-plugins");
 if (event.target == modfilemodal) { modfilemodal.style.display = "none"; }
 var libpathmodal = document.getElementById("LIBRARY_PATH");
 if (event.target == libpathmodal) { libpathmodal.style.display = "none"; }
 var rpathmodal = document.getElementById("--enable-rpath");
 if (event.target == rpathmodal ) { rpathmodal.style.display = "none"; }
 var standmodal = document.getElementById("--standalone-executable");
 if (event.target == standmodal) { standmodal.style.display = "none"; }
 var pymodal = document.getElementById("PYTHON_BIN");
 if (event.target == pymodal) { pymodal.style.display = "none"; }
 var patchmodal = document.getElementById("patch-options");
 if (event.target == patchmodal ) { patchmodal.style.display = "none"; }
}
</script>
<table> <tr>
<td rowspan="@nrows@" style="padding-right: 20px">
<div class="dropdown">
 <button class="dropbtn">How would you like to build PLUMED?</button>
 <div class="dropdown-content">
  <a onclick='showInstructions("local")'>I want to compile a single version of PLUMED on a local machine.</a>
  <a onclick='showInstructions("cluster")'>I want to compile PLUMED on a cluster so it can be used by several users.</a>
  <a onclick='showInstructions("multiple")'>I would like to install multiple versions of PLUMED on a single machine.</a>
  <a onclick='showInstructions("developer")'>I would like to install PLUMED and work on the development of new features.</a>
  <a onclick='showInstructions("cross")'>I would like to cross compile PLUMED.</a>
  <a onclick='showInstructions("macports")'>I would like to install PLUMED using MacPorts.</a>
  <a onclick='showInstructions("conda")'>I would like to install PLUMED using Conda.</a>
  <a onclick='showInstructions("python")'>I would like to call PLUMED from python.</a>
 </div>
</div>
</td>
@MDCODEBOXES@
</tr></table>
<div style="width: 100%; float:left" id="installdiv"></div>
<div id="CXX" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>CXX Environment Variable</h2>
</div>
<div class="modal-body">
<p>The <code>CXX</code> environmental variable controls the c++ compiler that will be used to build PLUMED. If you would like to compile PLUMED with 
MPI the compiler that is passed to this variable should be an MPI compiler. Consequently, if you work on a machine where <code>CXX</code> is
set to a serial compiler and <code>MPICXX</code> is set to a MPI compiler you should configure PLUMED using the command:</p>
<pre class="fragment">
&gt; ./configure CXX="$MPICXX"
</pre>
<p>If the <code>CXX</code> environmental variable is not set equal to an MPI compiler then MPI is <b>not enabled</b>. This behaviour differs
from what some other configure scripts do. To be clear, when PLUMED is complied variables such as <code>MPICXX</code> are completely ignored.</p>

<p> PLUMED does search for MPI and non-MPI compilers with common names even if the <code>CXX</code> environmental variable is not set.
The search is only a few of the possible compiler names, however. If your compiler is named "g++-mp-4.8" you will need to explicitly specify this 
by using the <code>CXX</code> environmental variable.</p>

<p>PLUMED 2.4 requires a compiler that supports C++11. The following compilers (or later versions) should be sufficient:</p>
<ul>
<li> gcc 4.8.1</li>
<li> clang 3.3</li>
<li> intel 15</li>
</ul>

<p>The <code>./configure</code> script will check whether or not your compiler supports C++11.
Some compilers that do not declare full C++11 support have several C++11 features enabled and can thus be used 
to compile PLUMED (this is the case for the intel 15 compiler for instance).</p>  

<p>If you see a warning about C++11 support during <code>./configure</code> it is then important to check that
PLUMED compiles correctly. The best way to do this to to execute the regtests using the command <code>make regtest</code>.</p> 

<p>We regularly test a number of compilers on <a href="https://github.com/plumed/plumed2/tree/master/.travis">travis-ci</a>. These compilers should thus always compile PLUMED correctly.</p>
</div>
</div>
</div>
<div id="CC" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>CC Environment Variable</h2>
</div>
<div class="modal-body">
<p>The <code>CC</code> environmental variable controls the c compiler that is used in the building of PLUMED. The majority of PLUMED is built using the c++ compiler that is 
specified using the <code>CXX</code> environmental variable. A C compiler must still be specified nonetheless as the wrapper that is used to interface PLUMED with the MD codes is 
written in C.</p>

<p>PLUMED does seach for MPI and non-MPI compilers with a number of standard names even if the <code>CC</code> environmental variable is not set.
Only a few of the possible compiler name are searched, however. If your compiler is named "gcc-mp-4.8" you will need to explicitly specify this 
by using the <code>CC</code> environmental variable.</p>
</div>
</div>
</div>
<div id="CXXFLAGS" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>CXXFLAGS Environment Variable</h2>
</div>
<div class="modal-body">
<p>The <code>CXXFLAGS</code> environmental variable allows you to control the compilation options just as you might expect.</p>
</div>
</div>
</div>
<div id="LDFLAGS" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>LDFLAGS Environment Variable</h2>
</div>
<div class="modal-body">
<p>The main purpose of the <code>configure</code> script is to locate all the libraries that are required to build PLUMED. PLUMED will by default check for these libraries in 
certain standard places. If the libraries are not in these standard places, however, then you will need to provide information the configure script on where to find them by using the
<code>LDFLAGS</code>, <code>CPPFLAGS</code> and <code>LIBS</code> flags. If suitable libraries are not found during the execution of the <code>configure</code> script then these 
features will be disabled. You should thus carefully check the log that is output by the <code>configure</code> command to ensure that all the libraries you require have been found. 
</p>

<p>The <code>LDFLAGS</code> variable tells the <code>configure</code> script the names of directories in which the libraries can be found. If the xdrfile libraries are 
in /opt/local (i.e. where MacPorts puts them) then you would use the command below to tell PLUMED where to find them.</p>
<pre class="fragment">
&gt; ./configure LDFLAGS=-L/opt/local/lib CPPFLAGS=-I/opt/local/include
</pre>
</div>
</div>
</div>
<div id="CPPFLAGS" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>CPPFLAGS Environment Variable</h2>
</div>
<div class="modal-body">
<p>The main purpose of the <code>configure</code> script is to locate all the libraries that are required to build PLUMED. PLUMED will by default check for these libraries in 
certain standard places. If the libraries are not in these standard places, however, then you will need to provide information the configure script on where to find them by using the
<code>LDFLAGS</code>, <code>CPPFLAGS</code> and <code>LIBS</code> flags. If suitable libraries are not found during the execution of the <code>configure</code> script then these 
features will be disabled. You should thus carefully check the log that is output by the <code>configure</code> command to ensure that all the libraries you require have been found.</p>

<p>The <code>CPPFLAGS</code> variable tells the <code>configure</code> script the names of directories in which the include files for the libraries can be found. If the xdrfile libraries are 
in /opt/local (i.e. where MacPorts puts them) then you would use the command below to tell PLUMED where to find them.</p>
<pre class="fragment">
&gt; ./configure LDFLAGS=-L/opt/local/lib CPPFLAGS=-I/opt/local/include
</pre>
</div>
</div>
</div>
<div id="STATIC_LIBS" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>STATIC_LIBS Environment Variable</h2>
</div>
<div class="modal-body">
<p>By adding libararies to the <code>STATIC_LIBS</code> environmental variable you ensure that they are used in all linking steps. Those libraries that are speicified in the <code>LIBS</code>
environmental variable are only used when the PLUMED kernel library is linked.</p>
</div>
</div>
</div>
<div id="LIBS" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>LIBS Environment Variable</h2>
</div>
<div class="modal-body">
<p>The main purpose of the <code>configure</code> script is to locate all the libraries that are required to build PLUMED. PLUMED will by default check for these libraries in 
certain standard places. If the libraries are not in these standard places, however, then you will need to provide information the configure script on where to find them by using the
<code>LDFLAGS</code>, <code>CPPFLAGS</code> and <code>LIBS</code> flags. If suitable libraries are not found during the execution of the <code>configure</code> script then these 
features will be disabled. You should thus carefully check the log that is output by the <code>configure</code> command to ensure that all the libraries you require have been found.</p>

<p>You are required to set the <code>LIBS</code> environmental variable if your libraries have non-standard names. If, for example, your xdrfile library is 
called /opt/local/lib/libmyxdrfile.so you can link it using:</p>
<pre class="fragment">
&gt; ./configure LDFLAGS=-L/opt/local/lib CPPFLAGS=-I/opt/local/include LIBS=-lmyxdrfile
</pre>
<p> This command works because PLUMED first tries to link a routine from the xdr library without adding any additional flag. If and only if this command fails the additional flag, 
"-lxdrfile" is added to the <code>LIBS</code> envrionmental variable. Consequently, if the user has specified the name of the xdr library using the <code>LIBS</code> environmental
flag this custom version will be used in place of the standardly-named "-lxdrfile" library.</p>

<p>If, for any reason, you would like more control over the libraries that are used to build PLUMED you can use the 
<a onclick='openModal("--disable-libsearch")'>--disable-libsearch</a> option.</p>
</div>
</div>
</div>
<div id="--disable-libsearch" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>--disable-libsearch option</h2>
</div>
<div class="modal-body">
<p> By default PLUMED uses a two step procedure when searching for a library. If, for example, PLUMED is searching for the xdr library it will first try to link a routine from this library 
without adding any additional flags. If this command fails the additional flag, "-lxdrfile" is then added to the <code>LIBS</code> envrionmental variable and the test is run again.
The <code>--disable-libsearch</code> command ensures that this second step in the search is ommitted. Consequently, when this option is used only those libraries that are 
explicitly passed through the <code>LIBS</code> environmental variable are linked. For example if the command:</p>
<pre class="fragment">
&gt; ./configure --disable-libsearch LIBS=-lxdrfile
</pre>
<p>is used then only xdrfile is linked. The BLAS and LAPACK libraries are not be linked and the internal versions of these libraries are used. This option is 
useful when installing PLUMED within package managers such as MacPorts as it ensures that only the desired libraries are linked. In other words, no spurious
dependencies are introduced. The only exception to this rule is <code>-ldl</code>, which is a system library on Linux.</p>
</div>
</div>
</div>
<div id="--enable-modules" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>--enable-modules option</h2>
</div>
<div class="modal-body">
<p>PLUMED is made up of modules. Some of these are used more frequently than others and so the ones that we believe that are used 
less often are disabled in a default compilation. If you configure using the option <code>--enable-modules=all</code> then PLUMED will compile 
all these modules. This option would be recommended if you are compiling the code for multiple users and if you are not sure of those 
users need.</p>

<p>If, however, you are compiling the code for yourself only then it may be better to only compile those modules you specifically need.
If you would like to compile with the core modules and the dimred and ves modules you could use the command <code>--enable-modules=dimred:ves</code>.</p>
</div>
</div>
</div>
<div id="--enable-asmjit" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>--enable-asmjit option</h2>
</div>
<div class="modal-body">
<p>If there are a lot of <a class="el" href="_c_u_s_t_o_m.html">CUSTOM</a> functions or <a class="el" href="switchingfunction.html">switching functions</a> 
in your input then you may be heavilty using the lepton library that is included in PLUMED.
The calls to this library can be made significantly faster by using a 
<a href="https://github.com/asmjit/asmjit.git">just-in-time compiler</a>. Furthermore, 
as of PLUMED 2.6, ASMJIT is embedded in PLUMED. To enable it you use the <code>--enable-asmjit</code> flag in configure</p>

<p>If, for any reason, you want to disable ASMJIT at runtime you can do by setting the <code>PLUMED_USE_ASMJIT</code> environment variable as shown below:</p>

<pre class="fragment">export PLUMED_USE_ASMJIT=no</pre>
</div>
</div>
</div>
<div id="--prefix" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>--prefix option</h2>
</div>
<div class="modal-body">
<p>The <code>--prefix</code> environmental variable controls the location in which PLUMED is installed when the <code>make install</code> command is run. By default
PLUMED is installed in <code>/usr/local/</code> and you will thus have to run <code>make install</code> with super-user permissions i.e. using <code>sudo make install</code>.
You can, however, specify the location in which to install PLUMED by using the <code>--prefix</code> environmental variable. For example, to install PLUMED in <code>$HOME/opt</code>
you would configure using:
</p>
<pre class="fragment">
&gt; ./configure --prefix=$HOME/opt
</pre>
<p>Upon install, the executable would then be copied to $HOME/opt/bin, the libraries to $HOME/opt/lib, the include files to $HOME/opt/include, and the
documentation to $HOME/opt/shared/doc/plumed. In addition, a directory called $HOME/opt/lib/plumed will be created containing files such as the
patch files and the object files (for static patches). Advanced users can further customise the install directories using 
using standard autoconf directories (e.g. <code>./configure --bindir=/usr/bin64</code>). Further information on these options can be obtained by running the command:
<pre class="fragment">
&gt; ./configure --help
</pre>
</div>
</div>
</div>
<div id="--enable-debug" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>--enable-debug option</h2>
</div>
<div class="modal-body">
This option is useful if you are implementing new funcationality. If it is enabled then PLUMED 
is compiled with with various debug flags in place that check if code is workign as it should during 
execution. These additional checks will slow down PLUMED, however. In addition, the 
full symbol tables will be written in the executable and the final executable will be much larger.
</div>
</div>
</div>
<div id="--enable-mpi" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>--enable-mpi option</h2>
</div>
<div class="modal-body">
<p>An MPI search is enabled by default and so the <code>--enable-mpi</code> flag is not strictly necessary. It is important to note, however, that
if MPI search is enabled then compilers named "mpic++" and "mpicxx" are searched for first, which may be confusing as many 
other programs behave in a different way.</p>  

<p>When PLUMED is configured with MPI enabled (as it is by default) autoconf checks if a code containing MPI calls can be compiled using the specified 
compiler. If this code can be compiled then MPI is enabled. If it cannot then MPI will be disabled. If you thus specify a non-MPI compiler using the
<code>CXX</code> environmental variable there is no reason to also specify <code>--disable-mpi</code>. The <code>--disable-mpi</code> command is only 
necessary if you specifed an MPI compiler to the <code>CXX</code> environmental variable but you don't want PLUMED to be compiled with MPI support.
In other words, the correct way to enable MPI is to pass <code>./configure</code> an MPI C++ compiler by using the <code>CXX</code> environmental variable.
If this is done you can then treat the MPI library in the same way that all the other libraries that PLUMED tries to link by default.</p>
</div>
</div>
</div>
<div id="--enable-external-lapack" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>--enable-external-lapack option</h2>
</div>
<div class="modal-body">
<p>PLUMED uses the blas and lapack libraries to perform linear algebra and the source code for these two libraries is included. <b>For production calculations,
however, we would strongly encourage you to link to libraries containing suitably-optimized versions of these libraries and to not rely on the versions of blas and 
lapack that are included within PLUMED.</b></p>

<p>You can use a specific version of BLAS or LAPACK availble to <code>configure</code> using the <code>LDFLAGS</code> and <code>LIBS</code> environment variables as shown
below:</p>
<pre class="fragment">
&gt; ./configure LDFLAGS=-L/path/to/blas/lib LIBS=-lnameoflib
</pre>
<p>The <code>configure</code> script will check if the functions that are required from this library are found in the usual way (i.e. by checking without additional 
-llapack and -lblas options added to the <code>LIBS</code> environmental variable first and then with these options added to these variabels). Furthermore, two checks 
for any required functions are performed. In the first, of these checks the name of the functions are searched for. If this search is unsucessful then a 
search for a version of the function with a final underscore added to the name is then performed. If no functions are found in these two steps then the versions of 
blas and lapack that are internal to PLUMED are used.</p>  

<p>If you prefer to disable the search for external versions of blas and lapack because, for example, the system libraries have problems you can configure using the 
command</p>
<pre class="fragment">
&gt; ./configure --disable-external-blas
</pre>
<p>You can also only disable the external LAPACK. In other words, you can use the internal version of LAPACK with an external version of BLAS by using the command:</p>
<pre class="fragment">
&gt; ./configure --disable-external-lapack
</pre>
<p> Using this option is quite sensible as one can typically only heavily optimize the BLAS library. Using the internal LAPACK with an external version of BLAS 
should, therefore, not slow your calculations down significantly. You could thus use this option on systems where the native LAPACK libraries have problems</p> 

<p>We have had numerous emails from users who have struggled to link PLUMED with BLAS and LAPACK. We have even struggled with this ourselves. From these experiences, however,
some sensible steps that you can use to check whether or not the configuration has been setup correctly have emerged. If you are having difficulties you should 
Open the Makefile.conf file that is output after <code>configure</code> has finished running. Check that the flags necessary for loading the 
BLAS and LAPACK libraries are included in the the <code>DYNAMIC_LIB</code> variable that is contained in this file. If BLAS and LAPACK have been detected correctly
-llapack, -lblas and, in some cases, -lgfortran will appear after this variable. On some machines full path specification with -L may be necessary in place of -llapack, 
-lblas anmd -gfortran. Depending on system configuration, your libraries may even not be called -llapack and -lblas. You will thus have to ensure that the correct names 
for these libraries are used. If the correct libraries appear after <code>DYNAMIC_LIB</code> in the Makefile.conf file and if PLUMED still doesn't compile some other things
to try include:
<ul>
<li> If the linker complains and suggests recompiling LAPACK with -fPIC, it means that you have static LAPACK libraries. Either install dynamic LAPACK libraries
or switch to static compilation of PLUMED by unsetting the SOEXT variable in the configuration file.</li>
<li> If the linker complains about other missing functions (typically functions with names that start with the prefix "for_") you need to link some additional Fortran libraries. 
 PLUMED is written in C++ and C++ linkers often do not include Fortran libraries by default, which is unfortunate as the fortran libraries are required for LAPACK and BLAS to 
 work. If you encounter this problem we woudl recommend checking the documentation for your compiler.</li>
<li> If the linker complains that dsyevr_ cannot be found, try adding the flag -DF77_NO_UNDERSCORE to <code>CPPFLAGS</code>. "./configure" will automatically try this solution but maybe you 
 have more luck.</li>
</ul>
</div>
</div>
</div>
<div id="--enable-molfile-plugins" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>--enable-molfile-plugins option</h2>
</div>
<div class="modal-body">
<p>The <code>--enable-molfile-plugins</code> flag ensures that PLUMED is linked with VMD plugins. If PLUMED is configured with this flag you will thus be able to 
read many more trajectory formats for analysis with plumed driver. This option is thus enabled by deafult so including the <code>--enable-molfile-plugins</code> flag  
is not strictly necessary. Furthermore, PLUMED includes the molfile sourcecode and you thus do not need to link to an external library. If, however, you wish to 
use your own version of the mofile plugins you can get them by downloading the SOURCE of VMD. This sourcecode contains
a plugins directory. You will need to adapt the build.sh script within this directory and then compile the code within the directory. Once this process is completed
you should get the molfile plugins compiled as a static library called <code>libmolfile_plugin.a</code>. The location of this library, the <code>libmolfile_plugin.h</code> and 
the <code>molfile_plugin.h</code> include files will then need to be passed to configure using a command similar to:</p>
<pre class="fragment">
&gt; ./configure LDFLAGS="-L/pathtovmdplugins/ARCH/molfile" CPPFLAGS="-I/pathtovmdplugins/include -I/pathtovmdplugins/ARCH/molfile"
</pre>
<p>It may also be necessary to add the location of the TCL interpreter to the <code>LDFLAGS</code> using a command such as:</p>
<pre class="fragment">
&gt; ./configure LDFLAGS="-ltcl8.5 -L/mypathtotcl -L/pathtovmdplugins/ARCH/molfile" \
         CPPFLAGS="-I/pathtovmdplugins/include -I/pathtovmdplugins/ARCH/molfile"
</pre>
<p>To be clear, however, the molfile plugin sourcecode is included in PLUMED and you can thus use these features even if you do not download VMD on your machine.</p>
</div>
</div>
</div>
<div id="LIBRARY_PATH" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>LIBRARY_PATH envrionmental variable</h2>
</div>
<div class="modal-body">
<p> In a typical environment configured using the <a href="http://modules.sourceforge.net">module framework</a> the <code>LIBRARY_PATH</code> environment variable 
contains the path to all the modules loaded at compilation time. When PLUMED is compiled using the <a onclick='openModal("--enable-rpath")'>--enable-rpath</a> 
option the paths defined in <code>LIBRARY_PATH</code> are automatically hard coded into the PLUMED library.</p>
<div style="width: 100%; float:left" id="librpath"></div>
</div>
</div>
</div>
<div id="--enable-rpath" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>--enable-rpath option</h2>
</div>
<div class="modal-body">
<p>On OSX it is common practice to hard code the full path to all the libraries that were used during the compilation within the compiled library. 
The consequence of this is that, if an executable is linked to shared library, then that executable can find all the libraries which are required to run the 
functions within the linked library. The operating system knows where to search for the libraries that were used when the linked libraries was compiled as the 
path to these libraries is stored within the linked library. When the <code>libplumed.dylib</code> library is compiled on a system running OSX there is thus no 
problem in finding external blas or lapack libraries at later times.</p>

<p>On Linux it is not common practice to hard code paths into libraries in the manner described above. Consequently, if you use the 
<code>LDFLAGS</code> option to specify the path to a library at compile time. In other words, if you configure PLUMED using the command:</p>
<pre class="fragment">
&gt; ./configure LDFLAGS="-L/opt/local/lib"
</pre>
<p>Then the libraries in /opt/local/lib may not be found at runtime. One visible symptom of this problem is that <code>src/lib/plumed-shared</code> will 
not be linked correctly. The consequence of this is that it is not possible to dynamically link PLUMED with an MD code later, which is only a problem if you 
would like to link PLUMED dynamically.</p>

<p>You can resolve this issue by using the command below when configuring PLUMED:</p>
<pre class="fragment">
&gt; ./configure LIBRARY_PATH=/path --enable-rpath
</pre>
<p>This command will hard code the paths defined in <code>LIBRARY_PATH</code> into the PLUMED library and is equivalent to the command:</p>
<pre class="fragment">
&gt; ./configure LDFLAGS="-L/path -Wl,-rpath,/path"
</pre>
</div>
</div>
</div>
<div style="display:none;" id="config-opts">
<p>You can control the compilers and compiler flags that are used to build PLUMED by setting the environmental variables. For example, if you wanted to use
the icpc and icc compilers you might use the command shown in the example below:</p>
         
@configure("./configure CXX=icpc CC=icc CXXFLAGS=-O3 LDFLAGS=-L/opt/local/lib CPPFLAGS=-I/opt/local/include LIBS=-lmyxdrfile" copts1)@

<p>You can enable and disable various features within PLUMED by using the --enable-FEATURE and --disable-FEATURE options during the configure step. We would
recommend using a configure command with the following options enabled.</p>

@configure("./configure LIBRARY_PATH=/path --enable-rpath --enable-modules=all --enable-mpi --enable-asmjit --enable-external-lapack --enable-molfile-plugins --prefix=/usr/local" copts2)@
</div>
<div style="display:none;" id="local-1">
<h1> Building PLUMED on a local machine </h1>
<h2>Configuring PLUMED</h2>

<p>The first step in building PLUMED is to configure the makefiles based on the setup of your computer. You can do this by issuing the following command:</p>

@configure("./configure" localinp1)@
</div>
<div style="display:none;" id="multiple-1">
<h1> Building multiple PLUMED versions </h1>
<h2>Configuring PLUMED</h2>

<p>The first step in building PLUMED is to configure the makefiles based on the setup of your computer. If you would like to have multiple versions of the code
installed on your computer, you will also need to give them all different names and place them in different locations. You will thus need to use a configure 
command such as the one below: </p>

@configure("./configure prefix=$HOME/opt --program-suffix=2.2 --program-prefix=mpi-" multiinp1)@

<p>By using this command you ensure that the PLUMED executible will be named <code>mpi-plumed2.2</code> and that other PLUMED files will be named similarly. To load this PLUMED library 
you will thus need to use the flag <code>-lmpi-plumed2.2</code>. PLUMED header files would then be included by using <code>#include <mpi-plumed_2.2/tools/Vector.h></code>.  
It is also possible to use arbitrary scripts to edit the name of the executable by adding the option <code>--program-transform-name=PROGRAM</code> to your configure command.
(see <a href="http://www.gnu.org/software/autoconf/manual/autoconf-2.69/html_node/Transformation-Examples.html#Transformation-Examples"> autoconf documentation </a> for more info).
This options is useful if you do not want to set up modules. As detailed below, however, we believe that using modules is more flexible.</p> 
</div>
<div style="display:none;" id="developer-1">
<h1> Building PLUMED for development purposes </h1>
<h2>Configuring PLUMED</h2>

<p>The first step in building PLUMED is to configure the makefiles based on the setup of your computer. If you are developing PLUMED, we recommend that you configure
using the options below:</p>

@configure("./configure --enable-modules=all --enable-mpi --enable-debug" devinp1)@

<p>Once you are confident that your feature is working correctly you can then reconfigure and compile an optimized version of the code.</p>

<p>You can control the compilers and compiler flags that are used to build PLUMED by setting the environment variables. For example, if you wanted to use
the icpc and icc compilers you might use the command shown in the example below:</p>

@configure("./configure CXX=icpc CC=icc CXXFLAGS=-O3 LDFLAGS=-L/opt/local/lib CPPFLAGS=-I/opt/local/include LIBS=-lmyxdrfile" devinp2)@
</div>
<div style="display:none;" id="developer-2">
<h2> Running PLUMED</h2>

<p>If you are developing PLUMED you are probably compiling the code regularly. You thus might choose to <b>not</b> run the <code>make install</code> command
every time you recompile the code. If you are working in this way and if you use the bash shell you can can run PLUMED from the compilation directory by using the <code>sourceme.sh</code>
file that was created by configure script. This file appears in the main PLUMED directory and can be "sourced" as shown below:</p>
<pre class="fragment">
&gt; source sourceme.sh
</pre>
<p>Running this command ensures that the "plumed" executable appears in your path. To test this executable you might run:</p>
<pre class="fragment">
&gt; plumed -h
</pre>
<h2> Testing PLUMED </h2>
<p>You can test if plumed has been compiled correctly by using the commands:</p>
<pre class="fragment">
&gt; make check
</pre>
<p>Alternatively, you can run all tests using the following command: </p>
<pre class="fragment">
&gt; source sourceme.sh
&gt; cd regtest
&gt; make
</pre>
<p> You can even run a particular test by changing to the directory that contains it and by running the command <code>make</code> within that directory. </p>
<h2> Merging your changes </h2>
<p>Before starting work on a new feature that you plan to share with the PLUMED community we would ask you to read 
<a hred="../../developer-doc/html/_how_to_contribute_to_plumed.html">this brief summary of best practise.</a></p>

<p>If you finish your coding a new feature in PLUMED and if you want to share it with the community then you can open a 
<a href="https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request" >pull request</a> on the PLUMED
repository. Before doing so we would ask that you:</p>
<ul>
<li> Write suitable <a href="../../developer-doc/html/using_doxygen.html">documentation</a> for your new feature</li>
<li> Add a new <a href="../../developer-doc/html/regtests.html">regtest or modified an existing regtest</a> to validate your changes</li> 
<li> <a href="../../developer-doc/html/_using_external_libs.html">Modify the configure script</a> so that any libraries your feature uses can be setup at configure time</li>
<li> Check the format of your code using <a href="../../developer-doc/html/_plumedcheck.html">plumedcheck</a> and 
<a href="../..//developer-doc/html/_code_formatting.html">astyle</a></li>
</ul>
</div>
<div style="display:none;" id="cluster-1">
<h1> Building PLUMED on a cluster </h1>
<h2>Some general advice</h2>
<p>If you are compiling PLUMED on a cluster and if several users will take advantage of the code we would recomment that you:</p>
<ul>
<li> Follow the advice given below about using modules and the module file that PLUMED provides to set up the environment. N.B you will need to edit this file to make it suitable for your environment.</li>
<li> Compile with all the modules enabled so that users can use the full range of features. In addition, ensure that features that are provided by specific libraries are enabled. For example, users employing gromacs will likely want to use libxdrfile to write trr/xtc files. If gromacs is installed before PLUMED, and if the location of libxdrfile is provided to PLUMED during the configure stage then this feature will be available.</li>
<li> Try to patch all MD codes using the <code>--runtime</code> option. When codes are patched in this way users are able to combine any of the installed gromacs/amber/etc versions with any of the 
installed PLUMED versions. It is sometimes claimed that statically linked codes are faster. In our experience, however, this is not true. Furthermore, we have found that non-trivial linking issues are often encountered when building static executables because PLUMED is written in C++ and thus requires the appropriate C++ library to be linked as well as other additional libraries (e.g. libxdrfile).</li>
<li> Keep track of the version of PLUMED you used to patch each of the MD codes. You might do this by giving the MD code modules names such as <code>gromacs/4.6.7p2.2.0</code> to indicate that this is gromacs patched with PLUMED 2.2.0. By keeping track of this information you make it more straightforward for users to report errors on our forum.</li>  
<li> Use the environment variables descsribed below to control the compiler flags and to build an optimized executable.</li>
</ul>
<p>Advice on specific topics and architectures is provided by the drop down menu below. If you have advice on compiling PLUMED on specific architecture please share your tricks! 
You can post information in your blog, or ask as to update the list of topics in this dropdown. The advice you can provide others on compiling PLUMED is very useful and we will happily 
post it here. </p>

@computer-data@

<h2>Configuring PLUMED</h2>
<p>As a bare minimum we would recommend using a command like the one shown below when configuring PLUMED on a cluster</p>

@configure("./configure --enable-modules=all --enable-rpath" clustp1)@

<p>The <code>--enable-rpath</code> flag is recommended as by including this flag you ensure that PLUMED will remember the locations of any runtime libraries that were used at compile time.
The locations of these libraries will thus not need to provided at runtime. This trick often does not work for fundamental libraries such as the C++ and MPI library. The PLUMED module 
should thus load the compiler and MPI modules.</p>
</div>
<div style="display:none;" id="cross-1">
<h1>Cross compiling PLUMED</h1>
<h2>Configuring PLUMED</h2>
<p>If you are cross compiling PLUMED on a different machine than you intend to run it on then you will need to begin by configuring using the usual command.</p>

@configure("./configure" cross1)@

<p>You will likely need to use the options detailed below from controlling the compilers and compiler flags as detailed below, however.</p>
</div>
<div style="display:none;" id="cross-testing">
<h2> Testing PLUMED </h2>
<p>If you have compiled the executable on a machine that is different from the one on which you will ultimately run it the 
<code>plumed</code> executable is not available in the compilation environment. You thus cannot perform the regtests on
the machine and you cannot compile the manual.</p>

<p>You might try to run the regtests on the computing nodes in this case. You may need to do some tweaks to get this 
to work as machines where people do cross-compiling often have architectures with limited capabilities on the compute nodes.</p> 

<h2> Patching MD codes with a cross compiled code </h2>
<p>As discussed in the previous section, if PLUMED has been cross compiled the <code>plumed</code> executable is not available in the compilation environment.
You thus cannot run <code>plumed patch</code> on the compiler nodes. It is also likely that it will not be possible to run this command on the compute nodes
as on machines where cross compilation is necessary it is often not possible to fork new processes from the compute nodes. To get around these issues you can 
use the command: </p>
<pre class="fragment">
&gt; plumed-patch 
</pre>
<p>instead. This script provides a "shell only" implementation of <code>plumed patch</code> and thus does not launch of the <code>plumed</code> executable.
You can thus run this command (and patch MD codes) on the compiler nodes. If you have installed PLUMED you can find the <code>plumed-patch</code> code
script (and some other similar shell scripts that allow you to run various plumed command line tools without launching the <code>plumed</code> executable)
in the directory <code>$prefix/plumed/plumed-</code>. The reason these files are not included in the execution path (prefix/bin) is to avoid clashes.</p>
</div>
<div style="display:none;" id="compiling">
<h2>Compiling PLUMED</h2>
<p>PLUMED can be compiled using the following command Once the <code>configure</code> script has finished running:</p>
<pre class="fragment">
&gt; make -j 4
</pre>
<p>This command compiles the entire code and produces a number of files in the <code>src/lib</code> directory, including the executable
<code>src/lib/plumed</code>. If the shared libraries are enabled a shared libraries called <code>src/lib/libKernel.so</code> should also be generated.
The extension for this file will be <code>.dylib</code> if you are building PLUMED on a Mac. There are some rare occasions when this command fails
because the configure script fails to find some of the required libraries. If this happens you can <a onclick='showData("makefileconf","makeconfdiv")'>edit the Makefile.conf file</a> so that the suitable
compilation options are set up.</p>
<div style="width: 100%; float:left" id="makeconfdiv"></div>
</div>
<div style="display:none;" id="makefileconf">
<p>The most important variables to edit in the Makefile.conf file are:</p>
<ul>
<li> DYNAMIC_LIB : the libraries (e.g. BLAS and LAPACK) that are linked when the the PLUMED library is compiled are listed here.  
Flags such as -L/path/to/xdrfile -lxdrfile that were specified using the environment variables <code>LDFLAGS</code> and <code>LIBS</code> 
in the call to <code>configure</code> will appear here. The fact that the flags specified using <code>LIBS</code> is perhaps confusing but this 
is necessary to keep the configuration files compatible with PLUMED 2.0. Notice that for the PLUMED shared library to be compiled correctly only dynamic libraries should be listed here.</li>
<li> LIBS : the libraries that are needed to patching the MD codes (typically only <code>-ldl</code> which is needed to have functions for dynamic loading) are specified here.</li> 
<li> CPPFLAGS : definitions that are needed to enable specific optional functions (e.g. -D__PLUMED_HAS_XDRFILE to enable the xdrfile library) are specified here.</li>
<li> SOEXT : the extension for shared libraries on your system is specified here. Typically an "so" extension is used on UNIX and a "dylib" extension is used on mac. 
If your system does not support dynamic libraries or if you would like a static executable you can just set this variable blank i.e. <code>SOEXT=</code>. </li>
</ul>
</div>
<div style="display:none;" id="installing">
<h2> Installing PLUMED </h2>
<p>You can install PLUMED using the command:</p>
<pre class="fragment">
&gt; make install
</pre>
<p>or the command:</p>
<pre class="fragment">
&gt; sudo make install
</pre>
<p>If you want to install PLUMED in a system directory. Unless modifications are made to the standard autoconf directories this command copies the 
executable to $prefix/bin, the libraries to $prefix/lib,
the include files to $prefix/include, and the documentation to $prefix/shared/doc/plumed. A directory called
$prefix/lib/plumed is also created by this command. This directory contains several other files, including
the patch files and the object files that are used for static patching.  
<p>$prefix here is the directory specified using the <a onclick='openModal("--prefix")'>--prefix</a> keyword of the configure script.</p>
<p>Once PLUMED has been installed using the <code>make install</code> command you can delete the original compilation directory 
or you can recompile a different PLUMED version in the same place. You should not delete any of the installed files, however, as 
<a onclick='openModal("--standalone-executable")'>PLUMED will not run</a> if there are files missing from these directories</p>
</div>
<div id="--standalone-executable" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>--standalone-executable option</h2>
</div>
<div class="modal-body" id="--standalone-executable-content">
<p>The PLUMED executable which relies on the resource files present in the compilation directory.
During installation these files are copied to $prefix/lib/plumed and the compilation directory can therefore be deleted.
If you do not install PLUMED, however, you need to ensure that none of the files in <code>src/lib/plumed</code> are not
moved or renamed.</p>
<p>The path to the PLUMED root directory is hard coded in the plumed executable as can be verified by using the command:</p>
<pre class="fragment">
&gt; plumed info --root
</pre>
<p>If you try to run a non-installed plumed executable and if the <code>src/lib/plumed</code>  
directory is not in place PLUMED will throw an error as shown below:</p>
<pre class="fragment">
&gt; plumed help
ERROR: I cannot find /xxx/yyy/patches directory
</pre>
<p>You can force plumed to ignore this error and to run anyway by using the option --standalone-executable as shown below</p>
<pre class="fragment">
&gt; plumed --standalone-executable help
</pre>
<p>Many features are not be available if you run in this way. This is currently the only way to use the PLUMED static executable on Windows, however.</p>
</div>
</div>
</div>
<div style="display:none;" id="testing">
<h2> Testing PLUMED </h2>
<p>It is important to test PLUMED every time you install PLUMED, as Even though we regularly perform tests on 
<a href="http://travis-ci.org/plumed/plumed2">[Travis-CI]</a>, it is possible that aggressive optimization or even architecture dependent features
trigger bugs that do not show up on travis. To test the PLUMED executable that you have installed you should, therefore, type</p>
<pre class="fragment">&gt; make installcheck </pre>
<p>Running this command tells PLUMED to run all the tests in the test suite using the version of PLUMED that was just installed. If you would like to test
PLUMED using some other PLUMED version then you can. The tests can be run using the PLUMED executable that is in your current path by using the commands: </p>
<pre class="fragment">
&gt; cd regtest
&gt; make
</pre>
<p>To determine the version of PLUMED that has been run using the commands above you can use the command:</p>
<pre class="fragment">
&gt; which plumed
</pre>
<p>When you do this you might get certain <a onclick='showData("testerrors","errdiv")'>innocous errors</a></p>
<div style="width: 100%; float:left" id="errdiv"></div>
</div>
<div style="display:none;" id="testerrors">
<p>If you see errors when you run the regression tests using the commands:</p>
<pre class="fragment">
&gt; cd regtest
&gt; make
</pre>
<p>The first thing you should try to work out is the version of PLUMED that is being used to run the tests. When the tests are run 
using the commands above the version of PLUMED that is found in the PATH is used. If the version of PLUMED in your path is not the same
as the version that you just downloaded and compiled then you may observe one of two innocuous errors:</p>
<ul>
<li> If the PLUMED executable in your path is older than the test suite you are running on the tests might fail. In this case the tests
fail as some feature was introduced to PLUMED in a newer version.</li>
<li> The tests might also fail if the PLUMED executable in your path is newer than the test suite. In this case some non-backward compatible change 
 was made in PLUMED. We try to keep the number of non-backward compatible changes small, but as you can see from the <a class="el" href="_change_log.html">ChangeLog</a> there
 are typically a few non-compatable changes at every new major release.</li>
</ul>
</div>
<div style="display:none;" id="modules-1">
<h2>Setting up your environment</h2>
<p>Once the installation and testing has been completed you are ready to use PLUMED. Running PLUMED will be more straightforward if you set up the environment correctly.
If the environment is setup correctly you should be able to:</p>
<ul>
<li> use the <code>plumed</code> executable from the command line. Notice that you can also use this executable before installing. </li>
<li> link against the PLUMED library using the <code>-lplumed</code> flag for the linker. If this flag can be employed then it is possible to use PLUMED library in general purpose programs </li>
<li> use PLUMED internal functionality (C++ classes) including header files such as <code>#include <plumed/tools/Vector.h></code> in general purpose programs </li>
</ul>
</div>
<div style="display:none;" id="modules-2">
<p>The easiest way to setup the environment is to use <a onclick='showData("modules-3","localmoddiv")'>the module framework</a>. If you have installed plumed into a system directory such 
as <code>/usr/local</code> then the environment should already be setup correctly.</p>
<div style="width: 100%; float:left" id="localmoddiv"></div>
</div>
<div style="display:none;" id="modules-3">
<p>The easiest way to setup the environment is to use <a href="http://modules.sourceforge.net">the module framework</a>.  
A suitable module file for PLUMED can be found in $prefix/lib/plumed/src/lib/modulefile after installation is completed.
You can edit this file or just put it into your modulefile directory directly. If you do so it is then straightforward to 
set up the environment. </p>
<p> Notice that if you have installed more than one version of PLUMED on your machine you can use the module framework to easily
switch between them. </p> 
<p> Lastly note that even if you do not want to use modules it may still be useful to look at the modulefile as this file will tell 
you which environment variables need to be set for PLUMED to work correctly. </p>
</div>
<div style="display:none;" id="macports">
<h1> Building PLUMED with macports</h1>
<p>If you are using PLUMED on a Mac, you can intall it using MacPorts. To take advantage of this option
you proceed as follows:</p>
<ul>
<li> Install <a href="https://www.macports.org/">MacPorts</a></li> 
<li> Type <code>sudo port install plumed</code></li>
</ul>
<p>The default variant that is installed using <code>sudo port install plumed</code> is shipped as a compiled
binary and is thus significantly faster to install. A number of different PLUMED versions can be installed using MacPorts, however.  
You can get a list of the versions that are available by using the command:</p>

<pre class="fragment">&gt; sudo port info plumed</pre>

<p>The various options that can be viewed using this command allow you to install PLUMED using multiple different compilers.
You can thus install plumed with mpich using the command:</p>

<pre class="fragment">&gt; sudo port install plumed +mpich</pre>

<p>We would recommend using the recent clang compiler instead of native compilers to take advantage of openMP. To install
using this option you can use the command:</p> 

<pre class="fragment">&gt; sudo port install plumed +mpich +clang50</pre>

<p>Notice that support for c++11 with gcc compilers in MacPorts is difficult as it is not possible to 
use the system's c++ library. For this reason, we, therefore, only support the clang compilers. If you are 
interested reading more about this issue we would direct you to 
<a href="https://github.com/macports/macports-ports/pull/1252">this discussion</a> </p> 

<p>In addition to the variants that allow you to use MacPorts with a range of compilers there are also options that 
allow you to compile with debug flags (<code>+debug</code>), to pick a linear algebra library
(e.g. <code>+openblas</code>) and to enable all the optional modules (<code>+allmodules</code>).
You can also install a developer version of PLUMED by using the command: </p> 

<pre class="fragment">&gt; sudo port install plumed-devel</pre>

<p>This developer version is typically a later version of PLUMED that is not yet considered stable. When installing
<code>plumed-devel</code> using MacPorts you can use all the variants that were available for <code>plumed</code> to 
customize the compilation. You cannot, however, install <code>plumed-devel</code> and <code>plumed</code> at the same time.</p>

<p>You can also use MacPorts to install a plumed-patched version of gromacs. To install gromacs patched with a stable version of PLUMED that
was compiled using the clang-5.0 compiler and mpich you would use the following sequence of commands:</p>

<pre class="fragment"> 
&gt; sudo port install plumed +mpich +clang50
&gt; sudo port install gromacs-plumed +mpich +clang50
</pre>

<p>If, by contract, the objective was to use a version of gromacs patched with the development version of PLUMED you would use the commands:</p>

<pre class="fragment">
&gt; sudo port install plumed-devel +mpich +clang50
&gt; sudo port install gromacs-plumed +mpich +clang50
</pre>

<p>Notice that the same compiler variants must be used for gromacs and PLUMED (in this example <code>+mpich +clang50</code>). If gromacs and PLUMED are built using different 
compiler variants then compilation will fail.</p>

<p>The patched version of gromacs that can be installed using MacPorts links PLUMED in runtime mode. The path for libplumedKernel.dylib in the MacPorts tree
is hardcoded, however. As a consequence:<p>

<ul>
<li>If gromacs is run and if the <code>PLUMED_KERNEL</code> environment variable is unset (or set to empty), then the version of PLUMED that was installed using MacPorts is used.</li>
<li>If gromacs is run and if the <code>PLUMED_KERNEL</code> environment variable points to another instance of the PLUMED library then this other instance of PLUMED is used in place of 
the version that was installed using MacPorts</li>
</ul>

<p>Having the <code>PLUMED_KERNEL</code> operate in this way is useful if you are developing PLUMED as you can install gromacs once using MacPorts and then combine it with any version of 
PLUMED on your computer.</p>
</div>
<div style="display:none;" id="conda">
<h1><a class="anchor" id="Installation-conda"></a>Installing PLUMED with conda</h1>
<p>If you use the conda package manager you can quictkly install a pre-compiled PLUMED binary by using the following command: </p>
<pre class="fragment">&gt; conda install -c conda-forge plumed </pre>
<p> Similarly, the python wrappers can be installed by using the command: </p>
<pre class="fragment">&gt; conda install -c conda-forge py-plumed </pre>
<p>These packages are part of <a href="https://anaconda.org/conda-forge">conda-forge</a>. They should thus be binary compatible with other codes from 
the same distribution. Furthermore, it should also be possible to combine the PLUMED kernel installed from conda with an MD code that has been compiled outside 
of conda (or within a different conda environment) if PLUMED is linked in runtime mode. The only variable that you need to set in order to access to the PLUMED kernel 
that is installed through conda is <code>PLUMED_KERNEL</code> (e.g., <code>export PLUMED_KERNEL=/conda/prefix/lib/libplumedKernel.so</code>).</p>

<p>Conda binaries are only available for Linux and macOS. The installed conda binaries were configured using the command:</p>

@configure-conda@

<p> so the features are limited accordindly. Notice that there additional conda packages are available on the <a href="https://anaconda.org/plumed/plumed">plumed</a> channel. 
These packages are for testing only.</p>
</div>
<div style="display:none;" id="python-1">
<h1>Calling PLUMED from python</h1>
<p>It is possible to call PLUMED 2.5 and later versions in a python script. If you would like to do so you simply add the command below in your python script:</p>
<div class="fragment"><div class="line"><a name="l00001"></a><span class="lineno">  1</span>&#160;<span class="keyword">import</span> plumed</div></div>
<p>The interface is designed with developers in mind rather than users as the interface from python to PLUMED is similar to the interface that is used to link PLUMED 
with the various MD codes. If you are only compiling one version of PLUMED you might configure using the command: </p>

@configure("./configure PYTHON_BIN=python3.6 --prefix=$HOME/opt" python1)@

<p>Once you have compiled and installed PLUMED by following the instructions that are laid out in the sections that follow, the python wrappers for PLUMED will be in 
<code>$HOME/opt/plumed/python</code>. To use plumed in a python script the location of these wrappers must be added <code>PYTHONPATH</code> environment variable using
a command such as the one shown below:</p>
<pre class="fragment">
&gt; PYTHONPATH="$HOME/opt/lib/plumed/python:$PYTHONPATH"
</pre>
<p>If you have installed muliple PLUMED versions you might find it easier to install the Python wrappers and PLUMED separately. The simplest way to install the python
wrappers by themselves is with <code>pip</code> as shown below:</p>
<pre class="fragment">
&gt; pip3.6 install --user plumed
</pre>
<p> In the command above the <code>--user</code> flag allows you to install the packages on your home. Notice also that you don't even need to download 
PLUMED in order to install the wrappers. The command above can be run before you configure, compile and install PLUMED. You will, however, need PLUMED 
in order to use it within a python script. You can tell the wrappers where PLUMED is by setting the <code>PLUMED_KERNEL</code> environment variable as shown
below:</p>
<pre class="fragment">
&gt; PLUMED_KERNEL=$HOME/opt/lib/libplumedKernel.so
</pre>
<p>When you install the wrappers in this manner described above you will download those that are packaged on <a href="https://pypi.org/project/plumed/">[Pypi]</a>.
If you want to install the development version of the wrappers using pip you should download the PLUMED repository and use the following commands:</p>
<pre class="fragment">
&gt; pip3.6 install --user cython 
&gt; cd plumed2/python
&gt; make pip
&gt; pip3.6 install --user .
</pre>
<p>You are highly recommended to use a virtualenv when installing the development version as you will then ensure that the code you install does not interfere with the released pacakages.</p>
<h2>Customising your complation environment</h2>
</div>
<div style="display:none;" id="mdcodeother">
<h2> Linking other codes with PLUMED </h2> 
<p> 
If your MD code is not listed among those that already have support, you can implement an interface for it by following the instructions in 
the <a href="../../developer-doc/html/index.html">developer manual</a>. 
</p>
</div>
<div id="PYTHON_BIN" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>PYTHON_BIN Environment Variable</h2>
</div>
<div class="modal-body">
<p>The <code>PYTHON_BIN</code> environment variable tells configure which version of python should be used to construct the PLUMED wrappers. Oftetimes this flag is unecessary
as configure will find any <code>python</code> executable that also has the <code>cython</code> module available automatically. If your python interpreter has a name 
that is different from python you will have to set this environment variable.</p>
</div>
</div>
</div>
<div id="patch-options" class="modal">
<div class="modal-content">
<div class="modal-header">
<h2>Patching options</h2>
</div> 
<div class="modal-body">
<p>The way PLUMED is linked is controlled by a flag that is given to the <code>plumed patch</code> command. This flag can be set to any one of the following three options:</p>
<ul>
<li> <b>--static</b> With this flag PLUMED is linked as a collection of object files. This is only suggested if you
 absolutely need a static executable. Notice that when this setting is used it is often more difficult to configure
 the MD code properly as all the libraries that PLUMED depends on need to be specified properly. The <code>./configure</code> script
 does its best to look after all this for you but it cannot solve all the problems you might encounter. For example, we have had reports that 
 this patching mode does not work properly on OSX. </li>
<li> <b>--shared</b> This is the default mode for linking PLUMED. When this option is specified PLUMED is linked as a shared library. 
 One consequence of this is that when PLUMED is updated, there is no need to recompile the MD code.
 Linking with the --shared option is superior to linking with --static as the libraries that PLUMED depends are linked automatically.
 If, however, you later remove the directory where the version of PLUMED that you linked with then then MD code will not run anymore. </li>
<li> <b>--runtime</b> If you patch your MD code using this option then you can choose the location of the PLUMED library at runtime by setting the environment variable PLUMED_KERNEL.
 This option is probably the most flexible of the three, and we encourage system administrators to use this option when installing
 PLUMED on shared facilities. When this setting is used it is possible to update the PLUMED library and the MD code separately. Users can, thereore, combine the MD code with different
 versions of PLUMED at will. If you are using this option we would also recommend using a modulefile to set the runtime environment.</li>
</ul>
<p>Notice that the precise behaviour of the <code>--static</code> flag depends on the PLUMED version. For versions of PLUMED 
earlier than 2.5 there was no possibility to link PLUMED as a static library. In PLUMED 2.5 onwards, however, the 
<code>./configure</code> script tries to set up the system so that a <code>libplumed.a</code> file is produced. If an MD
code is patched with PLUMED 2.5 or later and the <code>--static</code> option the MD code is linked against this static library.
If you wish to revert to the pre-version-2.5 behaviour and to just link with the object files you will need to configure PLUMED using the 
command:
<pre class="fragment">
./configure --disable-static-archive
</pre>
</div>
</div>
</div>
@MDCODEDATA@
@endhtmlonly
