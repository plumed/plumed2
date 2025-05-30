import os
from sys import platform as sys_platform
import json
import shutil
import requests
import subprocess
import yaml
import glob
import numpy as np
from pathlib import Path
from datetime import date 
from bs4 import BeautifulSoup
from contextlib import contextmanager
from PlumedToHTML import processMarkdown, processMarkdownString, get_javascript, get_css 
import networkx as nx
from multiprocessing import Pool, cpu_count
import functools



PLUMED="plumed"

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def getActionDocumentationFromPlumed(syntax) :
    docs = {}
    tags = {}
    for file in glob.glob("../src/*/*.cpp") : 
        with open(file,"r") as f :
           content = f.read()
        if "//+PLUMEDOC" not in content : continue
        actioname, founddocs,  indocs = "", "", False
        for line in content.splitlines() :
            if "//+ENDPLUMEDOC" in line :
               if not indocs : raise Exception("Found ENDPLUMEDDOC before PLUMEDOC in " + file)
               if len(actioname)>0 : 
                  print("Found documentation for ", actioname)
                  docs[actioname] = founddocs
               actioname, founddocs, indocs = "", "", False
            elif "//+PLUMEDOC" in line :
               if indocs : raise Exception("Found PLUMEDDOC before ENDPLUMEDOC in " + file) 
               indocs, taglist = True, ""
               for key in line.split() :
                   if key=="//+PLUMEDOC" : continue 
                   if key in syntax.keys() or key in plumed_syntax["cltools"].keys() :
                      if len(actioname)>0 : raise Exception("found more than one action name in input for " + key )
                      actioname = key
                   else : 
                      taglist = key + " "
               if actioname in syntax.keys() :
                  tags[actioname] = taglist
            elif indocs and not "/*" in line and not "*/" in line : founddocs += line + "\n"
    return docs, tags

def get_reference(doi):
    # initialize strings
    ref=""
    ref_url=""
    # retrieve citation from doi
    if(len(doi)>0):
      # check if unpublished/submitted
      if(doi.lower()=="unpublished" or doi.lower()=="submitted"):
        ref=doi.lower()
      # get info from doi
      else:
        # get citation
        cit = subprocess.check_output('curl -LH "Accept: text/bibliography; style=science" \'http://dx.doi.org/'+doi+'\'', shell=True).decode('utf-8').strip()
        if("DOI Not Found" in cit):
         ref="DOI not found"
        else:
         # get citation
         ref=cit[3:cit.find(", doi")]
         # and url
         ref_url="http://dx.doi.org/"+doi
    return ref,ref_url

def create_map( URL ) :
    page = requests.get(URL)
    soup = BeautifulSoup(page.content, "html.parser")
    script = soup.find("script",id="actionChart")
    xdata = {}
    ydata = {}
    for line in script.text.splitlines() :
        if "var" not in line or "for" in line :
           continue
        if "xValues" in line and "=" in line :
           xdata = json.loads( line.split("=")[1].replace(";","") )
        if "yValues" in line and "=" in line :
           ydata = json.loads( line.split("=")[1].replace(";","") )
    
    return dict(map(lambda i,j : (i,j) , xdata,ydata))

def printDataTable( f, titles, tabledata, tagdictionary={}, extraAttributes={} ) :
    if "Tags" in titles : 
        if titles[2]!="Tags" : raise Exception("found tag column in surprising location")
        # Get all tags in table
        mytags = set() 
        for line in tabledata :
            for tt in line[2].split() :
               mytags.add(tt)

        f.write("<div class=\"dropdown\">\n")
        f.write("<button class=\"dropbtn\" onclick=\'clearTagSelection()\'>Tags</button>\n")
        f.write("<div class=\"dropdown-content\">\n")
        for tag in mytags :
           f.write("<a href=\"#\" onclick=\'displayActionWithTags(\"" + tag + "\")\' onmouseover=\'displayTagData(\"" + tag + "\")\' onmouseleave=\'displayTagData(\"no selection\")\'>" + tag + "</a>\n")
        f.write("</div>\n")
        f.write("<span id=\"tagdisplay\"></span>\n")
        f.write("</div>\n")
        for tag in mytags : 
            if tag in tagdictionary.keys() :
               f.write("<span id=\"" + tag + "\" style=\"display:none;\"><b>" + tag + ":</b> " + tagdictionary[tag] + "</span>\n")
 
    f.write("<table id=\"browse-table\" class=\"display\">\n")
    f.write("<thead><tr>\n")
    for t in titles :
       usestyle = "text-align: left"
       myattributes = ""
       if t in extraAttributes.keys() :
          attributes=extraAttributes[t]
          if "style" in attributes:
            usestyle += f";{attributes['style']}"
          for attr in attributes.keys() :
            myattributes = f'{attr}="{attributes[attr]}" {myattributes}'

       myattributes = f'style="{usestyle}" {myattributes}'

       f.write(f'<th {myattributes}>{t}</th>\n')
    f.write("</tr></thead><tbody>\n")
    for r in tabledata : 
        if len(r)!=len(titles) :
           raise Exception("mismatch between number of columns in tabledata and number of titles")
        f.write("<tr>\n")
        for e in r :
           f.write("<td style=\"text-align: left\">" + e + "</td>\n")
        f.write("</tr>\n")
    f.write("</tbody></table>\n")

def createSummaryPage( broken, nodocs, undocumented, noexamples, unexempled_keywords ) :
    with open("docs/summary.md","w+") as f : 
       f.write("# Summary of manual coverage\n\n Data on this page was generated on " + date.today().strftime('%B %d, %Y') + "\n")

       if len(broken)>0 : 
          text = f"""
## List of pages with broken examples

There are {len(broken)} action pages with failing inputs.
"""
          f.write( text )
          printDataTable( f, ["Page","# broken"], broken )

       if len(nodocs)>0 :
          ## List of pages with broken examples
          text = f"""
## List of pages with no documentation

There are {len(nodocs)} pages with no in code documentation.
"""
          f.write( text )
          printDataTable( f, ["Page","Type"], nodocs )
     
       if len(undocumented)>0 : 
          text=f"""
## List of actions that have undocumented keywords

There are {len(undocumented)} actions with undocumented keywords
"""
          f.write( text )
          printDataTable( f, ["Action","Module"], undocumented )

       if len(noexamples)>0 :
          text=f"""
## List of actions that have no examples in manual

There are {len(noexamples)} action pages with no examples.
"""
          f.write( text )
          printDataTable( f, ["Action","Module"], noexamples )   

       if len(unexempled_keywords)>0 : 
          text=f"""
## List of actions for which there are not exemplers of all keywords.  

There are {len(unexempled_keywords)} action pages which don't have examples that illustrate how to use all the keywords. On these documentation pages keywords that do not appear in at least one example are shown in red in the syntax table.
"""
          f.write( text )
          printDataTable( f, ["Action","Module"], unexempled_keywords )     

def printIndexFile(of,version) :
    content=f"""
PLUMED Version {version}
------------------------

PLUMED is a community-developed code that can be used to incorporate additional functionality into multiple molecular dynamics codes and for analysing
trajectories. PLUMED is currently interfaced with the list of codes described [here](mdcodes.md).

PLUMED is a composed of a [modules](modules.md) that contain a variety of different functionalities but that share a common basic syntax. You can find
a list of the modules that are available within PLUMED [here](modules.md) or you can see a graphical view of the modules and the dependencies between them [here](modulegraph.md).

Each module contains implementations of a number of [actions](actions.md), [shortcuts](shortcuts.md) and [command line tools](module_cltools.md). 
You can find a list of all the commands that you can use in a PLUMED input file [here](actionlist.md) and descriptions of some of the tools you can use these commands with [here](module_cltools.md).

Please also note that some developers prefer not to include their codes in PLUMED.  To use functionality that has been written by these developed you can use the [LOAD](LOAD.md) command.

You can find instructions for installing PLUMED [here](https://www.plumed-tutorials.org/lessons/20/001/data/NAVIGATION.html).

To run PLUMED you need to provide one input file.  If you are completely unfamiliar with PLUMED we would recommend that you start by working through 
[the following tutorial](https://www.plumed-tutorials.org/lessons/21/001/data/NAVIGATION.html) or the following [10-minute video](http://www.youtube.com/watch?v=PxJP16qNCYs).

You can find many other tutorials for PLUMED [here](https://www.plumed-tutorials.org) and you can find examples of how PLUMED has been used in many academic research articles [here](https://www.plumed-nest.org).

If you would like to add new functionality to PLUMED you can find developer documentation [here](../../developer-doc/html/index.html) and a change log [here](changelog.md).

The documentation in this manual was built on [{date.today().strftime('%B %d, %Y')}](summary.md).
    """
    of.write(content)

def printChangeLog(clf) :
    content=f"""
Change Log
----------

This page contains the history of changes across different PLUMED versions.
The future releases are expected to follow more or less the pace
of the old release. This means:
- Approximately one release per year, after summer, a new release (2.X). These releases
  typically group together all the features that were contributed during the
  year.
- Approximately every three months, we announce a patch (e.g. 2.2.X).
  These releases typically contains bug fixes, and could occasionally contain a new feature.

A few months before each new release we provide a beta release.
We typically maintain release branches until the fifth patch release (2.X.5),
which should come out approximately 15 month after the original release (2.X).
After that, branches are not supported anymore.

Notice that occasionally we publish patches on the mailing list.
These patches are always included in the following release, but we encourage
users that want to be up to date to follow the mailing list.

Below you can find change logs for all the published releases.
We mostly add new features without breaking existing ones.
However, some of the changes lead to incompatible behavior.
In the Change Log we try to give as much visibility as possible to these changes
to avoid surprises.

We also log changes that are relevant if you are developing the code. However, these
change lists are not complete, and if you want to put your hands in the code
and maintain your own collective variables we suggest you to follow the development
on [github](https://github.com/plumed/plumed2).

"""
    clf.write(content)
    for version in glob.glob("../CHANGES/*.md") :
        shutil.copy(version, "docs/" + version.split("/")[-1] )
        clf.write("- Changes for [Version " + version.split("/")[-1].replace("v","").replace(".md","") + "](" + version.split("/")[-1] + ")\n")       

def printActionListPage(af,version,tabledata,tagdictionary) :
    content=f"""
Actions implemented in PLUMED Version {version}
-----------------------------------------------

The [actions](actions.md) that can be used within a PLUMED input file are listed below.

"""
    af.write(content)
    printDataTable(af,["Name", "Description", "Tags"],
                   tabledata,
                   tagdictionary,
                   extraAttributes={"Name":{"class":f"actionHeader"}})

def printModuleListPage(mf,version,tabledata) :
    content=f"""
Modules that make up PLUMED Version {version}
---------------------------------------------
       
The functionality in PLUMED is split up between a collection of modules.  Each of these modules contains a collection of
[actions](actions.md), [shortcuts](shortcuts.md) and [command line tools](module_cltools.md).  Some of these modules must be present
every time PLUMED is compiled, others are not essential but are still compiled by default unless you explicitly tell PLUMED not to compile by using:
          
```bash
./configure --disable-module=module-name
```
   
The remainder of the modules are not compiled unless you explicitly request PLUMED to include that module in your compilation by using the command:
   
```bash
./configure --enable-module=module-name
```

To enable or disable multiple modules one should provide them as a : separated list. Notice that `+modulename` and `modulename` both activate the module, whereas
`-modulename` deactivates it. E.g.

```bash
./configure --enable-modules=+adjmat:-colvar
```

will disable the colvar module and enable the adjmat module.  The : can, in fact, be ommitted when you use + and -.  In other words, the following command can 
be used in place of the previous one:

```bash
./configure --enable-modules=+adjmat-colvar
```

If you repeat the `--enable-modules` keyword only the last instance will be used. Thus `./configure --enable-modules=adjmat --enable-modules=-colvar` will _not_ do what you expect!

!!! note "old implementation"

    Until PLUMED 2.2, it was also possible to switch on or off modules by adding files
    in the plumed2/src directory. Since PLUMED 2.3 this is discouraged, since any choice made
    in this manner will be overwritten next time `./configure` is used.

There are also some shortcuts available:

- `./configure --enable-modules=all` can be used to enable all optional modules. This includes the maximum number of features in PLUMED, including modules that might not be properly functional.
- `./configure --enable-modules=none` or `./configure --disable-modules` can be used to disable all optional modules. This produces a minimal PLUMED which can be used as a library but which has no command line tools and no collective variables or biasing methods.
- `./configure --enable-modules=reset` or `./configure --enable-modules` can be used to enable the default modules.

The two kinds of syntax can be combined and, for example, `./configure --enable-modules=none:colvar` will cause a version of PLUMED with all the modules disabled with the exception of the colvar module to be compiled.
   
The table below lists all the available modules and tells you whether they are always compiled, on by default or off by default.  An alternative, graphical
view of this information is available [here](modulegraph.md).

"""
    mf.write(content)
    printDataTable(mf,
                   ["Name","Description","Authors","Type"],
                   tabledata)

def addSpecialGroupsToPage( file, groups ) :
    with open(file,"r") as f : lines = f.readlines() 
    with open(file,"w+") as of : 
         for l in lines :
             if "@MOLINFOGROUPS@" in l : 
                printDataTable( of, ["Name","Description"], groups )
             else : of.write( l ) 

def getModuleType(module): 
    if not Path("../../src/" + module + "/module.type").is_file() : return "always"
    with open("../../src/" + module + "/module.type") as f : return f.read().strip()

def drawModuleNode( index, key, of ) :
    of.write(  str(index) + "(\"" + key + "\")\n")
    if getModuleType(key)=="always" : of.write("style " + str(index) + " fill:blue\n")
    elif getModuleType(key)=="default-on" : of.write("style " + str(index) + " fill:green\n")
    elif getModuleType(key)=="default-off" : of.write("style " + str(index) + " fill:red\n")    
    else : raise Exception("don't know how to draw node of type " + getModuleType(key) )

def getModuleRequirementsFromMakefile(thismodule) : 
   with open("../src/" + thismodule + "/Makefile") as f :
        for line in f.readlines() :
            if "USE=" in line : return line.replace("USE=","").split()
   return []

def createModuleGraph( version, plumed_syntax ) :
   # Get all the module dependencies from plumed_syntax
   requires = {}
   for key, value in plumed_syntax.items() :
       if "module" not in value :
          continue
       thismodule = value["module"]
       if thismodule not in requires.keys() :
          requires[thismodule] = set()
       if "needs" in value :
          for req in value["needs"] :
              if plumed_syntax[req]["module"]!=thismodule :
                 requires[thismodule].add( plumed_syntax[req]["module"] )
   # Get dependencies from Makefiles
   for thismodule in requires.keys() :
       for conn in getModuleRequirementsFromMakefile(thismodule) :
           if conn in requires.keys() :
              requires[thismodule].add( conn ) 

   with cd("docs") :
      with open( "modulegraph.md", "w") as of :
        ghead = f"""
PLUMED Version {version}
------------------------

PLUMED is a composed of a modules that contain a variety of different functionalities but that share a common basic syntax. You can find 
all the modules that are available within PLUMED in the following graph. The graph also shows the interdependencies between the various modules. 
If you click on the modules in the graph module-specific information will open.  The colors of the nodes in the graph below indicate whether the module
is always compiled (blue), on by default (green) or off by default (red).  If you need a feature from a module that is by default off you need to explicitly tell
PLUMED to include it during the configure stage by using:

```bash
./configure --enable-module=module-name
```

Each module contains implementations of a number of [actions](actions.md). You can find a list of all the actions implemented in in PLUMED [here](actionlist.md).

You can view the information about the modules in the graph above in a table by following [the following link](modules.md). 

<pre class=\"mermaid\">
"""
        of.write(ghead + "\n")
        of.write("%%{init: {\"flowchart\": {\"defaultRenderer\": \"elk\"}} }%%\n")
        of.write("flowchart TD\n")
   
        k, translate, backtranslate = 0, {}, []
        for key, data in requires.items() :
            translate[key] = k
            backtranslate.append(key) 
            k = k + 1
        
        # And create the graph
        G = nx.DiGraph()
        for key, data in requires.items() :
            for dd in data : G.add_edge( translate[dd], translate[key] )

        # Find any closed loops in the graph and remove them
        cycles = list( nx.simple_cycles(G) )
        for cyc in cycles :
           for i in range(len(cyc)-1) : G.remove_edge( cyc[i], cyc[(i+1)%len(cyc)] )   

        # And create the graph showing the modules
        pG = nx.transitive_reduction(G)

        # Create a matrix with the connections
        graphmat = np.zeros([k,k])
        for edge in pG.edges() : graphmat[edge[0],edge[1]] = 1
        for cyc in cycles : 
            for i in range(len(cyc)-1) : graphmat[cyc[i], cyc[(i+1)%len(cyc)]] = 1

        drawn = np.zeros(k)
        for i in range(k) : 
            if backtranslate[i]=="core" : continue
           
            group = set([i])
            for j in range(k) :
                if np.sum(graphmat[:,i])>0 and np.all(graphmat[:,j]==graphmat[:,i]) and drawn[j]==0 : group.add(j)

            # This code ensures that if there are more than 2 nodes that have identical dependencies we draw them in 
            # a subgraph.  The resulting flow chart is less clustered with arrows       
            if len(group)>2 : 
               of.write("subgraph g" + str(i) + " [ ]\n")
               ncols, lgroup, row, col = 3, [], 0, 0 
               for j in group :  
                   lgroup.append(j)
                   if drawn[j]==0 :
                      drawModuleNode( j, backtranslate[j], of )
                      if row>0 :
                         ind = lgroup[(row-1)*ncols + col]
                         of.write( str(ind) + "~~~" + str(j) + ";\n")
                      col = col + 1
                      if col%ncols==0 : col, row = 0, row + 1 
                      drawn[j]=1
               of.write("end\n")
               for l in range(k) :
                   if graphmat[l,j]>0 :
                      if drawn[l]==0 :
                         drawModuleNode( l, backtranslate[l], of )
                         drawn[l]=1
                      of.write( str(l) + "--> g" + str(i) + ";\n" )
               for j in group : graphmat[:,j] = 0

        for i in range(k) :
            if drawn[i]==0 : drawModuleNode( i,  backtranslate[i], of ) 
            for j in range(k) :
                if graphmat[i,j]>0 : of.write( str(i) + "-->" + str(j) + ";\n" )

        # And finally the click stuff
        k=0
        for key, data in requires.items() :
            mod_dict = getModuleDictionary( key )
            of.write("click " + str(k) + " \"../module_" + key + "\" \"" + mod_dict["description"] + "[Authors:" + mod_dict["authors"] + "]\"\n" )
            k = k + 1
        
        of.write("</pre>\n")

def getModuleDictionary( modname ) :
    if os.path.exists("../../src/" + modname + "/module.yml") :
       with open("../../src/" + module + "/module.yml") as f : moddict = yaml.load(f,Loader=yaml.BaseLoader)
       return moddict
    return {"name": modname, "authors": "authors", "description": "Information about the module", "dois": [] }

def createModulePage( version, modname, mod_dict, neggs, nlessons, plumed_syntax, plumedtags, tagdictionary, broken_inputs ) :
    with open( "module_" + modname + ".md", "w") as f :
         f.write("# [Module](modules.md): " + modname + "\n\n")
         f.write("| Description    | Usage |\n")
         f.write("|:--------|:--------:|\n") 
         f.write("| " + mod_dict["description"] + "\n __Authors:__ " + mod_dict["authors"] + " | ")
         if nlessons>0 :
            f.write("[![used in " + str(nlessons) + " tutorials](https://img.shields.io/badge/tutorials-" + str(nlessons) + "-green.svg)](https://www.plumed-tutorials.org/browse.html?module=" + modname + ")")
         else : 
            f.write("![used in " + str(nlessons) + " tutorials](https://img.shields.io/badge/tutorials-0-red.svg)")
         if neggs>0 : 
            f.write("[![used in " + str(neggs) + " eggs](https://img.shields.io/badge/nest-" + str(neggs) + "-green.svg)](https://www.plumed-nest.org/browse.html?module=" + modname + ")")
         else : 
            f.write("![used in " + str(neggs) + " eggs](https://img.shields.io/badge/nest-0-red.svg)")
         f.write("|\n\n")
         f.write("## Details \n") 
         if os.path.exists("../../src/" + modname + "/module.md") :
            with open("../../src/" + modname + "/module.md") as iff : docs = iff.read()
            actions = set()
            ninp, nf = processMarkdownString( docs, "module_" + modname + ".md", (PLUMED,), (version,), actions, f, ghmarkdown=False ) 
            if nf[0]>0 : broken_inputs.append( ["<a href=\"../module_" + modname + "\">" + modname + "</a>", str(nf[0])] )
         foundaction=False
         for key, value in plumed_syntax.items() :
             if key=="modules" or key=="vimlink" or key=="replicalink" or key=="groups" or key=="cltools" or key!=value["displayname"] or value["module"]!=modname : continue
             foundaction=True
             break 

         if modname=="cltools" : 
            if foundaction : raise Exception("Found action in cltools module. There should be no actions in this module")

            f.write("## Available tools \n\n")
            f.write("The following command line tools are available within PLUMED\n\n")
            f.write("| Name | Module | Description |\n")
            f.write("|:-----|:----|:------------|\n")
            for key, value in plumed_syntax["cltools"].items() :
                f.write("| [" + key + "](" + key + ".md) | [" + value["module"] + "](module_" + value["module"] + ".md) | " + value["description"] + "|\n")  
         else : 
            if foundaction : 
               f.write("## Actions \n\n")
               f.write("The following actions are part of this module\n\n")
               tabledata = []
               for key, value in plumed_syntax.items() :
                   if key=="modules" or key=="vimlink" or key=="replicalink" or key=="groups" or key=="cltools" or key!=value["displayname"] or value["module"]!=modname : continue
                   alink = "<a href=\"../" + key + "\">" + key + "</a>" 
                   tabledata.append( [alink, str(value["description"]), str(plumedtags[key])] ) 
               printDataTable(f, ["Name", "Description", "Tags"], tabledata, tagdictionary)

            foundcltool=False
            for key, value in plumed_syntax["cltools"].items() :
                if value["module"]!=modname : continue
                foundcltool=True
                break

            if foundcltool :
               f.write("\n## Command line tools \n\n")
               f.write("The following command line tools are part of this module\n\n")
               f.write("| Name | Description |\n")
               f.write("|:-----|:------------|\n")
               for key, value in plumed_syntax["cltools"].items() :
                   if value["module"]!=modname : continue
                   f.write("| [" + key + "](" + key + ".md) |" + value["description"] + "|\n")

         dois = mod_dict["dois"] 
         if len(dois)>0 : 
            f.write("\n## References \n")
            f.write("More information about this module is available in the following articles:\n\n")
            for doi in dois :
                ref, ref_url = get_reference(doi)
                f.write("- [" + ref + "](" + ref_url + ")\n")         

def getKeywordDescription( docs ) :
    desc = docs["description"] 
    if "actionlink" in docs.keys() and docs["actionlink"]!="none" :
       desc = desc + ". Options for this keyword are explained in the documentation for [" + docs["actionlink"] + "](" + docs["actionlink"] + ".md)."
    return desc

def createCLToolPage( version, tool, value, plumeddocs, broken_inputs, undocumented_keywords, noexamples, nodocs ) :
    with open( "docs/" + tool + ".md", "w") as f :
         f.write("# [Command line tool](module_cltools.md): " + tool + "\n\n")
         f.write("| [**Module**](modules.md) | [**" + value["module"] + "**](module_" + value["module"]  + ".md) |\n")
         f.write("|:--------|:--------:|\n")
         f.write("| **Description**    | **Input** |\n")
         f.write("| " + value["description"] + " | " + value["inputtype"] + "|\n\n" )
         f.write("## Details \n")
         if tool in plumeddocs.keys() :
            if os.path.isfile("../src/" + value["module"] + "/module.yml") :
                actions = set()
                ninp, nf = processMarkdownString( plumeddocs[tool], "docs/" + tool + ".md", (PLUMED,), (version,), actions, f, ghmarkdown=False )
                if nf[0]>0 : broken_inputs.append( ["<a href=\"../" + tool + "\">" + tool + "</a>", str(nf[0])] )
                if ninp==0 : noexamples.append( ["<a href=\"../" + tool + "\">" + tool + "</a>", value["module"]] )
            else :
                raise Exception("could not find yml file for module " + value["module"])
         else : nodocs.append(["<a href=\"../" + tool + "\">" + tool + "</a>", "cltool"] )
         f.write("\n## Syntax \n")
         if value["inputtype"]=="file" : 
             f.write("The following table describes the keywords that should be used in the input file for this command line tool\n\n")
         else :
             f.write("The following table describes the command line options that are available for this tool\n\n")
         f.write("| Keyword     | Description |\n")
         f.write("|:------------|:-----------|\n")
         undoc = 0
         for key, docs in value["syntax"].items() :
             if len(docs["description"])==0 :
                undoc = undoc + 1 
             f.write("| " + key + " | " + docs["description"] + " |\n")
         if undoc>0 :
            undocumented_keywords.append( ["<a href=\"../" + action + "\">" + action + "</a>", value["module"]] )
         if "dois" in value and len(value["dois"])>0 :
            f.write("\n## References \n")
            f.write("More information about how this cltool can be used is available in the following articles:\n")
            for doi in value["dois"] :
                ref, ref_url = get_reference(doi)
                f.write("- [" + ref + "](" + ref_url + ")\n") 

def createActionPage( version, action, value, plumeddocs, neggs, nlessons) :
    undocumented_keyword = None
    broken_input = None
    noexample= None
    nodoc = None
    unexempled_keywords = None
    with open( action + ".md", "w") as f : 
         hasatoms= False
         hasargs = False
         for key, docs in value["syntax"].items() :
             if key=="output" :
                continue
             if docs["type"]=="atoms" :
                hasatoms=True
             elif "argtype" in docs.keys() :
                hasargs=True

         if "IS_SHORTCUT" in value["syntax"].keys() : 
            f.write("# [Shortcut](shortcuts.md): " + action + "\n\n")
         else :
            f.write("# [Action](actions.md): " + action + "\n\n")

         f.write("| [**Module**](modules.md) | [**" + value["module"] + "**](module_" + value["module"]  + ".md) |\n")
         f.write("|:--------|:--------:|\n")
         f.write("| **Description**    | **Usage** |\n")
         f.write("| " + value["description"] + " | ")
         if nlessons>0 : 
            f.write("[![used in " + str(nlessons) + " tutorials](https://img.shields.io/badge/tutorials-" + str(nlessons) + "-green.svg)](https://www.plumed-tutorials.org/browse.html?action=" + action + ")")
         else : 
            f.write("![used in " + str(nlessons) + " tutorials](https://img.shields.io/badge/tutorials-0-red.svg)")
         if neggs>0 : 
            f.write("[![used in " + str(neggs) + " eggs](https://img.shields.io/badge/nest-" + str(neggs) + "-green.svg)](https://www.plumed-nest.org/browse.html?action=" + action + ")")
         else : 
            f.write("![used in " + str(neggs) + " eggs](https://img.shields.io/badge/nest-0-red.svg)") 
         if "output" in value["syntax"] and "value" in value["syntax"]["output"] : 
            f.write("|\n | **output value** | **type** |\n")
            f.write("| " + value["syntax"]["output"]["value"]["description"] + " | " + value["syntax"]["output"]["value"]["type"] + " |\n\n" )
         elif "output" not in value["syntax"] and (hasatoms or hasargs) :
            f.write("|\n **This action outputs data to a file**. You can read more about how PLUMED manages output files [here](Files.md) | |\n\n" ) 
         else : 
            f.write(" | \n\n")

         depracated = False
         if "replacement" in value :
            depracated = True 
            f.write("!!! warning \"Deprecated\"\n\n")
            f.write("    This action has been deprecated and is no longer supported. Use [" + value["replacement"] + "](" + value["replacement"] + ".md) instead.\n\n") 

         # Build a list of the keywords for this action that we want to see in the documentation
         example_keywords = set({})
         for key, docs in value["syntax"].items() :
             if key=="output" or docs["type"]=="hidden" : continue
             example_keywords.add(key)

         f.write("## Details and examples \n")
         if action in plumeddocs.keys() :
            if os.path.isfile("../../src/" + value["module"] + "/module.yml") :
                actions = set()
                print("BEFORE EXAMPLE", action, example_keywords )
                _, nf = processMarkdownString( plumeddocs[action], action + ".md", (PLUMED,), (version,), actions, f, ghmarkdown=False, checkaction=action, checkactionkeywords=example_keywords )
                print("AFTER EXAMPLE", action, example_keywords )
                if len(example_keywords)>0 :
                   unexempled_keywords = ["<a href=\"../" + action + "\">" + action + "</a>", value["module"]]  
                if nf[0]>0 :
                   broken_input = ["<a href=\"../" + action + "\">" + action + "</a>", str(nf[0])]
                if not depracated and action not in actions :
                   noexample = ["<a href=\"../" + action + "\">" + action + "</a>", value["module"]]
            else :
                raise Exception("could not find documentation for action " + action )
         else :
            nodoc = ["<a href=\"../" + action + "\">" + action + "</a>", "action"]

         if hasatoms or hasargs :
            f.write("## Input\n\n")
            if hasatoms and hasargs : f.write("The [arguments](specifying_arguments.md) and [atoms](specifying_atoms.md) that serve as the input for this action are specified using one or more of the keywords in the following table.\n\n")
            elif hasatoms : f.write("The [atoms](specifying_atoms.md) that serve as the input for this action are specified using one or more of the keywords in the following table.\n\n")
            elif hasargs : f.write("The [arguments](specifying_arguments.md) that serve as the input for this action are specified using one or more of the keywords in the following table.\n\n")

            f.write("| Keyword |  Type | Description |\n")
            f.write("|:--------|:------:|:-----------|\n")
            for key, docs in value["syntax"].items() :
                if key=="output" : continue
                if docs["type"]=="atoms" : f.write("| " + key + " | atoms | " + docs["description"] + " |\n")
                elif "argtype" in docs.keys() : f.write("| " + key + " | " + docs["argtype"] + " | " + docs["description"] + " |\n")
            f.write("\n\n")

         if "output" in value["syntax"] and len(value["syntax"]["output"].keys())>1 :
            f.write("## Output components\n\n")
            if "value" in value["syntax"]["output"] and len(value["syntax"]["output"])==1 :
               pass
            else :
               onlydefault = True
               for key, docs in value["syntax"]["output"].items() :
                   if docs["flag"]!="default" :
                     onlydefault = False
               if onlydefault :
                  f.write("This action calculates the [values](specifying_arguments.md) in the following table.  These [values](specifying_arguments.md) can be referenced elsewhere in the input by using this Action's label followed by a dot and the name of the [value](specifying_arguments.md) required from the list below.\n\n")
                  f.write("| Name | Type | Description |\n")
                  f.write("|:-------|:-----|:-------|\n")
                  for key, docs in value["syntax"]["output"].items() :
                      if key=="value" :
                         continue 
                      f.write("| " + key + " | " + docs["type"] + " | " + docs["description"] + " | \n") 
                  f.write("\n\n")
               else : 
                  f.write("This action can calculate the [values](specifying_arguments.md) in the following table when the associated keyword is included in the input for the action. These [values](specifying_arguments.md) can be referenced elsewhere in the input by using this Action's label followed by a dot and the name of the [value](specifying_arguments.md) required from the list below.\n\n")
                  f.write("| Name | Type | Keyword | Description |\n")
                  f.write("|:-------|:-----|:----:|:-------|\n")
                  for key, docs in value["syntax"]["output"].items() :
                      if key=="value" :
                         continue 
                      f.write("| " + key + " | " + docs["type"] + " | " + docs["flag"] + " | " + docs["description"] + " | \n")
                  f.write("\n\n")

         f.write("\n## Full list of keywords \n")
         f.write("The following table describes the [keywords and options](parsing.md) that can be used with this action \n\n")
         f.write("| Keyword | Type | Default | Description |\n")
         f.write("|:-------|:----:|:-------:|:-----------|\n")
         undoc = 0
         for key, docs in value["syntax"].items() :
             if key=="output" : continue 
             if "argtype" in docs.keys() and "default" in docs.keys() : 
                if len(docs["description"])==0 :
                  undoc = undoc + 1
                if key in example_keywords : f.write("| <span style=\"color:red\">" + key + "</span> | " + docs["type"] + " | " + docs["description"] + " | \n")
                else : f.write("| " + key + " | input | " + docs["default"] + " | " + getKeywordDescription( docs ) + " |\n")
             elif docs["type"]=="atoms" or "argtype" in docs.keys() :
                if len(docs["description"])==0 :
                  undoc = undoc + 1 
                if key in example_keywords : f.write("| <span style=\"color:red\">" + key + "</span> | " + docs["type"] + " | " + docs["description"] + " | \n")
                else : f.write("| " + key + " | input | none | " +  getKeywordDescription( docs ) + " |\n") 
         for key, docs in value["syntax"].items() : 
             if key=="output" or "argtype" in docs.keys()  : continue
             if docs["type"]=="compulsory" and "default" in docs.keys()  : 
                if len(docs["description"])==0 :
                  undoc = undoc + 1
                if key in example_keywords : f.write("| <span style=\"color:red\">" + key + "</span> | " + docs["type"] + " | " + docs["description"] + " | \n")
                else : f.write("| " + key + " | compulsory | "  + docs["default"] + " | " + getKeywordDescription( docs ) + " |\n") 
             elif docs["type"]=="compulsory" : 
                if len(docs["description"])==0 :
                  undoc = undoc + 1
                if key in example_keywords : f.write("| <span style=\"color:red\">" + key + "</span> | " + docs["type"] + " | " + docs["description"] + " | \n")
                else : f.write("| " + key + " | compulsory | none | " + getKeywordDescription( docs ) + " |\n")
         for key, docs in value["syntax"].items() :
             if key=="output" or "argtype" in docs.keys() : continue
             if docs["type"]=="flag" : 
                if len(docs["description"])==0 :
                  undoc = undoc + 1
                if key in example_keywords : f.write("| <span style=\"color:red\">" + key + "</span> | " + docs["type"] + " | " + docs["description"] + " | \n")
                else : f.write("| " + key + " | optional | false | " + getKeywordDescription( docs ) + " |\n")
             if docs["type"]=="optional" :
                if len(docs["description"])==0 :
                  undoc = undoc + 1 
                if key in example_keywords : f.write("| <span style=\"color:red\">" + key + "</span> | " + docs["type"] + " | " + docs["description"] + " | \n")
                else : f.write("| " + key + " | optional | not used | " + getKeywordDescription( docs ) + " |\n")
         if undoc>0 :
            undocumented_keyword = ["<a href=\"../" + action + "\">" + action + "</a>", value["module"]]

         if "dois" in value and len(value["dois"])>0 :
            f.write("\n## References \n")
            f.write("More information about how this action can be used is available in the following articles:\n\n")
            for doi in value["dois"] :
                ref, ref_url = get_reference(doi)
                f.write("- [" + ref + "](" + ref_url + ")\n")
    return broken_input, undocumented_keyword, noexample, nodoc, unexempled_keywords

def actionPage(key,plumed_syntax,nest_map,school_map,plumedtags,version,plumeddocs) :
    value = plumed_syntax[key]
    # Now create the page contents
    neggs =  nest_map[key] if key in nest_map.keys() else 0
    nlessons = school_map[key] if key in school_map.keys() else 0
    print("Building action page", key )
    with cd("docs") :
       broken_input, undocumented_keyword, noexample, nodoc, unexempled_keywords =createActionPage( version, key, value, plumeddocs, neggs, nlessons,
                         )
    alink = f'<a href="../{key}">{key}</a>'
    return {"broken_input" :broken_input,
             "undocumented_keyword" :undocumented_keyword,
             "noexample" :noexample,
             "nodoc" :nodoc,
             "unexempled_keywords" :unexempled_keywords,
             "element": [alink, str(value["description"]), str(plumedtags[key])]}

class InvalidJSONError(Exception):
      """Raised when the JSON file is not valid."""
      def __init__(self, message):
         super().__init__(f"Invalid JSON: {message}")

if __name__ == "__main__" :    # Get the version of plumed that we are building the manual for
   version = subprocess.check_output(f'../src/lib/plumed info --version', shell=True).decode('utf-8').strip()

   # Create dictionaries that hold how often each action has been used in nest/tutorials
   nest_map = create_map("https://www.plumed-nest.org/summary.html")
   school_map = create_map("https://www.plumed-tutorials.org/summary.html")

   # Get list of plumed actions from syntax file
   keyfile = "../json/syntax.json"
   with open(keyfile) as f :
       try:
          plumed_syntax = json.load(f)
       except ValueError as ve:
          raise InvalidJSONError(str(ve))

   # Create the index file
   with open("docs/index.md", "w+") as of : 
      printIndexFile(of,version)
   # Create the changelog file
   with open("docs/changelog.md", "w+") as clf :
      printChangeLog(clf)
   # Copy the extra files we need to process all the inputs 
   shutil.copytree("extras","docs/extras")
   # Copy the figures 
   shutil.copytree("figures","docs/figures") 
   # Create the assets
   shutil.copytree("assets","docs/assets")
   # Create the javascript
   with open("docs/assets/plumedtohtml.js", "w+") as jf :
      jf.write( get_javascript() )
   # Create the css
   with open("docs/assets/plumedtohtml.css", "w+") as cf :
      cf.write( get_css() )

   # Holder for table with pages with broken input
   broken_inputs = []

   # Create the general pages
   actions = set()
   for page in glob.glob("*.md") :
       shutil.copy( page, "docs/" + page ) 
       with cd("docs") : 
          ninp, nf = processMarkdown( page, (PLUMED,), (version,), actions, ghmarkdown=False )
       if nf[0]>0 :
          broken_inputs.append( ["<a href=\"../" + page.replace(".md","") + "\">" + page + "</a>", str(nf[0])] )
       if os.path.exists("docs/colvar") :
          os.remove("docs/colvar")   # Do this with Plumed2HTML maybe

   # Create a page for each action
   plumeddocs, plumedtags = getActionDocumentationFromPlumed(plumed_syntax)

   actionKeys=[]
   notActionList = ["modules", "vimlink", "replicalink", "groups", "cltools"]
   
   biggestkey=[0,""]
   for key, value in plumed_syntax.items() :
      if key in notActionList or key!=value["displayname"] :
        continue
      actionKeys.append(key)
      if len(key)>biggestkey[0] :
         biggestkey = [len(key), key]

   # here for debugging purposes:
   with open("debugging.txt", "w") as bf :
      bf.write(biggestkey[1] + " " + str(biggestkey[0]) + "\n")
      bf.write(f"nest_map={nest_map}\n")
      bf.write(f"school_map={school_map}\n")
      print("The biggest action is", biggestkey[1], "with length", biggestkey[0])
   #create a function that will fix all the inputs but key for action page
   run_action_page = functools.partial(actionPage, plumed_syntax=plumed_syntax, 
                                            nest_map=nest_map, school_map=school_map,
                                            plumedtags=plumedtags, version=version,
                                            plumeddocs=plumeddocs)
   if sys_platform == "wasi" or sys_platform == "ios":
      #multiprocessing is not avaiable on WASI or iOS
      #https://docs.python.org/3/library/multiprocessing.html
      actiontable = [run_action_page(key) for key in sorted(actionKeys)]
   else :
      with Pool(cpu_count()) as pool:  
        actiontable = pool.map(run_action_page, sorted(actionKeys))
   tabledata = [x["element"] for x in actiontable]
   #append
   broken_inputs+=[x["broken_input"] for x in actiontable if x["broken_input"] is not None]
   undocumented_keywords = [x["undocumented_keyword"] for x in actiontable if x["undocumented_keyword"] is not None]
   noexamples = [x["noexample"] for x in actiontable if x["noexample"] is not None]
   nodocs = [x["nodoc"] for x in actiontable if x["nodoc"] is not None]
   unexempled_keywords = [x["unexempled_keywords"] for x in actiontable if x["unexempled_keywords"] is not None]
   
   # Read in all the module pages
   tagdictionary = {}
   for file in glob.glob("../src/*/module.yml") :
       with open(file) as f :
          moddict = yaml.load(f,Loader=yaml.BaseLoader)
       if "tags" in moddict : 
          for key, value in moddict["tags"].items() :
              if key in tagdictionary.keys() :
                 raise Exception("found duplicate definitions for tag " + key) 
              tagdictionary[key] = value

   # Find the list of tags in the tag list
   unfound_tags = set()
   for key, value in plumedtags.items() : 
       for tag in value.split() :
           if tag not in tagdictionary.keys() :
              tagdictionary[tag] = "place holder text for tag"
              unfound_tags.add(tag)
   print( "THESE ARE THE TAGS WE DIDN'T FIND", unfound_tags )

   # Create the page with the list of actions
   with open("docs/actionlist.md","w+") as actdb :
      printActionListPage( actdb, version, tabledata, tagdictionary )

   # Create a page for each cltool
   for key, value in plumed_syntax["cltools"].items() :
       print("Building CLtool page", key )
       createCLToolPage( version, key, value, plumeddocs, broken_inputs, undocumented_keywords, noexamples, nodocs )

   # Create a list of modules
   modules = {}
   for key in actionKeys :  
     value = plumed_syntax[key]
     if key in notActionList or key!=value["displayname"] :
       continue
     neggs = nest_map[key] if key in nest_map.keys() else 0
     nlessons = school_map[key] if key in school_map.keys() else 0
     if value["module"] not in modules.keys() :
        modules[value["module"]] = { "neggs": neggs, "nlessons": nlessons }
     else :
        modules[value["module"]]["neggs"] += neggs
        modules[value["module"]]["nlessons"] += nlessons

   # Take into account command line tools when building list of modules
   for key, value in plumed_syntax["cltools"].items() :
     neggs =  nest_map[key] if key in nest_map.keys() else 0
     nlessons = school_map[key] if key in school_map.keys() else 0
     if value["module"] not in modules.keys() :
        modules[value["module"]] = { "neggs": neggs, "nlessons": nlessons }
     else :
        modules[value["module"]]["neggs"] += neggs
        modules[value["module"]]["nlessons"] += nlessons  

   # And create each module page
   moduletabledata = []
   for module, value in modules.items() :
       mlink = "<a href=\"../module_" + module + "\">" + module + "</a>"
       print("Building module page", module )
       with cd("docs") :
          mod_dict = getModuleDictionary( module )
          moduletabledata.append( [mlink, mod_dict["description"], mod_dict["authors"], getModuleType(module) ] )
          createModulePage( version, module, mod_dict, value["neggs"], value["nlessons"], plumed_syntax, plumedtags, tagdictionary, broken_inputs )
       if not os.path.exists("../src/" + module + "/module.md") or not os.path.exists("../src/" + module + "/module.yml") : 
          nodocs.append(["<a href=\"../module_" + module + "\">" + module+ "</a>", "module"])
   # And the page with the list of modules
   with open("docs/modules.md","w+") as module_file :
      printModuleListPage( module_file, version, moduletabledata )
   # Create the graph that shows all the modules
   createModuleGraph( version, plumed_syntax ) 
   
   # Create a list containing all the special groups
   special_groups = []
   for key, value in plumed_syntax["groups"].items() :
      special_groups.append([ key, str(value["description"]) ])
   # Add tables with special groups to pages that need them
   addSpecialGroupsToPage( "docs/specifying_atoms.md", special_groups )
   addSpecialGroupsToPage( "docs/MOLINFO.md", special_groups ) 

   # And output the summary page
   createSummaryPage( broken_inputs, nodocs, undocumented_keywords, noexamples, unexempled_keywords )
