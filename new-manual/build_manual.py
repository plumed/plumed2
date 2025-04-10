import os
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

PLUMED="plumed"

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def getActionDocumentationFromPlumed() :
    docs = {}
    for file in glob.glob("../src/*/*.cpp") : 
        with open(file,"r") as f : content = f.read()
        if "//+PLUMEDOC" not in content : continue
        actioname, founddocs,  indocs = "", "", False
        for line in content.splitlines() :
            if "//+ENDPLUMEDOC" in line :
               if not indocs : raise Exception("Found ENDPLUMEDDOC before PLUMEDOC in " + file)
               print("Found documentation for ", actioname)
               docs[actioname] = founddocs
               actionname, founddocs, indocs = "", "", False
            elif "//+PLUMEDOC" in line :
               if indocs : raise Exception("Found PLUMEDDOC before ENDPLUMEDOC in " + file) 
               actioname, indocs = line.split()[2], True 
            elif indocs and not "/*" in line and not "*/" in line : founddocs += line + "\n"
    return docs

def get_reference(doi):
    # initialize strings
    ref=""; ref_url=""
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
    script_elements = soup.find_all("script")
    xdata, ydata = {}, {} 
    for script in script_elements : 
        lines = script.text.splitlines()
        for line in lines :
            if "var" not in line or "for" in line : continue
            if "xValues" in line and "=" in line :
               xdata = json.loads( line.split("=")[1].replace(";","") )
            if "yValues" in line and "=" in line :
               ydata = json.loads( line.split("=")[1].replace(";","") )
    
    return dict(map(lambda i,j : (i,j) , xdata,ydata))

def printDataTable( f, titles, tabledata ) :
    f.write("<table id=\"browse-table\" class=\"display\">\n")
    f.write("<thead><tr>\n")
    for t in titles : f.write("<th style=\"text-align: left\">" + t + "</th>\n")
    f.write("</tr></thead><tbody>\n")
    for r in tabledata : 
        if len(r)!=len(titles) : raise Exception("mismatch between number of columns in tabledata and number of titles")
        f.write("<tr>\n")
        for e in r : f.write("<td style=\"text-align: left\">" + e + "</td>\n")
        f.write("</tr>\n")
    f.write("</tbody></table>\n")

def createSummaryPage( broken, nodocs, undocumented, noexamples ) :
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

def printIndexFile(of,version) :
    content=f"""
PLUMED Version {version}
------------------------

PLUMED is a community-developed code that can be used to incorporate additional functionality into multiple molecular dynamics codes and for analysing
trajectories. PLUMED is a composed of a [modules](modules.md) that contain a variety of different functionalities but that share a common basic syntax. You can find
a list of the modules that are available within PLUMED [here](modules.md) or you can see a graphical view of the modules and the dependencies between them [here](modulegraph.md).

Each module contains implementations of a number of [actions](actions.md), [shortcuts](shortcuts.md) and [command line tools](module_cltools.md). 
You can find a list of all the commands that you can use in a PLUMED input file [here](actionlist.md) and descriptions of some of the tools you can use these commands with [here](module_cltools.md).

Please also note that some developers prefer not to include their codes in PLUMED.  To use functionality that has been written by these developed you can use the [LOAD](LOAD.md) command.

If you are completely unfamiliar with PLUMED we would recommend that you start by working through [the following tutorial](https://www.plumed-tutorials.org/lessons/21/001/data/NAVIGATION.html).

You can find many other tutorials for PLUMED [here](https://www.plumed-tutorials.org) and you can find examples of how PLUMED has been used in many academic research articles [here](https://www.plumed-nest.org).

If you would like to add new functionality to PLUMED you can find developer documentation [here](../../developer-doc/html/index.html).

The documentation in this manual was built on [{date.today().strftime('%B %d, %Y')}](summary.md).
    """
    of.write(content)

def printActionListPage(af,version,tabledata) :
    content=f"""
Actions implemented in PLUMED Version {version}
-----------------------------------------------

The [actions](actions.md) that can be used within a PLUMED input file are listed below.

"""
    af.write(content)
    printDataTable(af,["Name", "Description"], tabledata)

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
   
The table below lists all the available modules and tells you whether they are always compiled, on by default or off by default.  An alternative, graphical
view of this information is available [here](modulegraph.md).

"""
    mf.write(content)
    printDataTable(mf,["Name","Description","Authors","Type"],tabledata)

def addSpecialGroupsToPage( file, groups ) :
    with open(file,"r") as f : lines = f.readlines() 
    with open(file,"w+") as of : 
         for l in lines :
             if "@MOLINFOGROUPS@" in l : 
                printDataTable( of, ["Name","Description"], groups )
             else : of.write( l ) 

def getModuleType(module): 
    if not Path("../src/" + module + "/module.type").is_file() : return "always"
    with open("../src/" + module + "/module.type") as f : return f.read().strip()

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
       if "module" not in value : continue
       thismodule = value["module"]
       if thismodule not in requires.keys() : requires[thismodule] = set()
       if "needs" in value :
          for req in value["needs"] :
              if plumed_syntax[req]["module"]!=thismodule : requires[thismodule].add( plumed_syntax[req]["module"] )
   # Get dependencies from Makefiles
   for thismodule in requires.keys() :
       for conn in getModuleRequirementsFromMakefile(thismodule) :
           if conn in requires.keys() : requires[thismodule].add( conn ) 

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
                      drawn[j]==1
               of.write("end\n")
               for l in range(k) :
                   if graphmat[l,j]>0 :
                      if drawn[l]==0 :
                         drawModuleNode( l, backtranslate[l], of )
                         drawn[l]==1
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
    if os.path.exists("../../src/" + module + "/module.yml") :
       with open("../../src/" + module + "/module.yml") as f : moddict = yaml.load(f,Loader=yaml.BaseLoader)
       return moddict
    return {"name": modname, "authors": "authors", "description": "Information about the module", "dois": [] }

def createModulePage( version, modname, neggs, nlessons, plumed_syntax, broken_inputs ) :
    with open( "module_" + modname + ".md", "w") as f :
         f.write("# [Module](modules.md): " + modname + "\n\n")
         f.write("| Description    | Usage |\n")
         f.write("|:--------|:--------:|\n") 
         f.write("| " + getModuleDictionary(modname)["description"] + "\n __Authors:__ " + mod_dict["authors"] + " | ")
         if nlessons>0 :
            f.write("[![used in " + str(nlessons) + " tutorials](https://img.shields.io/badge/tutorials-" + str(nlessons) + "-green.svg)](https://www.plumed-tutorials.org/browse.html?search=" + modname + ")")
         else : 
            f.write("![used in " + str(nlessons) + " tutorials](https://img.shields.io/badge/tutorials-0-red.svg)")
         if neggs>0 : 
            f.write("[![used in " + str(neggs) + " eggs](https://img.shields.io/badge/nest-" + str(neggs) + "-green.svg)](https://www.plumed-nest.org/browse.html?search=" + modname + ")")
         else : 
            f.write("![used in " + str(neggs) + " eggs](https://img.shields.io/badge/nest-0-red.svg)")
         f.write("|\n\n")
         f.write("## Details \n") 
         if os.path.exists("../../src/" + modname + "/module.md") :
            with open("../../src/" + modname + "/module.md") as iff : docs = iff.read()
            actions = set()
            ninp, nf = processMarkdownString( docs, "module_" + modname + ".md", (PLUMED,), (version,), actions, f, ghmarkdown=False ) 
            if nf[0]>0 : broken_inputs.append( ["<a href=\"../module_" + modname + "\">" + modname + "</a>", str(nf[0])] )
         dois = getModuleDictionary(modname)["dois"] 
         if len(dois)>0 : 
            f.write("## References \n")
            f.write("More information about this module is available in the following articles:\n\n") 
            for doi in dois :
                ref, ref_url = get_reference(doi)
                f.write("- [" + ref + "](" + ref_url + ")\n")
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
               f.write("| Name | Description |\n")
               f.write("|:-----|:------------|\n")
               for key, value in plumed_syntax.items() :
                   if key=="modules" or key=="vimlink" or key=="replicalink" or key=="groups" or key=="cltools" or key!=value["displayname"] or value["module"]!=modname : continue
                   f.write("| [" + key + "](" + key + ".md) |" + value["description"] + "|\n")

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
            else :
                f.write("Text from manual goes here \n")
                noexamples.append( ["<a href=\"../" + tool + "\">" + tool + "</a>", value["module"]] )
         else : nodocs.append(["<a href=\"../" + tool + "\">" + tool + "</a>", "cltool"] )
         if "dois" in value and len(value["dois"])>0 :
            f.write("## References \n")
            f.write("More information about how this cltool can be used is available in the following articles:\n")
            for doi in value["dois"] :
                ref, ref_url = get_reference(doi)
                f.write("- [" + ref + "](" + ref_url + ")\n")
         f.write("\n## Syntax \n")
         if value["inputtype"]=="file" : 
             f.write("The following table describes the keywords that should be used in the input file for this command line tool\n\n")
         else :
             f.write("The following table describes the command line options that are available for this tool\n\n")
         f.write("| Keyword     | Description |\n")
         f.write("|:------------|:-----------|\n")
         undoc = 0
         for key, docs in value["syntax"].items() :
             if len(docs["description"])==0 : undoc = undoc + 1 
             f.write("| " + key + " | " + docs["description"] + " |\n")
         if undoc>0 : undocumented_keywords.append( ["<a href=\"../" + action + "\">" + action + "</a>", value["module"]] )

def createActionPage( version, action, value, plumeddocs, neggs, nlessons, broken_inputs, undocumented_keywords, noexamples, nodocs ) :
    with open( action + ".md", "w") as f : 
         hasatoms, hasargs = False, False
         for key, docs in value["syntax"].items() :
             if key=="output" : continue
             if docs["type"]=="atoms" : hasatoms=True
             elif "argtype" in docs.keys() : hasargs=True

         if "IS_SHORTCUT" in value["syntax"].keys() : f.write("# [Shortcut](shortcuts.md): " + action + "\n\n")
         else : f.write("# [Action](actions.md): " + action + "\n\n")

         f.write("| [**Module**](modules.md) | [**" + value["module"] + "**](module_" + value["module"]  + ".md) |\n")
         f.write("|:--------|:--------:|\n")
         f.write("| **Description**    | **Usage** |\n")
         f.write("| " + value["description"] + " | ")
         if nlessons>0 : 
            f.write("[![used in " + str(nlessons) + " tutorials](https://img.shields.io/badge/tutorials-" + str(nlessons) + "-green.svg)](https://www.plumed-tutorials.org/browse.html?search=" + action + ")")
         else : 
            f.write("![used in " + str(nlessons) + " tutorials](https://img.shields.io/badge/tutorials-0-red.svg)")
         if neggs>0 : 
            f.write("[![used in " + str(neggs) + " eggs](https://img.shields.io/badge/nest-" + str(neggs) + "-green.svg)](https://www.plumed-nest.org/browse.html?search=" + action + ")")
         else : 
            f.write("![used in " + str(neggs) + " eggs](https://img.shields.io/badge/nest-0-red.svg)") 
         if "output" in value["syntax"] and "value" in value["syntax"]["output"] : 
            f.write("|\n | **output value** | **type** |\n")
            f.write("| " + value["syntax"]["output"]["value"]["description"] + " | " + value["syntax"]["output"]["value"]["type"] + " |\n\n" )
         elif "output" not in value["syntax"] and (hasatoms or hasargs) :
            f.write("|\n **This action outputs data to a file**. You can read more about how PLUMED manages output files [here](Files.md) | |\n\n" ) 
         else : 
            f.write(" | \n\n")

         if "output" in value["syntax"] and len(value["syntax"]["output"].keys())>1 :
            f.write("## Output components\n\n")
            if "value" in value["syntax"]["output"] and len(value["syntax"]["output"])==1 :
               pass
            else :
               onlydefault = True
               for key, docs in value["syntax"]["output"].items() :
                   if docs["flag"]!="default" : onlydefault = False
               if onlydefault :
                  f.write("This action calculates the [values](specifying_arguments.md) in the following table.  These [values](specifying_arguments.md) can be referenced elsewhere in the input by using this Action's label followed by a dot and the name of the [value](specifying_arguments.md) required from the list below.\n\n")
                  f.write("| Name | Type | Description |\n")
                  f.write("|:-------|:-----|:-------|\n")
                  for key, docs in value["syntax"]["output"].items() :
                      if key=="value" : continue 
                      f.write("| " + key + " | " + docs["type"] + " | " + docs["description"] + " | \n") 
                  f.write("\n\n")
               else : 
                  f.write("This action can calculate the [values](specifying_arguments.md) in the following table when the associated keyword is included in the input for the action. These [values](specifying_arguments.md) can be referenced elsewhere in the input by using this Action's label followed by a dot and the name of the [value](specifying_arguments.md) required from the list below.\n\n")
                  f.write("| Name | Type | Keyword | Description |\n")
                  f.write("|:-------|:-----|:----:|:-------|\n")
                  for key, docs in value["syntax"]["output"].items() :
                      if key=="value" : continue 
                      f.write("| " + key + " | " + docs["type"] + " | " + docs["flag"] + " | " + docs["description"] + " | \n")
                  f.write("\n\n")
         
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

         ninp, nfail = 0, 0
         f.write("## Further details and examples \n")
         if action in plumeddocs.keys() :
            if os.path.isfile("../../src/" + value["module"] + "/module.yml") : 
                actions = set()
                ninp, nf = processMarkdownString( plumeddocs[action], action + ".md", (PLUMED,), (version,), actions, f, ghmarkdown=False )
                if nf[0]>0 : broken_inputs.append( ["<a href=\"../" + action + "\">" + action + "</a>", str(nf[0])] )
            else : 
                f.write("Text from manual goes here \n")
                noexamples.append( ["<a href=\"../" + action + "\">" + action + "</a>", value["module"]] )
         else : nodocs.append(["<a href=\"../" + action + "\">" + action + "</a>", "action"] )
         if "dois" in value and len(value["dois"])>0 : 
            f.write("## References \n")
            f.write("More information about how this action can be used is available in the following articles:\n")
            for doi in value["dois"] :
                ref, ref_url = get_reference(doi)
                f.write("- [" + ref + "](" + ref_url + ")\n")
         f.write("\n## Syntax \n")
         f.write("The following table describes the [keywords and options](parsing.md) that can be used with this action \n\n")
         f.write("| Keyword | Type | Default | Description |\n")
         f.write("|:-------|:----:|:-------:|:-----------|\n")
         undoc = 0
         for key, docs in value["syntax"].items() :
             if key=="output" : continue 
             if "argtype" in docs.keys() and "default" in docs.keys() : 
                if len(docs["description"])==0 : undoc = undoc + 1
                f.write("| " + key + " | input | " + docs["default"] + " | " + getKeywordDescription( docs ) + " |\n")
             elif docs["type"]=="atoms" or "argtype" in docs.keys() :
                if len(docs["description"])==0 : undoc = undoc + 1 
                f.write("| " + key + " | input | none | " +  getKeywordDescription( docs ) + " |\n") 
         for key, docs in value["syntax"].items() : 
             if key=="output" or "argtype" in docs.keys()  : continue
             if docs["type"]=="compulsory" and "default" in docs.keys()  : 
                if len(docs["description"])==0 : undoc = undoc + 1
                f.write("| " + key + " | compulsory | "  + docs["default"] + " | " + getKeywordDescription( docs ) + " |\n") 
             elif docs["type"]=="compulsory" : 
                if len(docs["description"])==0 : undoc = undoc + 1
                f.write("| " + key + " | compulsory | none | " + getKeywordDescription( docs ) + " |\n")
         for key, docs in value["syntax"].items() :
             if key=="output" or "argtype" in docs.keys() : continue
             if docs["type"]=="flag" : 
                if len(docs["description"])==0 : undoc = undoc + 1
                f.write("| " + key + " | optional | false | " + getKeywordDescription( docs ) + " |\n")
             if docs["type"]=="optional" :
                if len(docs["description"])==0 : undoc = undoc + 1 
                f.write("| " + key + " | optional | not used | " + getKeywordDescription( docs ) + " |\n")
         if undoc>0 : undocumented_keywords.append( ["<a href=\"../" + action + "\">" + action + "</a>", value["module"]] )

if __name__ == "__main__" : 
   # Get the version of plumed that we are building the manual for
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
          raise InvalidJSONError(ve)

   # Create the index file
   with open("docs/index.md", "w+") as of : printIndexFile(of,version)
   # Copy the extra files we need to process all the inputs 
   shutil.copytree("extras","docs/extras")
   # Create the assets
   shutil.copytree("assets","docs/assets")
   # Create the javascript
   with open("docs/assets/plumedtohtml.js", "w+") as jf : jf.write( get_javascript() )
   # Create the css
   with open("docs/assets/plumedtohtml.css", "w+") as cf : cf.write( get_css() )

   # Holder for table with pages with broken inputs
   broken_inputs = []

   # Create the general pages
   actions = set()
   for page in glob.glob("*.md") :
       shutil.copy( page, "docs/" + page ) 
       with cd("docs") : 
          ninp, nf = processMarkdown( page, (PLUMED,), (version,), actions, ghmarkdown=False )
       if nf[0]>0 : broken_inputs.append( ["<a href=\"../" + page.replace(".md","") + "\">" + page + "</a>", str(nf[0])] )
       if os.path.exists("docs/colvar") : os.remove("docs/colvar")   # Do this with Plumed2HTML maybe

   # Create a page for each action
   plumeddocs, tabledata, undocumented_keywords, noexamples, nodocs = getActionDocumentationFromPlumed(), [], [], [], []
   for key, value in plumed_syntax.items() :
      if key=="modules" or key=="vimlink" or key=="replicalink" or key=="groups" or key=="cltools" or key!=value["displayname"] : continue
      # Now create the page contents
      neggs, nlessons = 0, 0
      if key in nest_map.keys() : neggs = nest_map[key]
      if key in school_map.keys() : nlessons = school_map[key] 
      print("Building action page", key )
      with cd("docs") : 
         createActionPage( version, key, value, plumeddocs, neggs, nlessons, broken_inputs, undocumented_keywords, noexamples, nodocs )
      alink = "<a href=\"../" + key + "\">" + key + "</a>"
      tabledata.append( [alink, str(value["description"])] ) 
   # Create the page with the list of actions
   with open("docs/actionlist.md","w+") as actdb : printActionListPage( actdb, version, tabledata )

   # Create a page for each cltool
   for key, value in plumed_syntax["cltools"].items() :
       print("Building CLtool page", key )
       createCLToolPage( version, key, value, plumeddocs, broken_inputs, undocumented_keywords, noexamples, nodocs )

   # Create a list of modules
   modules = {}
   for key, value in plumed_syntax.items() :  
     if key=="modules" or key=="vimlink" or key=="replicalink" or key=="groups" or key=="cltools" or key!=value["displayname"] : continue
     nlessons, neggs = 0, 0
     if key in school_map.keys() : nlessons = school_map[key]
     if key in nest_map.keys() : neggs = nest_map[key]
     if value["module"] not in modules.keys() :
        modules[value["module"]] = { "neggs": neggs, "nlessons": nlessons }
     else : modules[value["module"]]["neggs"], modules[value["module"]]["nlessons"] = modules[value["module"]]["neggs"] + neggs, modules[value["module"]]["nlessons"] + nlessons

   # Take into account command line tools when building list of modules
   for key, value in plumed_syntax["cltools"].items() :
     nlessons, neggs = 0, 0
     if key in school_map.keys() : nlessons = school_map[key] 
     if key in nest_map.keys() : neggs = nest_map[key]
     if value["module"] not in modules.keys() :
        modules[value["module"]] = { "neggs": neggs, "nlessons": nlessons }
     else : modules[value["module"]]["neggs"], modules[value["module"]]["nlessons"] = modules[value["module"]]["neggs"] + neggs, modules[value["module"]]["nlessons"] + nlessons  

   # And create each module page
   moduletabledata = []
   for module, value in modules.items() :
       mlink = "<a href=\"../module_" + module + "\">" + module + "</a>"
       print("Building module page", module )
       with cd("docs") :
          mod_dict = getModuleDictionary( module )
          moduletabledata.append( [mlink, mod_dict["description"], mod_dict["authors"], getModuleType(module) ] )
          createModulePage( version, module, value["neggs"], value["nlessons"], plumed_syntax, broken_inputs )
       if not os.path.exists("../src/" + module + "/module.md") or not os.path.exists("../src/" + module + "/module.yml") : 
          nodocs.append(["<a href=\"../module_" + module + "\">" + module+ "</a>", "module"])
   # And the page with the list of modules
   with open("docs/modules.md","w+") as module_file : printModuleListPage( module_file, version, moduletabledata )
   # Create the graph that shows all the modules
   createModuleGraph( version, plumed_syntax ) 
   
   # Create a list containing all the special groups
   special_groups = []
   for key, value in plumed_syntax["groups"].items() : special_groups.append([ key, str(value["description"]) ])
   # Add tables with special groups to pages that need them
   addSpecialGroupsToPage( "docs/specifying_atoms.md", special_groups )
   addSpecialGroupsToPage( "docs/MOLINFO.md", special_groups ) 

   # And output the summary page
   createSummaryPage( broken_inputs, nodocs, undocumented_keywords, noexamples )
