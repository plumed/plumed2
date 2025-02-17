import subprocess
import json

def get_dois_from_syntax() :
  with open("syntax.json") as f :
    try:
       plumed_syntax = json.load(f)
    except ValueError as ve:
       raise InvalidJSONError(ve)

  doilist = []
  for key, value in plumed_syntax.items() :
      if not isinstance(value,dict) or "dois" not in value.keys() : continue
      for doi in value["dois"] : 
          if doi not in doilist : doilist.append(doi)
  return doilist

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
        try:
          # get citation
          cit = subprocess.check_output('curl -LH "Accept: text/bibliography; style=science" \'http://dx.doi.org/'+doi+'\'', shell=True).decode('utf-8').strip()
          if("DOI Not Found" in cit):
           ref="DOI not found"
          else:
           # get citation
           ref=cit[3:cit.find(", doi")]
           # and url
           ref_url="http://dx.doi.org/"+doi 
        except:
          ref="DOI not found"
    return ref

if __name__ == "__main__" :
   # Get all the dois that are listed in the syntax file
   doilist = get_dois_from_syntax()
   # Now get all the references using Max's magic script and create the map
   of = open("../src/tools/CitationMap.inc","w+")
   for doi in doilist : 
       reference = get_reference(doi)
       if doi==doilist[-1] : of.write("{\"" + doi + "\",\"" + get_reference(doi) + "\"}") 
       else : of.write("{\"" + doi + "\",\"" + get_reference(doi) + "\"}")
   of.close()
