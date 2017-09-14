{
  if(match($0,"BEGIN_PLUMED_FILE")){
    inside=1;
    endplumed=0;
    sub("BEGIN_PLUMED_FILE","");
    print;
    next;
  }
  if(inside && match($0,"</pre>")){
    inside=0;
    print;
    next;
  }
  if(!inside){
    print;
    next;
  }

# DRAFT LINK TO DOC:
 copy=$0
 sub("#.*","",copy);
 if(endplumed) copy="";
 nw=split(copy,words);
 if(match(words[1],".*:$")){
   action=words[2];
 } else {
   action=words[1];
 }
 if(action=="__FILL__") action=""
 if(action=="ENDPLUMED"){
   endplumed=1;
 }
 actionx="";
 for(i=1;i<=length(action);i++){
   letter=substr(action,i,1);
   if(match(letter,"[A-Z]")) letter = "_" tolower(letter);
   actionx=actionx letter;
 }
 if(incontinuation) action="";
 if(incontinuation && words[1]=="...") incontinuation=0;
 else if(!incontinuation && words[nw]=="...") incontinuation=1;

 if(length(action)>0){
   actionfile="html/" actionx ".html";
   if(getline tmp < actionfile < 0) print "WARNING: file " actionfile " does not exist" > "/dev/stderr";
   else {
     p=match($0,action);
     $0=substr($0,1,p-1) "<a href=\"./" actionx ".html\" style=\"color:green\">" action "</a>" substr($0,p+length(action));
   }
 }

 gsub("__FILL__","<span style=\"background-color:yellow\">__FILL__</span>");

# comments:
#  sub("#","<span style=\"color:blue\">#");
#  if(match($0,"span style=")) $0=$0 "</span>";
  sub("#.*$","<span style=\"color:blue\">&</span>");
  if(endplumed) sub("^.*$","<span style=\"color:blue\">&</span>");
  print
}
