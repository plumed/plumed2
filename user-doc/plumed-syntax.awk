BEGIN{inside=0}{
  if(match($0,"BEGIN_PLUMED_FILE")){
    if(inside) {
       print "Found BEGIN_PLUMED_FILE before exit of previous" > "/dev/stderr"
       exit 1; 
    }
    inside=1;
    nreplicas=1;
    number=number+1;
    endplumed=0;
    status="broken";
    if(match($0,"working")){
       status="working";
    } else if(match($0,"loads")){
       status="loads";
    } else if(match($0,"incomplete")){
       status="incomplete";
    }
    datadir="";
    for(i=1;i<=NF;++i) {
       if( match($i,"DATADIR=") ) { datadir=$i; sub("DATADIR=","",datadir); }
    }
    if( datadir!="" ) { system("cp ../../"datadir"/* ."); }
    next;
  }
  if(inside && match($0,"</pre>")){
    inside=0; close("example"number".dat"); 
    if( nreplicas==1 ) { system("plumed gen_example --plumed example"number".dat --out example"number".html --name eg"number" --status " status "> /dev/null"); } 
    else { system("mpirun -np " nreplicas " plumed gen_example --plumed example"number".dat --out example"number".html --name eg"number" --status " status " --multi " nreplicas " > /dev/null");}
    system("cat example"number".html"); 
    system("rm example"number".dat example"number".html");
    sub("</pre>","");
    print;
    next;
  }
  if(!inside){
    print;
    next;
  }

# DRAFT LINK TO DOC:
  fname="example"number".dat";
  print $0 >> fname;
# Find replicas to use
  if(match($1,"SETTINGS")) {
     for(i=1;i<=NF;++i) {
         if( match($i,"NREPLICAS=") ) { nreplicas=$i; sub("NREPLICAS=","",nreplicas); }
     }
  }
}
