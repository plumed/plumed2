{
  if(match($0,"BEGIN_PLUMED_FILE")){
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
    inside=0; close("example.dat"); 
    if( nreplicas==1 ) { system("plumed gen_example --plumed example.dat --out example.html --name eg"number" --status " status "> /dev/null"); } 
    else { system("mpirun -np " nreplicas " plumed gen_example --plumed example.dat --out example.html --name eg"number" --status " status " --multi " nreplicas " > /dev/null");}
    system("cat example.html"); 
    system("rm example.dat example.html");
    sub("</pre>","");
    print;
    next;
  }
  if(!inside){
    print;
    next;
  }

# DRAFT LINK TO DOC:
  print $0 >> "example.dat";
# Find replicas to use
  if(match($1,"SETTINGS")) {
     for(i=1;i<=NF;++i) {
         if( match($i,"NREPLICAS=") ) { nreplicas=$i; sub("NREPLICAS=","",nreplicas); }
     }
  }
}
