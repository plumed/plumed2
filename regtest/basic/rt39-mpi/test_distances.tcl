mol new diala_traj_nm.xyz  type xyz first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set nn1 [  molinfo top get numframes ] 
mol new all.pdb  type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
set nn2 [  molinfo 1 get numframes ] 
set sel1 "index 0 4 5 6 7 8 9 10 14 15 16 17 18  "
set sel2 " all   "
#for { set i 0 } { $i < $nn1 } { incr i } {
for { set i 0 } { $i < 1 } { incr i } {
puts "FRAME $i"
 #scale all traj 
 set runsel_all [  atomselect 0 "all" frame $i ]
 set c  [ $runsel_all get { x y z } ]
 set ci [ $runsel_all get index ]
 for  { set j 0 } { $j < [llength $c ]} { incr j } {
         set idx [ lindex $ci $j  ]
         set aa [ atomselect 0 "index $idx"  frame $i ]
         $aa set x [ expr [ lindex [ lindex $c $j ] 0 ] *10. ]
         $aa set y [ expr [ lindex [ lindex $c $j ] 1 ] *10. ]
         $aa set z [ expr [ lindex [ lindex $c $j ] 2 ] *10. ]
 } 
 set runsel [ atomselect 0 $sel1 frame $i    ]
 for { set j 0 } { $j < $nn2  } { incr j } {
	set refsel [ atomselect 1 $sel2 frame $j   ]
        set mymat  [ measure fit  $runsel $refsel ]
        [ atomselect 1 "all " ]  move $mymat
        puts  "[expr pow([measure rmsd $runsel $refsel],1)/100.]" 
 } 
}
puts "END"
