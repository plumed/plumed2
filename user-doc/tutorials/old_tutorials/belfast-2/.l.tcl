set name "wrong_path.dat" ; puts -nonewline  ""
mol new $name type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
#
# reset com
#
proc reset_com { molid id } {
        set selocc [ atomselect $molid  "occupancy > 0." frame $id  ] ; puts -nonewline  ""
        set occval  [ $selocc get occupancy  ] ; puts -nonewline  "" 
        set occindex  [ $selocc get index  ] ; puts -nonewline  ""
	set gc [veczero]
    	set k 0
    	set totw 0. 
        foreach coord [ $selocc  get {x y z}] {
    		set myw [ lindex $occval $k ] 
    	        set gc [vecadd $gc [ vecscale $myw $coord ] ] 
    		set totw [ expr $totw+ $myw ]
    		incr k
        }
    	set gc  [vecscale [expr 1.0 /$totw ] $gc]
    	set newcoord []
        foreach coord [ [ atomselect $molid all  frame $id ] get {x y z}] {
    	    lappend newcoord [ vecsub $coord $gc ]	
    	}
    	[  atomselect $molid all frame $id ] set  { x y z } $newcoord 

}

proc calc_msd_and_dist { molid id1 id2 } {

	set plumedver 1 
	#puts "calculating msd and dist vector ... "
	#
	# save coordinates	
	#
        set store_id1 [ [ atomselect $molid  all frame $id1 ]  get { x y z } ] 
        set store_id2 [ [ atomselect $molid  all frame $id2 ]  get { x y z } ] 
	# 
	# center the com 
	# 
	#puts "reset com of $id1 and $id2 "	
	reset_com $molid $id1
	reset_com $molid $id2
	#puts "rescale the coordinate according the alignment weights "
        #puts "rescaling values ..."
	# save before the expansion but before reset of com
        set allcoor_id1 [ [ atomselect $molid  all frame $id1 ]  get { x y z } ] 
        set allcoor_id2 [ [ atomselect $molid  all frame $id2 ]  get { x y z } ] 
        #
        set selocc [ atomselect $molid  "occupancy > 0."  ] ; puts -nonewline  ""
        set occval  [ $selocc get occupancy  ] ; puts -nonewline  "" 
        set occindex  [ $selocc get index  ] ; puts -nonewline  ""
        set coord [ [  atomselect $molid "index $occindex" frame $id1 ] get { x y z } ] 
        for { set j 0 } { $j <  [llength $coord ] } { incr j } { 
                set myind [ lindex $occindex $j ]  ; puts -nonewline  ""
		if { $plumedver == 2 } { 
                	[  atomselect $molid "index $myind" frame $id1 ] set  { x y z }  [ list [ vecscale [ expr sqrt( [ lindex $occval $j ] ) ] [ lindex $coord $j ] ]    ] 
		} elseif { $plumedver == 1 } { 
                	[  atomselect $molid "index $myind" frame $id1 ] set  { x y z }  [ list [ vecscale  [ lindex $occval $j ]  [ lindex $coord $j ] ]    ] 
		}
        
        }
        set coord [ [  atomselect $molid "index $occindex" frame $id2 ] get { x y z } ] 
        for { set j 0 } { $j <  [llength $coord ] } { incr j } { 
                set myind [ lindex $occindex $j ]  ; puts -nonewline  ""
   	        if { $plumedver == 2 } {
             	   [  atomselect $molid "index $myind" frame $id2 ] set  { x y z }  [ list [ vecscale [ expr sqrt ( [ lindex $occval $j ] ) ]  [ lindex $coord $j ] ]    ] 
   	        } elseif { $plumedver == 1 } {
                   [  atomselect $molid "index $myind" frame $id2 ] set  { x y z }  [ list [ vecscale  [ lindex $occval $j ]   [ lindex $coord $j ] ]    ] 
		}
        }
       # now that it is scaled is time to calculate the alignment
       set mymat  [ measure fit [ atomselect $molid "index $occindex" frame $id2 ] [  atomselect $molid "index $occindex" frame $id1 ]  ]	
       # restore the value post resetcom	
       [ atomselect $molid all frame $id1 ] set  { x y z } [ lindex $allcoor_id1  ] ; puts -nonewline  ""
       [ atomselect $molid all frame $id2 ] set  { x y z } [ lindex $allcoor_id2  ] ; puts -nonewline  ""
       # now move id2 onto id1 
       [ atomselect $molid  all frame $id2 ] move $mymat 
       # now calculate the weighted msd
       set selbeta [ atomselect top  "beta > 0."  ] ; puts -nonewline  ""
       set betaval [ $selbeta get beta  ] ; puts -nonewline  ""
       set betaindex  [ $selbeta get index  ] ; puts -nonewline  ""
       set totw 0. 
       set totv 0. 
       for { set k 0 } { $k < [ llength $betaindex ] } { incr k } {
             set myind [ lindex $betaindex $k] ; puts -nonewline  ""	
             set ww [ lindex $betaval $k ]	 ; puts -nonewline  ""
             set pi [ lindex [ [  atomselect $molid "index $myind" frame $id2 ] get { x y z }   ] 0 ] ; puts -nonewline  ""
             set pj [ lindex [ [  atomselect $molid "index $myind" frame $id1 ] get { x y z }   ] 0 ]		 
	     if { $plumedver == 1 } {
	             set v [ expr [  veclength2 [ vecsub $pi $pj ] ] *  $ww * $ww ]  ; puts -nonewline  ""	
        	     set totw [ expr $totw + $ww*$ww ]  ; puts -nonewline  ""
             	     set totv [ expr $totv + [expr $v ] ]  ; puts -nonewline  ""
	     } elseif { $plumedver == 2 } {
	        set v [ expr [  veclength2 [ vecsub $pi $pj ] ]]  ; puts -nonewline  ""	
        	set totv [ expr $totv + ($v)*$ww  ]  ; puts -nonewline  ""
       	        set totw [ expr $totw + $ww ]  ; puts -nonewline  ""
	     }
             unset pi; unset pj ; unset v ; unset ww ; unset myind
       } 
       if { $plumedver == 1 } {
         # a-la-plumed 1
         set totv [ expr sqrt($totv / $totw) ]		
       } elseif { $plumedver == 2 } {
         set totv [ expr sqrt($totv / $totw) ]		
       }
       #puts "MSD between $id2  and $id1 : $totv "
       # 
       # give back also id1 and id2-id1 
       #
       set pi [ [  atomselect $molid all frame $id2 ] get { x y z }  ] ; puts -nonewline  "" 
       set pj [ [  atomselect $molid all frame $id1 ] get { x y z }  ]	; puts -nonewline  ""	 
       set delta [] ; puts -nonewline  ""
       for  { set k 0 } { $k < [ llength $pj ] } { incr k } {	
	 	lappend delta [ vecsub [ lindex $pi $k ] [ lindex $pj $k ] ]	; puts -nonewline  ""	
       }	

	return [ list  $totv $pj $delta  ] 
}
proc eval_lambda {  molid    } {
	set nframes [ molinfo top get numframes ]
	set delta [] ;
	for { set i  0  } { $i < $nframes  } { incr i } {
			# calculate the distance
	        	set rr [ calc_msd_and_dist $molid $i [expr  $i +1 ] ] ; puts -nonewline  "" 
			#puts "DELTA from $i to $j :  [lindex $rr 0 ] "	
			lappend delta  [lindex $rr 0 ]
	}

	set avg 0. ; 
	set devstd 0. ; 
	for { set i 0 } { $i < [ llength $delta ]} { incr i } {
		set avg [ expr $avg + [lindex $delta  $i ] **2  ]  ;
		set devstd [ expr $devstd + [lindex $delta  $i ] **4  ]  ;
	}
	set avg [ expr $avg / [ llength $delta ] ] 
	set devstd [ expr $devstd / [ llength $delta ] ] 
	set devstd [  expr sqrt( $devstd -$avg*$avg ) ] 
	puts "******************************"
	puts ">>>> MSD        $avg  Ang^2 " 
	puts ">>>> MSD_STDEV  $devstd  Ang^2 " 
	set safe [ expr 2.3/( $devstd+$avg  ) ] ;puts -nonewline "" 
	set lambda [ expr  2.3/$avg  ]  ;puts -nonewline "" 
	puts ">>>> LAMBDA $lambda  \[Ang ^ -2\] "
	puts ">>>> SAFE_LIMIT $safe \[Ang ^ -2\] "
	puts ">>>> Note: if you use gromacs you have to increase these values of 100 since units are \[ nm ^-2 \]"
	puts ">>>> Note2: LAMBDA and SAFE_LIMIT are truly different your matrix could be messed up: check MSD_MATRIX FILE "
	puts "******************************"
}

proc eval_matrix {  molid fname  } {
	set bak1 0
	set bak2 1 
	set nframes [ molinfo top get numframes ]
	for { set i  0  } { $i < $nframes  } { incr i } {
		puts "MATRIX_CALCULATION element $i "	
		for { set j $i } { $j < $nframes  } { incr j } {
			# calculate the distance
	        	set rr [ calc_msd_and_dist $molid $i $j ] ; puts -nonewline  "" 
			#puts "DELTA from $i to $j :  [lindex $rr 0 ] "	
			set arr($i,$j) [lindex $rr 0 ]
			set arr($j,$i) [lindex $rr 0 ]
			# set back the coordinates
		}
	}
   	set fp [ open $fname w ] ; puts -nonewline  ""
	for { set i  0  } { $i < $nframes  } { incr i } {
		for { set j 0 } { $j < $nframes } { incr j } {
		   puts $fp "[expr $i+1] [expr $j+1]  $arr($i,$j)"
		}
   		puts $fp  "  "
   	}
	close $fp
}
proc eval_diag {  molid  fname } {
	set nframes [ molinfo top get numframes ]
	for { set i  0  } { $i < [ expr $nframes -1 ]  } { incr i } {
		puts "DIAGONAL_CALCULATION element $i "	
		# calculate the distance
	        set rr [ calc_msd_and_dist $molid $i [expr  $i +1 ] ] ; puts -nonewline  "" 
		#puts "DELTA from $i to $j :  [lindex $rr 0 ] "	
		set arr($i) [lindex $rr 0 ]
	}
   	
   	set fp [ open $fname w ] ; puts -nonewline  ""
	for { set i  0  } { $i < [ expr $nframes -1 ] } { incr i } {
	   puts $fp "[expr $i+1]   $arr($i)"
   	}
	close $fp
	
}
eval_matrix  0  "RMSD_MATRIX" ;  puts -nonewline  ""
eval_diag  0  "RMSD_DIAGONAL"  ;  puts -nonewline  ""
eval_lambda  0   ;  puts -nonewline  ""
puts "BYE..."
quit
