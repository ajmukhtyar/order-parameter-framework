set fid1 [open "clus.txt" r]
set count 0
while {[gets $fid1 line]>=0} {
set sel [atomselect top "index $count"]
set val [lindex $line 0]
$sel set beta $val
set count [expr $count + 1] 
}
