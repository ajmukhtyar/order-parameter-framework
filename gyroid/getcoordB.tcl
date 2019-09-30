[atomselect top "name 2"] writelammpstrj B.lammpstrj
set sel [atomselect top "name 2"]
set slist [$sel get index]
set cutf 3.0
[atomselect top "name 1"] set type 1
[atomselect top "name 2"] set type 2
####Get coordinates of all particles######
set fid2 [open "B_atomtype.txt" w]
set fid3 [open "B_sameneigh.txt" w]
foreach i $slist {
set sel2 [atomselect top "index $i"]
set coord [$sel2 get {x y z}]
set j [expr ($i-6912)]

###Get atom type####
set seln [$sel2 get type]
puts $fid2 "$j $seln"
$sel2 delete

###Get same neighbors###
if ($seln==1) {
set selA [atomselect top "name 1 and not index $i and pbwithin $cutf of index $i"]
set selni [$selA get index]
$selA delete
}
if ($seln==2) {
set selB [atomselect top "name 2 and not index $i and pbwithin $cutf of index $i"]
set selni [$selB get index]
$selB delete
}
puts $fid3 "$j $selni"  
}

close $fid2
close $fid3

 
