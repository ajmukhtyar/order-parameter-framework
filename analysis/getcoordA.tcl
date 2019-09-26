set k SSSS
animate goto 1
[atomselect top "name 1"] writelammpstrj A$k.lammpstrj
#set nfr [molinfo top get numframes]
set sel [atomselect top all]
set slist [$sel get index]
#animate goto $nfr
set cutf 3.0
[atomselect top "name 1"] set type 1
[atomselect top "name 2"] set type 2
#set ii SSSS
#set jj [expr $ii + 1]
#for {set k RRRR} {$k < $jj} {incr k TTTT} {
#animate goto $k
####Get coordinates of all particles######
#set fid [open "coord$k.xyz" w]
#set fid1 [open "neigh$k.txt" w]
set fid2 [open "atomtype$k.txt" w]
set fid3 [open "sameneigh$k.txt" w]
foreach i $slist {
set sel2 [atomselect top "index $i"]
set coord [$sel2 get {x y z}]
#puts $fid "$i [lindex $coord 0]"

#####Get neighbor atoms#####
#set sel3 [atomselect top "pbwithin $cutf of index $i"]
#set seli [$sel3 get index]
#puts $fid1 "$i $seli"
#$sel3 delete

###Get atom type####
set seln [$sel2 get type]
puts $fid2 "$i $seln"
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
puts $fid3 "$i $selni"  
}

#close $fid
#close $fid1
close $fid2
close $fid3
#}

 
