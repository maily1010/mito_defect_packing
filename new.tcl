# Uses VMD version 1.9.2 
# Written by Ky Wildermuth (May 2016)
# updated (Aug 27 2016): add protein min/max datafile for cutoff analysis. add PBC dimension datafile for granger timeseries
# updated (Oct 17 2016): configured such that protein images do not include membrane elements
# updated (Jan 01 2017): New wrapping command seems to work better (-centersel).

# NOTES AND DIRECTIONS:
# Script assumes membrane is centered at origin (z-coord)
# Modify dcd and psf to be read, and home directory that contains these files
# Modify lipid selection in VMD representation section of script to species in your system (default is DOPC and DOPS)
# Modify representation type ('modstyle' variable. default is 'quicksurf' but could use 'VDW') -- quicksurf is very slow if you have many images. VDW spheres might be sufficient, but no tests done yet (April 2017) 
# If dcd has wrapping issue (bonds stretched when viewed in VMD), run the pbc wrap command (default commented out)
# Modify or remove coordinate data collection section (default: activated). Used in analysis with associated Python scripts
# Can optionally make a video for qualitative viewing of peptide binding from side view -- see loop labeled "side" (default commented out)
# Script is designed for use of scavenger partition on DT1: will not write image files if they already exist!
# Script will build output storage architecture in the home directory that is used by associated Python analysis scripts

#######################
# SET USER PARAMETERS #
#######################
# Select the x and y dimensions of the images (NOTE: recommended to maintain square ratio)
set xsize 600
set ysize 600
# Select main working directory (directory containing psf and dcd files)
set home "/data/mai/MITO/PK/M1"
# move to the working directory
cd $home
# Select names of DCD and PSF files and load the molecule and all frames; configure 'step' parameter to be frames skipped per interval
mol delete all
mol new {coor_files/coor_1900.dcd} type {dcd} first 0 last -1 step 1 waitfor -1
mol addfile {coor_files/step5_assembly.psf} type {psf} first 0 last -1 step 1 waitfor -1 0
# define variable for common lipid polar motif
set polar_motif "P O13 O14 O12 O11 C1 HA HB C2 HS O21 C21 O22 C22 H2R H2S C3 HX HY O31 C31 O32 C32 H2X H2Y"
# List species present and associated atoms in polar region 
# if your lipid is not here, add to the list: Polar considered 2nd carbon atom and above
# Must also edit the species called in the "CREATE REPRESENTATION" section of the script below
set polar_PS "N HN1 HN2 HN3 C12 H12A C13 O13A O13B C11 H11A H11B $polar_motif"
set polar_PI "C12 H2 O2 HO2 C13 H3 O3 HO3 C14 H4 O4 HO4 C15 H5 O5 HO5 C16 H6 O6 HO6 C11 H1 $polar_motif"
set polar_PA "H12 $polar_motif"
set polar_PC "N C12 C13 C14 C15 H12A H12B H13A H13B H13C H14A H14B H14C H15A H15B H15C C11 H11A H11B $polar_motif"
set polar_PE "N HN1 HN2 HN3 C12 H12A H12B C11 H11A H11B $polar_motif"
##########################

# make sub_directories
cd $home
file mkdir results
# set results "/results"
cd results
file mkdir bottom top data_files side
file mkdir [file join [pwd] top/top_lipid] [file join [pwd] top/top_prot]
file mkdir [file join [pwd] bottom/bot_lipid] [file join [pwd] bottom/bot_prot]

# move working directory to data_files
cd data_files
# set filenames to write data to
set file1 [open "pbc_scale_dimensions.dat" w]
set file2 [open "top_phos_coord.dat" w]
set file3 [open "bot_phos_coord.dat" w]
set file4 [open "custom_coord.dat" w]
#set file5 [open "prot_coord.dat" w]
#set file6 [open "prot_minmax_coord.dat" w]
set file7 [open "PBC_coord.dat" w]

# move back to main directory
cd $home
cd results
# subdirectories to output images to
set top_lipid [file join [pwd] top/top_lipid]
set bot_lipid [file join [pwd] bottom/bot_lipid]
#set top_prot [file join [pwd] top/top_prot]
#set bot_prot [file join [pwd] bottom/bot_prot]
set side [file join [pwd] side]

# SETUP GENERAL IMAGE DISPLAY
# select lights
light 0 on
light 1 on
light 2 off
light 3 off
# turn on orthographic view
display projection Orthographic
# turn off axes in corner
axes location Off
# minimize near clip
display nearclip set 0.010000
# reset view angle 
display resetview
package require pbctools
# if protein (or any other molecule) is wrapping around such that bonds are stretched, run the following line
# pbc wrap -centersel "protein" -center com -compound residue -all
# Find number of frames
set frames [molinfo top get numframes]

# CREATE REPRESENTATIONS
# POLAR REP (blue)
molinfo top get numreps
mol modselect 0 0 (resname DOPE and name $polar_PE) or (resname POPC and name $polar_PC) 
mol modcolor 0 0 ColorID 0
mol modstyle 0 0 QuickSurf 1.000000 0.500000 1.000000 1.000000
mol showperiodic 0 0 xyXY
mol numperiodic 0 0 1
mol modmaterial 0 0 Goodsell
# NONPOLAR REP (yellow)
mol addrep 0
mol modselect 1 0 (resname DOPE and not name $polar_PE) or (resname POPC and not name $polar_PC) 
mol modcolor 1 0 ColorID 4
mol modstyle 1 0 QuickSurf 1.000000 0.500000 1.000000 1.000000
mol showperiodic 0 1 xyXY
mol numperiodic 0 1 1
mol modmaterial 1 0 Goodsell
# HMMM REP (yellow) 
mol addrep 0
mol modselect 2 0 resname DCLE and z < 20 and z > -20
mol modcolor 2 0 ColorID 4
mol showperiodic 0 2 xyXY
mol numperiodic 0 2 1
mol selupdate 2 0 1
mol modstyle 2 0 QuickSurf 1.000000 0.500000 1.000000 1.000000
mol modmaterial 2 0 Goodsell
# PROTEIN REP (red) with no PBC so that only in the main cell
mol addrep 0
mol modcolor 3 0 ColorID 1
mol modselect 3 0 segname PROA
mol modstyle 3 0 QuickSurf 1.000000 0.500000 1.000000 1.000000
mol modmaterial 3 0 Goodsell
# see output file for list of reps to check
mol list

# COLLECT COORDINATE DATA OF PROTEIN AND PHOSPHATES OF BILAYER
# select all phosphates on top leaflet of bilayer
set phos_sel_top [atomselect top "name P and z > 40"]
# select all phosphates on bottom leaflet of bilayer
set phos_sel_bot [atomselect top "name P and z < 40"]
# OPTIONAL: Custom residue selection (default is alpha carbon of TRP) 
set custom_sel [atomselect top "resname TRP and name CA"]
# select all protein residues for additional data collection
#set prot_sel [atomselect top "protein"]
# loop over all frames and make selections defiend above
for {set i 0 } {$i < $frames } {incr i } {
    $custom_sel frame $i
    $phos_sel_top frame $i
    $phos_sel_bot frame $i
#    $prot_sel frame $i
    # collect coordinate informatin for selections defined above
    set custom_coord [$custom_sel get {x y z}]
    set phos_com_top [measure center $phos_sel_top]
    set phos_com_bot [measure center $phos_sel_bot]
#    set prot_com [measure center $prot_sel]
#    set prot_minmax [measure minmax $prot_sel]
    set PBC_cell [pbc get -now]
    # put collected coordinates in respective data files
    puts $file2 "$i $phos_com_top"
    puts $file3 "$i $phos_com_bot"
    puts $file4 "$i $custom_coord"
#    puts $file5 "$i $prot_com"
#    puts $file6 "$i $prot_minmax"
    puts $file7 "$i $PBC_cell"
}
# close the coordinate file
close $file2
close $file3
close $file4
#close $file5
#close $file6
close $file7

# PREPARE AND RENDER IMAGES
#
# Not recommended to change this variable. If you want to change interval of frames to be read from DCD, change when importing dcd.
set int 1
# reset view angle 
display resetview
#
# TOP LIPID
# draw a green pbc box around COM of protein
#pbc box -center com -centersel protein -shiftcenter {0 0 60} -color green
pbc box -center origin -shiftcenter {60 60 40} -color green
# measure dimensions of pbc box of first frame (to be used for scaling in analysis)
#pbc set {134.466293 134.466293 83.306702} -all
set cell [pbc get -first 0 -last 0]
# record pbc dimensions to data file and close the file
puts $file1 "$cell"
close $file1
# turn on membrane reps, turn of protein reps
mol showrep 0 0 on
mol showrep 0 1 on
mol showrep 0 2 on
mol showrep 0 3 off
# move to the first frame 
set frame 0
# move to top_lipid directory 
cd $top_lipid
# loop through the frames
for {set i 0} {$i < $frames} {incr i $int} {
    # go to the given frame
    animate goto $i
    # take the picture
	set filename snap.[format "%04d" $i]
    set flag [file exists $filename.bmp]
    if {$flag == 0} {
        pwd
	render Tachyon $filename "/ibpc/vrubel/mai/VMD/vmd-1.9.3/lib/tachyon/tachyon_LINUXAMD64" $filename -lowshade -res 600 600 -format BMP -o $filename.bmp
	file delete snap.[format "%04d" $i] 
    }
}
#
# TOP PROTEIN
# turn on protein rep
#mol showrep 0 0 off
#mol showrep 0 1 off
#mol showrep 0 2 off
#mol showrep 0 3 on
#cd $top_prot
# loop through the frames
#for {set i 0} {$i < $frames} {incr i $int} {
    # go to the given frame
#    animate goto $i
    # take the picture
#    set filename snap.[format "%04d" $i]
#    set flag [file exists $filename.bmp]
#    if {$flag == 0} {
#        pwd
#	render Tachyon $filename "/ibpc/vrubel/mai/VMD/vmd-1.9.3/lib/tachyon/tachyon_LINUXAMD64" $filename -lowshade -res 600 600 -format BMP -o $filename.bmp
#	file delete snap.[format "%04d" $i] 
#    }
#}
#
# BOTTOM PROTEIN
#rotate x by 180
#pbc box -center com -centersel protein -shiftcenter {0 0 -60} -color green
#cd $bot_prot
# loop through the frames
#for {set i 0} {$i < $frames} {incr i $int} {
    # go to the given frame
#    animate goto $i
    # take the picture
#    set filename snap.[format "%04d" $i]
#    set flag [file exists $filename.bmp]
#    if {$flag == 0} {
#        pwd
#	render Tachyon $filename "/ibpc/vrubel/mai/VMD/vmd-1.9.3/lib/tachyon/tachyon_LINUXAMD64" $filename -lowshade -res 600 600 -format BMP -o $filename.bmp 
#    	file delete snap.[format "%04d" $i] 
#	}
#}
#
# BOTTOM LIPID
mol showrep 0 0 on
mol showrep 0 1 on
mol showrep 0 2 on
mol showrep 0 3 off
cd $bot_lipid
# loop through the frames
for {set i 0} {$i < $frames} {incr i $int} {
    # go to the given frame
    animate goto $i
    # take the picture
    set filename snap.[format "%04d" $i]
    set flag [file exists $filename.bmp]
    if {$flag == 0} {
        pwd
	render Tachyon $filename "/ibpc/vrubel/mai/VMD/vmd-1.9.3/lib/tachyon/tachyon_LINUXAMD64" $filename -lowshade -res 600 600 -format BMP -o $filename.bmp 
   	file delete snap.[format "%04d" $i] 
	}
}
