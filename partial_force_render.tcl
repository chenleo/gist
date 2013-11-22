##
## Partial Force Render 36
##  
## Copy Not Right Chen Chen B040-DKFZ
## 
## A script to render trajectory according to scalar force calculate from partial_force
## 
## Major way to use:
## Step 1: 	mol load gro *.gro xtc *.xtc 	#loaded *.gro and correspond *.xtc into vmd
## Step 2.1: 	package require partial_force_render
## Step 2.2:	partial_force_render "*.csv"	#use .csv file to render the whole structure
##
## Display:
## 	Color map correspond force value, together with a scale bar.
##
## Note:
## *.csv file is the output from partial_force calculatiion program. 
##	The first line is the head (time, residues)
##	Start from the second line, is the force value use to render trajcetories.
##	Notice the script only checks the number of residue corrspond to main molecule.
##	User should take the responsibility to load the right file from certain calculation.
##

package provide partial_force_render 36

proc partial_force_render { csv_file } {

	#Initialize:
	puts "Start to load csv file......"
	
	#set color to user
	mol modcolor 0 top User
	mol colupdate 0 top 1

	#Set force color to RGB
	puts "set force color to RGB"
	color scale method "RGB"
	
	#Get number of frames to render
	puts "set frame_number to total frames"
	set frame_number [molinfo top get numframes]

	#Select top atoms
	puts "select atoms to top all"
	set sel [atomselect top all]

	#Now open the .csv data file
	puts "open file"
	set force_info [open $csv_file r]

	#set color max and min
	puts "set color min and max"
	set color_max 1000	
	#Largest force 1000, above will count into max
	
	set color_min 0.0	
	#Of course

	#read out head
	puts "get rid of head"
	gets $force_info head_info

	#Now loop over the frames, starts from 1 because 0 is the base.
	puts "looping frames"
	for {set current_frame 1} {$current_frame <= $frame_number} {incr current_frame} {
		
		#set current_time [lindex $force_value 0]
		#puts "rednering for frame $current_frame, for time $current_time"

		#real start to read in force_information
		gets $force_info frame_force

		#according to csv format, seperate frame_force by "," then put it to force_value
		set force_value [split $frame_force ","]
		
		set current_time [lindex $force_value 0]
	        puts "rednering for frame $current_frame, for time $current_time"

		#Then loop through force_value for each residue
		for {set id_res 1}  { $id_res < [llength $force_value]} {incr id_res} {
		
			set current_res [atomselect top "resid $id_res"]	
			#since id_res start from 1
			
			$current_res frame $current_frame
			$current_res set user [lindex $force_value $id_res]
			$current_res delete
		}	
		
	}

	puts "Passed"
	mol scaleminmax top 0 $color_min $color_max
	puts "should showed"

	close $force_info 
	puts "Done!"
}
