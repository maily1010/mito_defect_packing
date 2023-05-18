# MEMBRANE PACKING DEFECT ANALYSIS 
# WRITTEN FOR: PYTHON (v 2.7.8) and openCV (v 3.1.0) 
# USES MODULES: numpy, os, cv2
# CREATED BY: KY WILDERMUTH (05/2016)
#
# This script will conduct various packing defect analysis on a single leaflet of a lipid membrane.
# It will produce three output files that provide the following data for each frame in the image set: 
# (1) 	"local_defect.dat" provides the total normal area of defect exposure directly below the protein, 
# 		the normal area of the protein, and a comparison value that describes the magnitude of difference
#		between the two defect shapes.
# (2)	"global_defect.dat" provides the total normal area of the entire leaflet followed by the total normal
#		area of defect exposure
# (3)	for each frame in the image set, "defect_freq_area.dat" provides a list of the normal area of all individual 
#		defects as well as a comparison value that describes the magnitude of difference between the selected 
#		defect and the protein shape.
#
# To run this script, an appropriate image set must be prepared. 
# See the associated tcl script for outputting packing defect images from VMD
#
# All inputs that must be defined by the user are found in Step 1 below. I recommend visualizing various 
# steps using the cv2.imshow() command (see openCV documentation) when using a new image set to check for expected results
#
# Area measurements in this script are converted to Angstroms from pixels via an in-program conversion factor.
# The conversion factor is based on a VMD measurement of the PBC box dimensions in the first frame 
# of the image set and is assumed to be in Angstroms.
# The conversion factor is reported in the "pbc_scale_dimensions.dat" file
#
# This script can be used on a membrane simulation that does not have a protein. If this is the case, ignore
# the output file "local_defect.dat"
#
 
import numpy as np 
import cv2
import os, os.path
 
###############################
### STEP 1: SET USER INPUTS ###
###############################
 
# set a limit on defect selection based on area of contour (removes very small defects)
def is_contour_bad(c):
	# calculate the contour area in A^2 
	area = cv2.contourArea(c)*conv
	# the contour is "bad" if its area is smaller than 1 A^2
	return not area > 1
 
# Define the color range used for the masks of each component in the system
# NOTE: cv2 uses BGR color space 
# 1 is lipid core, 2 is protein, 3 is pbc box (default color ranges are for yellow, red, and green respectively)
boundary_1 = [([0, 200, 200], [25, 255, 255])]
boundary_2 = [([0, 0, 200], [25, 25, 255])]
boundary_3 = [([0, 150, 0], [60, 255, 60])]
 
parent_dir = "/ibpc/vrubel/mai/Desktop/Try10"
 
conv_path = os.path.join(parent_dir, "data_files/pbc_scale_dimensions.dat")
local_path = os.path.join(parent_dir, "data_files/local_defect_bot.dat")
global_path = os.path.join(parent_dir, "data_files/global_defect_bot.dat")
freq_path = os.path.join(parent_dir, "data_files/defect_freq_area_bot.dat")
 
# select image input directories:
lipid_path = os.path.join(parent_dir, "bottom/bot_lipid")
#prot_path = os.path.join(parent_dir, "bottom/bot_prot")
 
# enter x and y dimensions of the pbc box (in A^2) from the first frame used in the image set (obtained from VMD)
box_x, box_y = 134.466293, 134.466293
 
 
################################
### STEP 2: PROTEIN ANALYSIS ###
################################
 
# start the frame number counter
count = 0
#iterate over all *.bmp files in the selected directories
lst = os.listdir(lipid_path)
srt = sorted(lst)
for i in srt:
	count = count + 1
	print('%d' % count)
	#os.chdir(prot_path)
	#if i.endswith(".bmp"):
		# import protein image
	#	prot_image = cv2.imread(i)
		# loop over the boundary color range for the protein
	#	for (lower, upper) in boundary_2:
			# create NumPy arrays from the boundaries
	#		lower = np.array(lower, dtype = "uint8")
	#		upper = np.array(upper, dtype = "uint8")
		# apply the protein mask based on the color range
	#	prot_select = cv2.inRange(prot_image, lower, upper)
		# define a contour for the protein
	#	(_,prot_c, _) = cv2.findContours(prot_select, cv2.RETR_EXTERNAL,
	#	cv2.CHAIN_APPROX_SIMPLE)
		# draw the contour on the original mask and fill it in
	#	prot_mask = cv2.drawContours(prot_select, prot_c,-1,255,-1)
 
################################
### STEP 3: DEFINE UNIT CELL ###
################################
 
	# import lipid image
	os.chdir(lipid_path)
	if i.endswith(".bmp"):
		image_in = cv2.imread(i)
		# add filter with edge conservation to smooth masks (run only one of the following filters - try MSF if BF is not sufficient)
		#image = image_in 	# enable this if you you don't need to use a filter 
		#image = cv2.pyrMeanShiftFiltering(image_in, 21, 51)
		image = cv2.bilateralFilter(image_in,10,75,75)
 
		# for filter optimization, check the original and filtered images side by side:
		#cv2.imshow("", np.hstack([image, image_in]))
		#cv2.waitKey(0)
 
		# loop over the boundary color range of the pbc box
		for (lower, upper) in boundary_3:
			# create NumPy arrays from the boundaries
			lower = np.array(lower, dtype = "uint8")
			upper = np.array(upper, dtype = "uint8")	 
		# apply the pbc box mask based on the color range
		pbc_select = cv2.inRange(image, lower, upper)
		# select the contour corresponding to the unit cell box
		(_,box_c) = cv2.findContours(pbc_select, cv2.RETR_EXTERNAL,
		cv2.CHAIN_APPROX_SIMPLE)
		# draw the contour on the original mask and fill it in
		box_mask = cv2.drawContours(pbc_select, box_c,-1,255,-1)
		box = cv2.bitwise_and(image, image, mask = box_mask)
 
		# use the area of the pbc box in the first image as a reference for conversion
		if count < 2:
			# define area of first frame pbc box in A^2 from VMD data
			box_ang = box_x*box_y
			print(box_ang)
			#define area of first frame pbc box in picels from mask size
			box_pix = np.count_nonzero(box_mask)
			print(box_pix)
			# define the conversion ratio from pixels to A^2
			conv = box_ang / box_pix
			# remove existing output files if they exist
			try:
				os.remove(local_path)
				os.remove(global_path)
				os.remove(freq_path)
			except OSError:
				pass
			# add data to output file
			f0 = open(conv_path, 'a')
			f0.write('conversion factor is %f A^2/pixel. All values are reported in A^2 \n' % conv)
			f1 = open(local_path, 'a')
			f1.write('frame_# area_prot_defect area_protein comparison_value \n')
			f2 = open(freq_path, 'a')
			f2.write('frame_# area_defect comparison_value \n')
			f3 = open(global_path, 'a')
			f3.write('frame_# area_box area_total_defect \n')
 
#################################
### STEP 4: MEMBRANE ANALYSIS ###
#################################
 
		box_pix = np.count_nonzero(box_mask)
		box_ang = conv*box_pix
		print('area of box %d is %d square angstrom' % (count, box_ang))
		f3 = open(global_path, 'a')
		f3.write('%d %f ' % (count, box_ang))
 
		# loop over the boundary color range for defects
		for (lower, upper) in boundary_1:
			# create NumPy arrays from the boundaries
			lower = np.array(lower, dtype = "uint8")
			upper = np.array(upper, dtype = "uint8")
		# apply the mask based on color range
		mask = cv2.inRange(box, lower, upper)
		output_defect = cv2.bitwise_and(box, box, mask = mask)
		#cv2.imshow("defects", output_defect)
		#cv2.waitKey(0)
 
		# compute total area of defect mask and convert to square angstroms
		defectt_area_pix = np.count_nonzero(mask)
		defectt_area_ang = defectt_area_pix*conv
		print('area of defects: %d square angstrom' % defectt_area_ang)
		# add data to defect area file
		f3.write('%f \n' % defectt_area_ang)
 
###################################
### STEP 5: REMOVE BAD CONTOURS ###
###################################
 
		# Prepare image for contour analysis
		# convert defect mask to grayscale
		gray = cv2.cvtColor(output_defect, cv2.COLOR_BGR2GRAY)
		#add threshold to image
		(T, thresh) = cv2.threshold(gray, 5, 255, cv2.THRESH_BINARY)
		# find all contours
		(_,cnts_all, _) = cv2.findContours(thresh, cv2.RETR_EXTERNAL,
			cv2.CHAIN_APPROX_SIMPLE)
		# find all bad contours
		mask = np.ones(image.shape[:2], dtype="uint8")*255
		# loop over all contours
		for c in cnts_all:
			# if the contour is bad, draw it on the mask
			if is_contour_bad(c):
				cv2.drawContours(mask, [c],-1,0,-1)
		# remove bad contours from the image
		good_defects = cv2.bitwise_and(thresh, thresh, mask=mask)
		# find the remaining good contours
		(_,cnts_good, _) = cv2.findContours(good_defects, cv2.RETR_EXTERNAL,
			cv2.CHAIN_APPROX_SIMPLE)
		print ("There are a total of %d selected contours" % (len(cnts_good)))
 
###################################
### STEP 6: DEFECT ANALYSIS #######
###################################
 
		prot_defect = np.ones(image.shape[:2], dtype="uint8")*0
		# loop over the countours and compute the following
		for c in cnts_good:
			#define the single contour object mask
			mask = np.ones(image.shape[:2], dtype="uint8")*0
			single_contour = cv2.drawContours(mask, [c],-1,255,-1)
			#cv2.imshow("single contour", single_contour)
			#cv2.waitKey(0)
			prot_overlap = cv2.bitwise_and(prot_mask, prot_mask, mask = single_contour)
			prot_defect = prot_overlap + prot_defect
			# compare protein contour to defect contour
			compare = cv2.matchShapes(single_contour, prot_mask, 1, 0)
			# calculate the area of the contour
			defectf_area_pix = cv2.contourArea(c)
			defectf_area_ang = defectf_area_pix*conv
			# add area and distance data to defect frequency file
			f2.write("%d %f %f \n" % (count, defectf_area_ang, compare))
 
		# compare the protein shape to the shape of the defects directly below it
		prot_defect_area = np.count_nonzero(prot_defect)*conv
		prot_area = np.count_nonzero(prot_mask)*conv
		compare = cv2.matchShapes(prot_mask, prot_defect, 1, 0)
		f1.write("%d %f %f %f \n" % (count, prot_defect_area, prot_area, compare))
