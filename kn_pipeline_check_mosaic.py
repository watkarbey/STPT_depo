#!/usr/bin/env python
import ntpath
from scipy import io as sio
from os import listdir
import os
import sys
from os.path import isfile, join, isdir
import json
import re
import math
import argparse #https://docs.python.org/2/library/argparse.html
#from kn_pipeline_common import *
#from kn_common import pipetools,pipedef

tissue_cyte_exclude_sub_dirs={'additional_data','Fused'};
tissue_cyte_exclude_subsub_dirs={'trigger'};

parser = argparse.ArgumentParser(description='tissuce cyte check script')
parser.add_argument('--folder')
parser.add_argument('--debug',type=bool,default=False)
parser.add_argument('--mask',default="")
parser.add_argument('--channels',default=3,type=int)
parser.add_argument('--optical_sections_missing',action='store_true')

args = parser.parse_args()

debug=args.debug;
mask=args.mask;

channels = args.channels

firstfail=False;
last_img_id=-1;

def mosaic_read(filename):
		try:
			num_coords=0;
			if (os.path.isfile(filename)): 
					data = {};
					with open(filename) as data_file:    
							for line in data_file:
									comp=line.replace('\r','').replace('\n','').split(':',1);
									data[comp[0]]=comp[1];		
									if "XPos" in comp[0]:
										num_coords+=1
					isfile=True;
			else:
					data={}
					isfile=False;
					
			if "layers" in data:
				if not args.optical_sections_missing:
					num_coords=num_coords*int(data["layers"])

			else:
				sys.stderr.write("#W: no layers attribute in mosaic file "+filename+"\n")
		except: 
			sys.stderr.write("#E: error parsing "+filename+"\n")
			data={}
			isfile=False;
		return data, isfile, num_coords


def get_mosaic_files(folder,d=1,FID="",mosaic_expected=False):	    
	global firstfail
	global last_img_id
	global debug
	global mask
	
	#mosaicfiles=[f for f in listdir(folder) if isfile(join(folder, f)) and ( (f[len(f)-3:len(f)]=='txt') and ("Mosaic" in f))]

	allfiles_and_folders=listdir(folder);
	allfiles=[f for f in allfiles_and_folders if isfile(join(folder, f))];
	mosaicfiles=[f for f in allfiles if ( (f[len(f)-3:len(f)]=='txt') and ("Mosaic" in f))]
	#alltiff=[f for f in allfiles if ( (f[len(f)-3:len(f)]=='tif') )]
	alltiff=[f for f in allfiles if ( (f[len(f)-3:len(f)]=='tif')or (f[len(f)-3:len(f)]=='png')  )]	
	
	
		
	mfiles={};
	ok=True;

	if (not mosaic_expected) and (d>1):
		return mfiles,  not (len(mosaicfiles)>0), alltiff

	# ********************
	# we only accept 1! mosaic file per folder and files layer
	# ********************
	if d>1:
		if len(mosaicfiles)!=1:
				#sys.stderr.write("#E: found more or less Mosaic files than expected ("+format(len(mosaicfiles))+" instead of 1) in folder\n#E: "+folder+"\n")
				sys.stderr.write("#E: found more or less Mosaic files than expected ("+format(len(mosaicfiles))+" instead of 1) in folder\n")
				sys.stderr.write("#E: "+folder+"\n")
				ok=False;
							
				for mfile  in  mosaicfiles:	
					sys.stderr.write("#I: file:"+folder+'/'+mfile+"\n")

	# ********************
	# read the header mosaic in folder layer and calculate expected number of images
	# ********************
	total_num_expected_files=0;
	if d==2:
		for mfile  in  mosaicfiles:	
			mf,succ,numc=mosaic_read(folder+'/'+mfile);
			if not succ:
				sys.stderr.write("#E: could not read the Mosaic file "+folder+'/'+mfile+" :\n")
				sys.stderr.write("#E: cannot compute the number of expected files.\n")
				sys.stderr.write("#E: -> please fix this error first and run again before fixing other errors.\n")
				ok=False;
			else:
				mfiles[folder+'/'+mfile]=mf;
				total_num_expected_files=int(mf['mrows'])*int(mf['mcolumns'])
				if not args.optical_sections_missing:
					total_num_expected_files*=int(mf['layers'])


	# ********************
	# check mosaic file name in file layer
	# check for world coordinate vectors
	# ********************
	if (d==3):
		for mfile  in  mosaicfiles:	
			if not FID in mfile:
				sys.stderr.write("#E: name of mosaic file "+folder+'/'+mfile+" is incorrect:\n")
				ok=False;
			mf,succ,numc=mosaic_read(folder+'/'+mfile);
			if not succ:
				sys.stderr.write("#E: could not read the Mosaic file "+folder+'/'+mfile+" :\n")
				sys.stderr.write("#E: -> please fix this error first and run again before fixing other errors.\n")
				ok=False;
			else:	
			  alltiff_C01=[f for f in alltiff if "_01." in f]
			  if len(alltiff_C01) != (numc):
			    sys.stderr.write("#E: number of files differs from number of  coordinates in "+folder+'/'+mfile+":\n")
			    sys.stderr.write("#E: number files: "+format(len(alltiff_C01))+"\n")
			    sys.stderr.write("#E: number coordinates: "+ format((numc))+"\n")
			    ok=False;
			
	# ********************
	# check for file id numbers
	# ********************

			  if len(alltiff)>0:
				  total_num_expected_files=int(mf['mrows'])*int(mf['mcolumns'])
				  if not args.optical_sections_missing:
					  total_num_expected_files*=int(mf['layers'])
				  for channel in range(1,channels+1):
				      #tiff_f_c=[int(f.split("-")[2].split("_")[0]) for f in alltiff if "_0"+format(channel)+".tif" in f];
				      tiff_f_c=[int(f.split("-")[2].split("_")[0]) for f in alltiff if (("_0"+format(channel)+".tif" in f) or ("_0"+format(channel)+".png" in f))];

				      tiff_f_c=sorted(tiff_f_c)
				      
				      startnum=int(mf['startnum'])
				      endnum=startnum+total_num_expected_files-1
				      
				      if len(tiff_f_c)>0:
                                        lower_indx=0;
                                        higher_indx=0;
                                        for ncount in range(1,len(tiff_f_c)):
                                                f=tiff_f_c[ncount]
                                                if f<startnum:
                                                    lower_indx+=1;
                                                if f>endnum:    
                                                    higher_indx+=1;
                                                if ncount<len(tiff_f_c)-1:
                                                    f1=tiff_f_c[ncount+1]
                                                else:                                                    
                                                    f1=endnum
                                                step=f1-f;
                                                if step>1:
                                                    sys.stderr.write("#E: In folder "+folder+",\n")
                                                    sys.stderr.write("#E: channel  "+format(channel)+":\n")
                                                    sys.stderr.write("#E: missing "+format(step-1)+" imgs, range  ["+format(f+1)+" : "+format(f1-1)+"]\n")
                                                    ok=False;
                                                        
				      if  lower_indx>0:
                                                    sys.stderr.write("#E: In folder "+folder+",\n")
                                                    sys.stderr.write("#E: channel  "+format(channel)+":\n")
                                                    sys.stderr.write("#E: unexpected img ids, "+format(lower_indx)+" images with ids <  "+format(startnum)+"\n")
                                                    ok=False;
                                          
				      if  higher_indx>0:
                                                    sys.stderr.write("#E: In folder "+folder+",\n")
                                                    sys.stderr.write("#E: channel  "+format(channel)+":\n")
                                                    sys.stderr.write("#E: unexpected img ids, "+format(higher_indx)+" images with ids >  "+format(endnum)+"\n")
                                                    ok=False;
                                                                                           
				      	      
			  
	
	# ********************
	# get FID (mosaic file id)
	# ********************
	MFNAME="";
	if (len(mosaicfiles)==1) and (d==2):
			for mfile in mosaicfiles:
				#FID=mfile[7:len(f)-4]
				FID=mfile[7:len(mfile)-4]
				MFNAME=mfile;
				
	# ********************
	# get folders (sessions) (d==1) in tissucyte folder
	# ********************
	if len(FID)==0:
			onlydirs=[f for f in listdir(folder) if isdir(join(folder, f)) and (f not in tissue_cyte_exclude_sub_dirs) and ((len(mask)==0) or (mask in f))]	
	# ********************
	# check folder (d==2) and file layer (d==3) for folder conistency
	# ********************
	else:
			#print FID
			if debug and (d<3):
				print("DEBUG: FID: "+FID)

			onlydirs=[f for f in listdir(folder) if isdir(join(folder, f)) and 
					(	
					 	(FID in f) and 
					 	(f not in tissue_cyte_exclude_subsub_dirs) or (d==3))
					 ]	
			crosscheck=[f for f in listdir(folder) if isdir(join(folder, f)) and 
						(
							  (FID not in f) and 
							((f not in tissue_cyte_exclude_subsub_dirs) or (d==3)) ) ]


			if (len(crosscheck)>0):
						sys.stderr.write("#E: A missmatch between folder name and Mosaic file?!\n")
						sys.stderr.write("#E: ID given by Mosaic file name: \""+FID+"\". \n")
						sys.stderr.write("#E: problematic folder names in "+folder+":\n")
						for f in crosscheck:
								sys.stderr.write("#E: "+f+"\n")
						ok=False;

	if (d==1) and (len(onlydirs)==0):
			sys.stderr.write("#E: no subdirectories at all!\n")
			ok=False;
			
	# ********************
	# check if number of folders in (d==2) layer is identical to the number osections
	# ********************			
	if (d==2) and (len(mosaicfiles)==1):
			num_sections=int(mf['sections'])
			if debug:
				print("DEBUG: sections: "+format(num_sections))
			if len(onlydirs)!=num_sections:
			  sys.stderr.write("#E: the number of sections defined in \n")
			  sys.stderr.write("#E: "+folder+'/'+mfile+" \n")
			  sys.stderr.write("#E: differs from the number of section folders ("+format(num_sections)+" vs "+format(len(onlydirs))+") \n")
			  ok=False;
	
	alltiff_files=[];

	subfolder_id=0;
	for subfolder  in  onlydirs:	
		# ********************
		#mosic folder layer
		# ********************
		if d==1:
			sub_mfiles, subok ,tiff_f = get_mosaic_files(folder+'/'+subfolder,d+1,FID,True);
			ok = subok and ok;
			mfiles.update(sub_mfiles);
		# ********************
		#mosic files layer
		# ********************
		if d==2:
			if debug:
				print("DEBUG: FID: "+FID)
			subfolder=folder+'/'+subfolder;
			if debug:
				print("DEBUG: checking folder:"+subfolder)
			sub_mfiles, subok , tiff_f = get_mosaic_files(subfolder,d+1,FID,True);	
			alltiff_files=alltiff_files+tiff_f;
			ok = subok and ok;
			subfolder_id=subfolder_id+1;
			

			# ********************
			#checking if tiff file number is consisntent with the mosaic info
			# ********************
			if len(tiff_f)>0:
				for channel in range(1,channels+1):
					channel_count=0;
					for f in tiff_f:
						sig=[int(f) for f in re.split("(\d+)", f) if (re.match("(\d+)",f))];
						if sig[3]==channel:
							channel_count=channel_count+1;

					if total_num_expected_files!=channel_count:
						sys.stderr.write("#E: "+subfolder+",\n#W: CHANNEL "+str(channel)+":\n")
						sys.stderr.write("#E: # of expected tiff files ("+format(total_num_expected_files)+") differs from their actual number ("+format(channel_count)+")\n")
						sys.stderr.write("#E: -> the final results might be corrupted\n")
						ok=False;
						

	

	if (not ok) and (not firstfail):
		sys.stderr.write("#E: \n")
		sys.stderr.write("#E: ERROR ! \n")
		sys.stderr.write("#E: tissue data would not pass the pipeline\n")
		sys.stderr.write("#E: continue checking for further problems \n")
		sys.stderr.write("#E: \n")
		firstfail=True;
	#else:
	#	sys.stderr.write("tissue data would pass the pipeline (so far)\n")
	return mfiles, ok, alltiff


folder=args.folder#"/disk/waikiki/KAKUSHIN-NOU-DATA-SUB/raw/CM696/tissuecyte/" 
print("#I: checking folder "+folder)
mfiles, ok , tiff_f =get_mosaic_files(folder)

#print "ok ?: "+format(ok)

for key, value in mfiles.items():
	if not value['Sample ID'] in ntpath.basename(key):
		sys.stderr.write("#W: name of mosaic file "+ntpath.basename(key)+' differs from Sample ID'+value['Sample ID']+"\n")
		sys.stderr.write("#W: Full path: "+key+"\n")
	# ********************
	#check resolution
	# ********************
	if not float(value['zres']) <= 50:
		ok = False;
		sys.stderr.write("#E: zres of  mosaic file "+ntpath.basename(key)+' is larger than 50 :'+format(value['zres'])+"\n")
	# ********************
	#check sectionres
	# ********************
	if not float(value['sectionres']) == 50:
		ok = False;
		sys.stderr.write("#E: sectionres of  mosaic file "+ntpath.basename(key)+' is not 50 :'+format(value['sectionres'])+"\n")


	# ********************
	#????
	# ********************
if  (len(mfiles)==0):
		sys.stderr.write("#E: no mosaic files at all!\n")
		ok=False;


if not ok:
	sys.stderr.write("#E: Tissue cyte test failed\n");
	sys.exit(1);
else:
	sys.stderr.write("#I: Tissue cyte test passed\n");
	sys.exit(0);


