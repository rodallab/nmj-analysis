#@ File (label="Choose a folder to analyze", style = "directory") path

MAXIP=false; //Whether to make a maxIP and analyze instead of 3D stack
CHANS=newArray("Brp", "Dyn", "Nwk"); // Channel labels, if don't want to use, uncomment the next lineCHANS=newArray("Brp", "Dyn", "Nwk"); // Channel labels, if don't want to use, uncomment the next line
//CHANS=newArray(0);
IMAGEEXT = "tif"; //extension for images
SEGMENT=0; //This is the channel to use to make a mask. Put 0 to combine multiple channels (see COMBINE parameter below).
COMBINE=newArray(1,2,3); // Which channels to add together to make mask IF segment==0.
NORMVAL="Max"; // Either Mean or Max; used to normalize channels before addding together
SIGMA=4; //This is the amount to blur your image to make the mask. Put 0 for no blur.

MED_RADIUS = 8;
MED_ITER = 0;
GAUS_RADIUS = 3;
GAUS_ITER= 0;
SIZEFILTER=0; 

ERODE=1; //Set to number of times to erode mask after thresholding. Put 0 for none. Put negative for dilation.
BACKSUB=50; //Radius for background subtraction. Put 0 for no background subtraction.
METHOD="Otsu"; // Threshold algorithm for segmentation
FILLHOLES = true;
SAVEMASK = true;
MAKEISO=false;
ipflag="";
if(MAXIP==true)
	ipflag="_2D";
SAVEPATH=path+File.separator+"WholeNMJ_RawData"+ipflag;
RECURSIVE=false;


run("CLIJ2 Macro Extensions", "cl_device=[Quadro P2000]");
Ext.CLIJ2_clear();
run("Options...", "iterations=1 count=1 black");
run("Colors...", "foreground=white background=black selection=yellow");
if(roiManager("count")>0)
	clearRoi();
if(nImages>0)
	close("*");
	
setBatchMode(true);
Table.create("CumulativeData");
if(File.exists(SAVEPATH) != true)
	File.makeDirectory(SAVEPATH);
folders=0;
Measure(path);

function Measure(path)
	{

	cumRow=Table.size("CumulativeData");
	count=getImgCount(path);
	fl=getFileList(path);
	Table.create(path);
	imgdone=0;
	
	
	for(i=0; i<fl.length; i++)
		{
		if(File.isDirectory(path+File.separator+fl[i]) && RECURSIVE==true)
			{
			Measure(path+File.separator+fl[i]);
			}

		if(endsWith(fl[i], IMAGEEXT)==true)
			{
			open(path+File.separator+fl[i]);
			getDimensions(width, height, channels, oslices, frames);
			getVoxelSize(vwidth, vheight, vdepth, vunit);
			if(MAKEISO==true)
				makeIso(fl[i]);
			if(MAXIP==true)
				{
				run("Z Project...", "projection=[Max Intensity]");
				close(fl[i]);
				rename(fl[i]);
				}
				
			Table.set("File", imgdone, fl[i], path);
			Table.set("File", cumRow+imgdone, fl[i], "CumulativeData");
			run("Select All");
			Stack.getDimensions(w,h,c,s,f);
			if((CHANS.length>0) && (CHANS.length != c))
				exit("Length of channel list doesn't match number of channels");
			segment(fl[i]);

			selectWindow(fl[i]);
			run("Select All");
			if(BACKSUB>0)
				run("Subtract Background...", "rolling="+BACKSUB);
			
			selectWindow("mask");
			run("Select All");
		
			if((SAVEMASK==true) && (MAKEISO==true))		
				{		
				// Save subset of mask slices as ROIs to hopefully match original image
				whratio=vdepth/vwidth;
				selectImage("mask");
				run("Select All");
				run("Duplicate...", "title=forROIs duplicate");
				setThreshold(1, 65535, "raw");
				run("Convert to Mask", "background=Dark black");
				
				for(s=0; s<oslices; s++)
					{
					oslice=round(s*whratio);
					Stack.setSlice(oslice);
					run("Create Selection");
					if(selectionType()>-1)
						{
						Roi.setPosition(1, s+1, 1);
						//run("Properties... ", "name=10 position=1,2,15 group=none stroke=yellow width=0 fill=none");
						roiManager("add");
						}
					}
		
				roiManager("deselect");
				roiManager("save", SAVEPATH+File.separator+stripX(fl[i])+"_masks.zip");
				roiManager("reset");
			}
				
			selectImage("mask");
			makeROISnicenice();
			
			if((SAVEMASK==true)&&(MAKEISO==false))
				saveROIs(SAVEPATH, fl[i], "mask"); 
			//exit();
			if(roiManager("count")>0)
				{	
				for(ch=1; ch<c+1; ch++)
					{
					selectImage(fl[i]);
					if(c>1)
						Stack.setChannel(ch);
					
					c1pix=getPixValues(ch);
					totals=measurePre(fl[i], ch);  // 0 is for Pre measurements
					Table.set("PreArea", imgdone, totals[1], path);
					Table.set("PreArea", cumRow+imgdone, totals[1], "CumulativeData");
					
					if(CHANS.length == 0)
						{
						Table.set("Ch-"+ch+"_Mean", imgdone, totals[0], path);
						Table.set("Ch-"+ch+"_CoV", imgdone, totals[2], path);		
						Table.set("Ch-"+ch+"_Mean", cumRow+imgdone, totals[0], "CumulativeData");
						Table.set("Ch-"+ch+"_CoV", cumRow+imgdone, totals[2], "CumulativeData");		
						
						}
						else {
						Table.set(CHANS[ch-1]+"_Mean", imgdone, totals[0], path);
						Table.set(CHANS[ch-1]+"_CoV", imgdone, totals[2], path);
						Table.set(CHANS[ch-1]+"_Mean", cumRow+imgdone, totals[0], "CumulativeData");
						Table.set(CHANS[ch-1]+"_CoV", cumRow+imgdone, totals[2], "CumulativeData");
						}
	
				if(c>1)
					{
					for(c2=ch+1; c2<c+1; c2++)
						{
						Stack.setChannel(c2);
						c2pix=getPixValues(c2);
						pearson=PCC(c1pix, c2pix);
						if(CHANS.length == 0)
							{
							Table.set("PCC_Ch"+ch+"-Ch"+c2, imgdone, pearson, path);
							Table.set("PCC_Ch"+ch+"-Ch"+c2, cumRow+imgdone, pearson, "CumulativeData");		
							
							}
						else {
							Table.set("PCC_"+CHANS[ch-1]+"-"+CHANS[c2-1], imgdone, pearson, path);
							Table.set("PCC_"+CHANS[ch-1]+"-"+CHANS[c2-1], cumRow+imgdone, pearson, "CumulativeData");
							}
						}
					}
				}
				roiManager("reset");
				}
			else
				{
				for(ch=1; ch<c+1; ch++)
					Table.set("PreArea", imgdone, NaN, path);
					Table.set("PreArea", cumRow+imgdone, NaN, "CumulativeData");
					if(CHANS == "")
						{
						Table.set("Ch-"+ch+"_Mean", imgdone, NaN, path);
						Table.set("Ch-"+ch+"_CoV", imgdone, NaN, path);
						Table.set("Ch-"+ch+"_Mean", cumRow+imgdone, NaN, "CumulativeData");
						Table.set("Ch-"+ch+"_CoV", cumRow+imgdone, NaN, "CumulativeData");
						Table.set("PCC_Ch"+ch+"-Ch"+c2, imgdone, NaN, path);
						Table.set("PCC_Ch"+ch+"-Ch"+c2, cumRow+imgdone, NaN, "CumulativeData");		
						}
						else 
						{
						Table.set(CHANS[ch-1]+"_Mean", imgdone, NaN, path);
						Table.set(CHANS[ch-1]+"_CoV", imgdone, NaN, path);
						Table.set(CHANS[ch-1]+"_Mean", cumRow+imgdone, NaN, "CumulativeData");
						Table.set(CHANS[ch-1]+"_CoV", cumRow+imgdone, NaN, "CumulativeData");
						Table.set("PCC_"+CHANS[ch-1]+"-"+CHANS[c2-1], imgdone, NaN, path);
						Table.set("PCC_"+CHANS[ch-1]+"-"+CHANS[c2-1], cumRow+imgdone, NaN, "CumulativeData");
						}
				}
			close("*");	
			imgdone+=1;
			}
		//print(imgdone);
		Table.save(SAVEPATH+File.separator+"Presynaptic-Intensity_CoV.csv");	
		}
		
	}

if(folders>1)	
	Table.save(SAVEPATH+File.separator+"CumulativeData.csv");

function measurePre(img, ch)
	{
	nROI=roiManager("count");
	TotalArea=0;
	TotalID=0;
	WeightedCoV=0;
	selectWindow(img);
	
	//run("Median...", "radius=1");
	if(nROI==1)
		{
		roiManager("Select", 0);
		Stack.setChannel(ch);
		getStatistics(area, mean, min, max, std); //THIS IS THE PART WITH THE MEASUREMENT
		TotalArea=area;
		TotalID=area*mean;
		WeightedCoV=std/mean*area;
		
		}
		else
			{
			for(r=0; r<nROI; r++)
				{
				roiManager("Select", r);
				Stack.setChannel(ch);
				getStatistics(area, mean, min, max, std); //THIS IS THE PART WITH THE MEASUREMENT
				TotalArea+=area;
				TotalID+=area*mean;
				if(mean>0)
					WeightedCoV=WeightedCoV+std/mean*area;
				//print("r:"+r+"  WCoV:"+WeightedCoV+"  std:"+std);
				}
			}
			
	TotalMean=TotalID/TotalArea;
	TotalCoV=WeightedCoV/TotalArea;
	//print("totalCoV:"+TotalCoV);
	return newArray(TotalMean, TotalArea, TotalCoV);
	}

function getImgCount(path)
	{
	count=0;
	fl=getFileList(path);
	for(i=0; i<fl.length; i++)
		{
		if(endsWith(fl[i], IMAGEEXT)==true)
			count=count+1;
		}
	return count;
	}
	
function segment(imp)
	{
	selectImage(imp);
	Stack.getDimensions(width, height, channels, slices, frames);
	if(SEGMENT==0)
		{
		run("Duplicate...", "title=tocombine duplicate channels=1-"+channels); //mask for whole NMJ
		if(BACKSUB>0)
			run("Subtract Background...", "rolling="+BACKSUB+ " stack");
		if(slices>1)
			{
			run("Z Project...", "projection=[Max Intensity]");
			rename("forNorm");
			}
		else
			run("Duplicate...", "title=forNorm duplicate channels=1-"+channels); 
		
	
		selectImage("tocombine");
		run("32-bit");
		run("Split Channels");
		newImage("mask", "32-bit black", width, height, slices);
		
		// normalize each channel to mean to combine evenly
		for(c=1; c<channels+1; c++)
			{
			selectImage("forNorm");
			Stack.setChannel(c);
			run("Select All");
			overallmean=getValue(NORMVAL);
			
			selectImage("C"+c+"-tocombine");
			run("Divide...", "value="+overallmean+" stack");
			}
	
		//add channels to mask
		for(c=0; c<COMBINE.length; c++)
			{
			imageCalculator("Add stack", "mask", "C"+COMBINE[c]+"-tocombine");
			rename("mask");
			}
		
		selectImage("mask");
		
		getMinAndMax(min, max);
		if(slices>1)
			{
			
			for(s=1; s<=slices; s++)
				{
				Stack.setSlice(s);
				run("Select All");
				max0 = getValue("Max");
				min0 = getValue("Min");
				
				if(min0<min)
					min=min0;
				if(max0>max)
					max=max0;
				}
			}
			
		setMinAndMax(min, max);
		run("16-bit");
		
		}
		else run("Duplicate...", "title=mask duplicate channels="+SEGMENT); //mask for whole cell
	
	if(BACKSUB>0)
		run("Subtract Background...", "rolling="+BACKSUB+ " stack");
		
	if(SIGMA>0)
		run("Gaussian Blur...", "sigma="+SIGMA+" stack");
	if(MED_ITER>0)
		recMedFilter(MED_ITER, MED_RADIUS);	
		
	if(GAUS_ITER>0)
		recGausFilter(GAUS_ITER, GAUS_RADIUS);
	

	setAutoThreshold(METHOD+" dark stack");	
	run("Convert to Mask", "method="+METHOD+" background=Dark black stack");	

	if(FILLHOLES==true)
		run("Fill Holes", "stack");

	rename("mask");


	if(ERODE>0)
		{
		for(e=0; e<ERODE; e++)
			run("Erode", "stack");
		}
	
	if(ERODE<0)
		{
		for(e=0; e<abs(ERODE); e++)
			run("Dilate", "stack");
		}
	
	if(SIZEFILTER>0)
		{
		run("Select All");
		run("Analyze Particles...", "size="+SIZEFILTER+"-Infinity add");
		nRoi=roiManager("count");
		rois=Array.getSequence(nRoi);
		roiManager("select", rois);
		roiManager("Combine");
		run("Clear Outside");
		roiManager("reset");
		}
	
	}

function saveROIs(path, imp, suffix)
	{
	rois=newArray(roiManager("count"));
	for (n=0; n<rois.length; n++)
		{
		rois[n]=n;
		}
	roiManager("select", rois);
	roiManager("Save", path+File.separator+stripX(imp)+suffix+".zip");
	}

function clearRoi()
	{
	while(roiManager("count")>0)
		{
		roiManager("Select", 0);
		roiManager("Delete");
		}
	}

function stripX(string)
	{
	// This is because Macro language doesn't have a general use name without extension
	return substring(string, 0, lastIndexOf(string, "."));
	}

function makeROISnicenice()
	{
	getDimensions(width, height, channels, slices, frames);
	if(frames>slices)
		slices=frames;
	Stack.getPosition(c, slicepos, f);
	for(s=1; s<slices+1; s++)
		{
		Stack.setPosition(c, s, f);
		run("Create Selection");
		if(selectionType()>-1)
			{
			Roi.setPosition(c, s, f);
			roiManager("add");
			}
		}
	}

function getPixValues(chan)
	{
	// Get pixel intensity values for all points within an ROI and return as an array
	allvals=newArray(0);
	for(r=0; r<roiManager("count"); r++)
		{
		roiManager("select", r);	
		Roi.getContainedPoints(xpoints, ypoints);
		vals=newArray(xpoints.length);
		Stack.setChannel(chan);
		for(n=0; n<xpoints.length; n++)
			vals[n]=getPixel(xpoints[n], ypoints[n]);
		allvals = Array.concat(allvals, vals);
		}
	return allvals;
	}
	
function PCC(X,Y)
	{
	// Calculate the Pearson Correlation coefficient between two arrays
	if(X.length != Y.length)
		exit("array lengths must match in PCC");
	nan=false;
	// First check whether either array contains nan and remove these elements from both
	// do check first bc this is slow. The built in imagej 'filter' is useless for this :(
	for(i=0; i<X.length; i++)
		if(X[i]!=X[i])
			nan=true;
	if(nan==true)
		{
		y_filt = filterNaN(X, Y);
		x_filt = filterNaN(X, X);
		}
		else
		{
		x_filt = X;
		y_filt = Y;
		}
		
	nan=false;
	for(i=0; i<y_filt.length; i++)
		if(y_filt[i]!=y_filt[i])
			nan=true;
	if(nan==true)
		{
		x_filt = filterNaN(y_filt, x_filt);
		y_filt = filterNaN(y_filt, y_filt);	
		}
	Array.getStatistics(x_filt, Xmin, Xmax, Xmean, XstdDev);
	Array.getStatistics(y_filt, Ymin, Ymax, Ymean, YstdDev);

	xy=0; xsq=0; ysq=0;
	for(n=0; n<x_filt.length; n++)
		{
		x=x_filt[n]-Xmean;
		y=y_filt[n]-Ymean;
		xy+=(x*y);
		xsq+=(x*x);
		ysq+=(y*y);
		}

	pearson=(xy/sqrt(xsq*ysq));

	return pearson;
	}
	
function makeIso(img)
	{
	selectImage(img);
	Stack.getDimensions(width, height, channels, slices, frames);
	getVoxelSize(vwidth, vheight, vdepth, vunit);
	bitdepth=""+bitDepth()+"-bit";
	if(channels>1)
		{
		run("Split Channels");
		mergeString="";
		for(c=1; c<channels+1; c++)
			{
			currentimg="C"+c+"-"+fl[i];
			Ext.CLIJ2_push(currentimg);
			close(currentimg);
			iso="iso-C"+c+"-"+fl[i];
			Ext.CLIJ2_makeIsotropic(currentimg, iso, vwidth, vheight, vdepth, vwidth);		
			Ext.CLIJ2_pull(iso);
			mergeString+="c"+c+"=iso-C"+c+"-"+stripX(fl[i])+".tif "; 
			}
		mergeString+="create";
		run("Merge Channels...", mergeString);
		}
	else 
		{
		Ext.CLIJ2_push(img);
		close(img);
		iso = "iso";
		Ext.CLIJ2_makeIsotropic(img, iso, vwidth, vheight, vdepth, vwidth);		
		Ext.CLIJ2_pull(iso);
		}
	rename(fl[i]);
	setVoxelSize(vwidth, vheight, vwidth, vunit);
	}
	
function recGausFilter(iter, radius)
	{
	Ext.CLIJ2_clear();
	currentimg=getTitle();
	Ext.CLIJ2_push(currentimg);
	close(currentimg);
	
	for(ri=0; ri<iter; ri++)
		{
		nextimg=currentimg+toString(ri);
		Ext.CLIJ2_gaussianBlur3D(currentimg, nextimg, radius, radius, 0);
		currentimg=nextimg;
		}
	Ext.CLIJ2_pull(nextimg);
	}

function recMedFilter(iter, radius)
	{
	Ext.CLIJ2_clear();
	currentimg=getTitle();
	Ext.CLIJ2_push(currentimg);
	close(currentimg);
	
	for(ri=0; ri<iter; ri++)
		{
		nextimg=currentimg+toString(ri);
		Ext.CLIJ2_medianSliceBySliceSphere(currentimg, nextimg, radius, 1);
		currentimg=nextimg;
		}
	Ext.CLIJ2_pull(nextimg);
	}
	
/* ------ VERSION HISTORY ---------
 *  v5 
 *  Fixed Analyze particles and ROI bug 
 *  
 *  v8
 *  fixed bug in PCC
 *  
 *  v9
 *  added option to do recursive blurs and choose max normalization before adding channels to mask
 */