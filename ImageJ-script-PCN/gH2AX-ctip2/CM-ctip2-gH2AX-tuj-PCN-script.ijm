//To analyze all images 

dir="G:/PCN/CM-PCN-wt-D0-052821/gH2AX/gH2AX-ctip2/Data/"; 

sub=getFileList(dir);
//Array.print(sub);

//Loop through each individual folder
for(i=0;i<lengthOf(sub);i++) {
	//Choose an individual genotype/timepoint combination
	dir1=dir+sub[i];
	//Decide if continuing or not
	completed=getFileList(dir1+"data/");
	//print(lengthOf(completed)); }

	//if(lengthOf(completed)>=4) {} else {

		images=getFileList(dir1+"proj_max/");
		//Array.print(images);
		images_length=lengthOf(images);
		//Loop through each proj_max image
		for(j=0;j<images_length;j++) { 
			//open image and get title and rename - Hoescht image
			//print(images[j]);
			//print(matches(images[j],"[A-Z][0-9]{2}_w2_s[0-9]{1}_.tif")); } }
			
			//Run twice, switching which line is used; top for sites 1-9, bottom for sites>=10
			//if(matches(images[j],"[A-Z][0-9]{2}_w3_s[0-9]{1}_.tif")==1) {
			if(matches(images[j],"[A-Z][0-9]{2}_w2_s[0-9]{2}.tif")==1) {
			
			//print(images[j]); } } }
			open(dir1+"proj_max/"+images[j]);
			name=getTitle();
			run("Clear Results");
			run("Set Measurements...", "area limit redirect=None decimal=2");
			run("Duplicate...", "title=ctip2");
			run("Properties...", "channels=1 slices=1 frames=1 unit=um pixel_width=0.17 pixel_height=0.17 voxel_depth=0.17");

			//open Hoescht 
			image2=replace(images[j],"_w2_","_w1_");
			open(dir1+"proj_max/"+image2);
			run("Duplicate...", "title=Hoescht");
			run("Properties...", "channels=1 slices=1 frames=1 unit=um pixel_width=0.34 pixel_height=0.34 voxel_depth=0.34");

			//open gH2AX
			image3=replace(images[j],"_w2_","_w3_");
			open(dir1+"proj_max/"+image3);
			run("Duplicate...", "title=gH2AX");
			run("Properties...", "channels=1 slices=1 frames=1 unit=um pixel_width=0.34 pixel_height=0.34 voxel_depth=0.34");

//			//open proj
//			image4=replace(images[j],"_w2_","_w4_");
//			open(dir1+"proj_max/"+image4);
//			run("Duplicate...", "title=proj");
//			run("Properties...", "channels=1 slices=1 frames=1 unit=um pixel_width=0.34 pixel_height=0.34 voxel_depth=0.34");


///Identify dead cells from Hoescht image; make Hoescht and dead Hoescht mask 
			
			selectWindow("Hoescht");
			run("Duplicate...","title=Hoescht1");
			run("Median...", "radius=3");
			selectWindow("Hoescht");
			run("Duplicate...","title=Hoescht2");
			run("Gaussian Blur...", "sigma=80");
			imageCalculator("Subtract create","Hoescht1","Hoescht2");

			run("Duplicate...", "title=Hoescht-live");
			run("Duplicate...", "title=Hoescht-dead");
	     	setThreshold(15000, 65535); 
			run("Make Binary");
			run("Analyze Particles...", "size=0-infinity circularity=0-1.00 show=Masks");
			run("Invert");
			run("Duplicate...", "title=Hoescht-dead-1");
			run("Options...", "iterations=8 count=1 black do=Dilate");
			run("Invert");
			
			selectWindow("Hoescht-live");
			setThreshold(1900, 65535); 
			run("Make Binary");
			run("Analyze Particles...", "size=0-infinity circularity=0-1.00 show=Masks");
			close("Hoescht-live"); close("Hoescht-dead"); close("Hoescht1"); close("Hoescht2"); 
			close("Result of Hoescht1");

			selectWindow("Mask of Hoescht-live");
			run("Invert");
			run("Fill Holes");
			run("Options...", "iterations=2 count=1 black do=Dilate");
			run("Options...", "iterations=1 count=1 black do=Close");
			run("Watershed");
			run("Options...", "iterations=2 count=1 black do=Erode");
			run("Options...", "iterations=1 count=1 black do=Close");
			run("Duplicate...", "title=Hoescht-live-1");
			run("Invert");

			imageCalculator("Subtract create","Hoescht-live-1","Hoescht-dead-1");
			run("Invert");
			run("Options...", "iterations=1 count=1 black do=Dilate");
			
			run("Analyze Particles...", "size=49-500 circularity=0.3-1.00 show=Masks");
			run("Duplicate...", "title=TujNucl4");
			run("Invert");
			run("Fill Holes");
			close("Result of Hoescht-live-1");close("Mask of Result of Hoescht-live-1");
			close("Hoescht-live-1"); close("Mask of Hoescht-live"); close("Mask of Hoescht-dead");
			close("Hoescht-dead-1");
		
			run("Analyze Particles...", "size=0-infinity circularity=0-1 summarize");
			saveAs("Results", dir1+"data/Hoescht/"+name+"_hoescht_count.csv"); ///
			close(name+"_hoescht_count.csv"); ///
			run("Clear Results");close("Summary");


///Identify dead cells and CTIP2+ cells from CTIP2 image		
			selectWindow("ctip2");
			run("Duplicate...","title=ctip2-1");
			run("Median...", "radius=3");
			selectWindow("ctip2");
			run("Duplicate...","title=ctip2-2");
			run("Gaussian Blur...", "sigma=80");
			imageCalculator("Subtract create","ctip2-1","ctip2-2");
			close("ctip2-1"); close("ctip2-2"); 

			run("Duplicate...", "title=ctip2-live");
			run("Duplicate...", "title=ctip2-dead");
			setThreshold(4000, 65535); 
			run("Make Binary");
			run("Analyze Particles...", "size=0-infinity circularity=0-1.00 show=Masks");
			run("Invert");
			run("Duplicate...", "title=ctip2-dead-1");
			run("Options...", "iterations=8 count=1 black do=Dilate");
			run("Invert");
			
			selectWindow("ctip2-live");
			setThreshold(500, 65535); 
			run("Make Binary");
			run("Analyze Particles...", "size=0-infinity circularity=0-1.00 show=Masks");
			close("ctip2-live"); close("ctip2-dead"); 
			close("Result of ctip2-1");

			selectWindow("Mask of ctip2-live");
			run("Invert");
			run("Fill Holes");
			run("Options...", "iterations=2 count=1 black do=Dilate");
			run("Options...", "iterations=1 count=1 black do=Close");
			run("Watershed");
			run("Options...", "iterations=2 count=1 black do=Erode");
			run("Options...", "iterations=1 count=1 black do=Close");
			run("Duplicate...", "title=ctip2-live-1");
			run("Invert");

			imageCalculator("Subtract create","ctip2-live-1","ctip2-dead-1");
			run("Invert");
			run("Options...", "iterations=1 count=1 black do=Dilate");
			
			run("Analyze Particles...", "size=25-infinity circularity=0.3-1.00 show=Masks");
			run("Duplicate...", "title=ctip2pos-ori");
			run("Invert");
			run("Fill Holes"); close("ctip2-dead-1");
			close("Result of ctip2-live-1");close("Mask of Result of ctip2-live-1");
			close("ctip2-live-1"); close("Mask of ctip2-live"); close("Mask of ctip2-dead");

/// Make final CTIP2+/- masks for nuclei
			
			// Make final CTIP- nuclei mask
			selectWindow("ctip2pos-ori");
			run("Duplicate...","title=ctip2pos-ori1");
			run("Options...", "iterations=10 count=1 black do=Dilate");
			run("Invert LUT"); run("Invert");	
			
			selectWindow("TujNucl4");
			run("Invert LUT"); run("Invert");
			imageCalculator("Subtract create", "TujNucl4", "ctip2pos-ori1");
			run("Analyze Particles...", "size=20-infinity circularity=0.6-1.00 show=Masks");		
			run("Duplicate...","title=ctip2neg-nucl-fin"); run("Invert LUT");
			run("Analyze Particles...", "size=0-infinity circularity=0-1 summarize");
			saveAs("Results", dir1+"data/ctip2_tuj/"+name+"_ctip2_NEG_nucl_count.csv"); ///
			close(name+"_ctip2_NEG_nucl_count.csv"); ///
			run("Clear Results");close("Summary");
			close("Result of TujNucl4"); close("Mask of Result of TujNucl4");
			
			// Make final CTIP+ nuclei mask 
			selectWindow("ctip2neg-nucl-fin");
			run("Duplicate...","title=ctip2neg-nucl-1");
			run("Options...", "iterations=10 count=1 black do=Dilate");
			imageCalculator("Subtract create", "TujNucl4", "ctip2neg-nucl-fin");
			run("Analyze Particles...", "size=20-infinity circularity=0.6-1.00 show=Masks");		
			run("Duplicate...","title=ctip2pos-nucl-fin"); run("Invert LUT");
			run("Analyze Particles...", "size=0-infinity circularity=0-1 summarize");
			saveAs("Results", dir1+"data/ctip2_tuj/"+name+"_ctip2_POS_nucl_count.csv"); ///
			close(name+"_ctip2_POS_nucl_count.csv"); ///
			run("Clear Results");close("Summary");
			close("Result of TujNucl4"); close("Mask of Result of TujNucl4");

			close("ctip2neg-nucl-1");close("ctip2pos-ori1");

///////// gH2AX Analysis: 

	///Quantify gH2AX intensity in cytobody and nucleus for CTIP2+ cells 

			/// Remove bright signal and measure intensity of gH2AX in CTIP2+ nuclei & cyto ("ctip2pos-body-fin" mask)	
			selectWindow("gH2AX");
			run("Duplicate...","title=gH2AX-b"); run("8-bit");
			setThreshold(180,255);
			run("Make Binary");
			run("Options...", "iterations=10 count=1 black do=Dilate");
			run("Invert");
			imageCalculator("AND create","gH2AX-b","ctip2pos-nucl-fin");
			run("Analyze Particles...", "size=21-infinity circularity=0.6-1.00 show=Masks");
			run("Invert LUT");	
			run("Duplicate...", "title=ctip2pos-nucl-1"); close("Result of gH2AX-b");
		  
			run("Set Measurements...", "area mean standard min feret's integrated limit redirect=None decimal=2");
			run("Analyze Particles...", "size=0-infinity pixel circularity=0-1 show=Outlines add in_situ");
			
			selectWindow("gH2AX");
			run("Duplicate...","title=gH2AX-a");
			run("Gaussian Blur...", "sigma=100");
			imageCalculator("Subtract create","gH2AX","gH2AX-a");
			run("Duplicate...", "title=gH2AX-1");
	        setThreshold(900,65500); close("gH2AX-a"); 

	        if(roiManager("Count") > 0){
			run("From ROI Manager");
			}
			roiManager("Measure");
			saveAs("Results", dir1+"data/h2ax_tuj/"+name+"_h2ax_ctip2_POSnucl_intensity_count.csv"); ///
			close(name+"_h2ax_ctip2_POSnucl_intensity_count.csv"); ///
			run("Clear Results");close("Summary");close("ROI Manager");

			close("gH2AX-1"); close("ctip2pos-nucl-1");
			close("Mask of Result of gH2AX-b");


	///Quantify gH2AX intensity in cytobody and nucleus for CTIP2- cells 

			/// Remove bright signal and measure intensity of gH2AX in CTIP2- nuclei & cyto ("ctip2neg-body-fin" mask)	
			selectWindow("gH2AX-b");
			imageCalculator("AND create","gH2AX-b","ctip2neg-nucl-fin");
			run("Analyze Particles...", "size=21-infinity circularity=0.6-1.00 show=Masks");
			run("Invert LUT");	
			run("Duplicate...", "title=ctip2neg-nucl-1"); close("gH2AX-b"); close("Result of gH2AX-b");
		  
			run("Set Measurements...", "area mean standard min feret's integrated limit redirect=None decimal=2");
			run("Analyze Particles...", "size=0-infinity pixel circularity=0-1 show=Outlines add in_situ");
			
			selectWindow("Result of gH2AX");
			run("Duplicate...", "title=gH2AX-1");
	        setThreshold(900,65500); 

	        if(roiManager("Count") > 0){
			run("From ROI Manager");
			}
			roiManager("Measure");
			saveAs("Results", dir1+"data/h2ax_tuj/"+name+"_h2ax_ctip2_NEGnucl_intensity_count.csv"); ///
			close(name+"_h2ax_ctip2_NEGnucl_intensity_count.csv"); ///
			run("Clear Results");close("Summary");close("ROI Manager");
			close("gH2AX-1"); close("ctip2neg-nucl-1");
			close("Mask of Result of gH2AX-b");


///////////////// close all remaining windows 
			close("ctip2pos-nucl-fin"); 
			close("ctip2neg-nucl-fin"); 
			close("ctip2pos-ori");
	
			close("TujNucl4");close("body"); 
			close("Hoescht"); close("proj"); close("ctip2"); close("gH2AX");
			
			run("Close All");

		}
		 
	}

}

run("Close All");

