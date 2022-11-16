//To analyze all images in timecourse folder 

dir="E:/time_course-NGN2/CM-NGN2-example/D10/sting/Data/"; 

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
		
			//Run twice, switching which line is used; top for sites 1-9, bottom for sites>=10
			//if(matches(images[j],"[A-Z][0-9]{2}_w3_s[0-9]{1}_.tif")==1) {
			if(matches(images[j],"[A-Z][0-9]{2}_w2_s[0-9]{2}.tif")==1) {
			
			//print(images[j]); } } }
			open(dir1+"proj_max/"+images[j]);
			name=getTitle();
			run("Clear Results");
			run("Set Measurements...", "area limit redirect=None decimal=2");
			run("Duplicate...", "title=proj");
			run("Properties...", "channels=1 slices=1 frames=1 unit=um pixel_width=0.17 pixel_height=0.17 voxel_depth=0.17");

			//open Hoescht 
			image2=replace(images[j],"_w2_","_w1_");
			open(dir1+"proj_max/"+image2);
			run("Duplicate...", "title=Hoescht");
			run("Properties...", "channels=1 slices=1 frames=1 unit=um pixel_width=0.17 pixel_height=0.17 voxel_depth=0.17");

			//open sting
			image4=replace(images[j],"_w2_","_w4_");
			open(dir1+"proj_max/"+image4);
			run("Duplicate...", "title=sting");
			run("Properties...", "channels=1 slices=1 frames=1 unit=um pixel_width=0.17 pixel_height=0.17 voxel_depth=0.17");


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
	       		
			setThreshold(140, 65535); 
			run("Make Binary");
			run("Analyze Particles...", "size=0-infinity circularity=0-1.00 show=Masks");
			run("Invert");
			run("Duplicate...", "title=Hoescht-dead-1");
			run("Options...", "iterations=8 count=1 black do=Dilate");
			run("Invert");
			
			selectWindow("Hoescht-live");
			setThreshold(15, 65535); 
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
			
			run("Analyze Particles...", "size=50-infinity circularity=0.3-1.00 show=Masks");
			run("Duplicate...", "title=TujNucl4");
			run("Invert");
			run("Fill Holes");
			close("Result of Hoescht-live-1");close("Mask of Result of Hoescht-live-1");
			close("Hoescht-live-1"); close("Mask of Hoescht-live"); close("Mask of Hoescht-dead");
			close("Hoescht-dead-1");

///MAKE a mask of TUJ + Hoescht = TOTAL CELL to remove non specific TUJ cells
  // to remove bright small object in the TUJ signal
			selectWindow("proj");
			run("Duplicate...","title=projbis");
			run("Duplicate...","title=projbis1");
			run("Gaussian Blur...", "sigma=80");
			imageCalculator("Subtract create","projbis","projbis1");
			
	        setAutoThreshold("Intermodes dark");
			run("Make Binary"); close("projbis1");
			run("Options...", "iterations=3 count=1 black do=Dilate");
			
			selectWindow("Result of projbis"); run("Invert");
			imageCalculator("AND create","Result of projbis","projbis");
			close("Result of projbis");close("projbis");
			selectWindow("Result of Result of projbis");
			run("Duplicate...","title=proj1");
			close("Result of Result of projbis");

	//script to have cell bodies
			selectWindow("proj1");
			run("Duplicate...", "title=neurites_body");
			run("Gaussian Blur...", "sigma=6");
			run("Median...", "radius=0.3");
			setAutoThreshold("Moments dark");
			run("Convert to Mask"); run("8-bit");
			run("Options...", "iterations=7 count=1 black do=Erode");
			run("Options...", "iterations=5 count=1 black do=Dilate");
			run("Analyze Particles...", "size=70-infinity circularity=0.0-1.00 show=Masks");
			run("Invert"); run("Fill Holes");
			run("Analyze Particles...", "size=70-infinity circularity=0-1.00 show=Masks");
		  	close("neurites_body");	close("Mask of neurites_body");

	//to add live Hoescht nuclei to have cytoplasmic and nuclei signal
			selectWindow("TujNucl4");
			run("Invert"); run("8-bit");
			imageCalculator("Add create","Mask of Mask of neurites_body","TujNucl4");
			run("Invert LUT");
			run("Analyze Particles...", "size=90-infinity circularity=0-1 show=Masks");
			
			//final Hoescht: remove nuclei without tuj signal
			run("Invert");
			imageCalculator("And create","Mask of Result of Mask of Mask of neurites_body","TujNucl4");
			run("Invert");
			run("Options...", "iterations=1 count=1 black do=Dilate");
			imageCalculator("And create","Result of Mask of Result of Mask of Mask of neurites_body","TujNucl4");
			close("Result of Mask of Result of Mask of Mask of neurites_body");
			close("Result of Mask of Mask of neurites_body");
			close("TujNucl4");
			
	// Quantification of nuclei (count final Hoescht+ Tuj+ cells)
			selectWindow("Result of Result of Mask of Result of Mask of Mask of neurites_body");
			run("Duplicate...", "title=TujNucl4");
			run("Invert");
			run("Analyze Particles...", "size=0-infinity circularity=0-1 summarize");
			saveAs("Results", dir1+"data/Hoescht/"+name+"_hoescht_count.csv"); ///
			close(name+"_hoescht_count.csv"); ///
			run("Clear Results");close("Summary");
			close("Result of Result of Mask of Result of Mask of Mask of neurites_body");

			selectWindow("Mask of Result of Mask of Mask of neurites_body");
			run("Options...", "iterations=2 count=1 black do=Dilate");
			run("Fill Holes"); 
			run("Invert"); run("Invert LUT"); 

			selectWindow("proj1");
			run("Duplicate...", "title=neurites2");
			run("Gaussian Blur...", "sigma=5");
			run("Subtract Background...", "rolling=5");
			setAutoThreshold("Huang dark"); run("Make Binary");
			run("Analyze Particles...", "size=10-infinity circularity=0-1 show=Masks");
			run("Invert LUT"); run("Options...", "iterations=6 count=1 black do=dilate");
			imageCalculator("Substract create","Mask of Result of Mask of Mask of neurites_body","Mask of neurites2");
			
			run("Fill Holes"); run("Options...", "iterations=2 count=1 black do=erode");
			run("Analyze Particles...", "size=120-infinity circularity=0-1 show=Masks");
			run("Invert LUT");
			run("Options...", "iterations=2 count=1 black do=dilate"); run("Fill Holes");
			run("Options...", "iterations=3 count=1 black do=erode");
			
			run("Duplicate...", "title=body");
			run("Analyze Particles...", "size=0-infinity circularity=0-1 summarize");
			saveAs("Results", dir1+"data/area_tuj/"+name+"_Area_Tuj_body_count.csv"); ///
			close(name+"_Area_Tuj_body_count.csv"); ///
			run("Invert LUT");
			close("Mask of Result of Mask of Result of Mask of Mask of neurites_body");
			close("Result of Mask of Result of Mask of Mask of neurites_body");
			close("Mask of neurites2"); close("neurites2");
			close("Mask of Result of Mask of Mask of neurites_body");
			close("Mask of Mask of neurites_body");			
			run("Clear Results");close("Summary");
		
//neurites 
	//duplicate raw image and identify all neurites
			selectWindow("proj");
			run("Duplicate...", "title=neurites1");
			run("Gaussian Blur...", "sigma=0.7"); 
			run("Subtract Background...", "rolling=2");
			setThreshold(20,1000); 
			run("Convert to Mask");
		
			run("Options...", "iterations=1 count=4 black do=dilate");
			run("Options...", "iterations=2 count=4 black do=close");
			run("Duplicate...", "title=neurites_full");
			run("Analyze Particles...", "size=75-Infinity pixel circularity=0.00-0.9 show=Masks"); 			
			run("Invert LUT");
			
			imageCalculator("Substract create","Mask of neurites_full","body");

			selectWindow("Result of Mask of neurites_full");
			run("Duplicate...", "title=neurite");
			run("Invert");
			close("Result of Mask of neurites_full");
		    close("Mask of neurites_full");close("neurites_full");close("neurites1");

	// final NEURITE image
			selectWindow("neurite");
			run("Invert");
			run("Make Binary");
			run("Duplicate...", "title=neurites");
			run("Measure");
			saveAs("Results", dir1+"data/area_tuj/"+name+"_Area_Tuj_neurite_count.csv"); ///
			saveAs("PNG", dir1+"data/masks/"+name+"_full.png"); ///
			close(name+"_Area_Tuj_neurite_count.csv"); ///
			run("Clear Results");close("Summary");close("neurites");close("Results");

	//combine NEURITE and cell body
			selectWindow("body");
			run("Invert LUT");
			imageCalculator("Add create","neurite","body");
			
			selectWindow("Result of neurite");
			run("Duplicate...", "title=neurite_body_full");
			close("Result of neurite");
			close("proj1");


/////// STING Analysis: 
	/// Quantify STING (number/intensity) in cytobody, nucleus (to subtract in R), neurites,
	/// and perinuclear space

		/// Measure intensity of STING in nucleus + cytoplasm ("body" mask)	
			selectWindow("body"); 
	        run("Duplicate...","title=body-1");
			run("Set Measurements...", "area mean standard min feret's integrated limit redirect=None decimal=2");
			run("Analyze Particles...", "size=0-infinity pixel circularity=0-1 show=Outlines add in_situ");

			selectWindow("sting");
			run("Duplicate...","title=sting-a");
			run("Gaussian Blur...", "sigma=100");
			imageCalculator("Subtract create","sting","sting-a");
			run("Duplicate...","title=sting-1");
			setThreshold(1500,65500); close("sting-a");

			if(roiManager("Count") > 0){
			run("From ROI Manager");
			}
			roiManager("Measure");
			saveAs("Results", dir1+"data/sting_tuj/"+name+"_STING_Tuj_Total_intensity_count.csv"); ///
			close(name+"_STING_Tuj_Total_intensity_count.csv"); ///
			run("Clear Results");close("Summary");
			close("sting-1");close("ROI Manager");close("body-1");

		/// Measure number of STING particles in nucleus + cytoplasm 
			selectWindow("Result of sting"); 
			run("Duplicate...", "title=sting-2");
			setThreshold(1500,65500);
			run("Make Binary");
			run("Invert");
			selectWindow("body"); 
	        imageCalculator("Substract create","body","sting-2");
			run("Analyze Particles...", "size=0-infinity pixel circularity=0-1 summarize");
            saveAs("Results", dir1+"data/sting_tuj/"+name+"_STING_Tuj_Total_particle_count.csv"); ///
			close(name+"_STING_Tuj_Total_particle_count.csv"); ///
			run("Clear Results");close("Summary");
			close("ROI Manager");close("Result of body");

	/// Measure intensity of perinuclear STING (dilated "TujNucl4" mask)		
		selectWindow("TujNucl4"); 
		run("Duplicate...","title=TujNucl4-2");
		run("Options...", "iterations=12 count=1 black do=Dilate");
		run("Watershed");
		run("Set Measurements...", "area mean standard min feret's integrated limit redirect=None decimal=2");
		run("Analyze Particles...", "size=7-infinity pixel circularity=0-1 show=Outlines add in_situ");

		selectWindow("Result of sting");
		run("Duplicate...","title=sting-1"); 
		setThreshold(1500,65500); 

		if(roiManager("Count") > 0){
		run("From ROI Manager");
		}
		roiManager("Measure");
		saveAs("Results", dir1+"data/sting_tuj/"+name+"_STING_Tuj_Peri_Nuclear_intensity_count.csv"); ///
		close(name+"_STING_Tuj_Peri_Nuclear_intensity_count.csv"); ///
		run("Clear Results");close("Summary");
		close("ROI Manager");close("sting-1"); close("TujNucl4-2");

	/// Measure intensity of STING in neurites ("neurite" mask)
			selectWindow("Result of sting");
			run("Duplicate...","title=sting-b"); run("8-bit");
			setThreshold(180,255);
			run("Make Binary");
			run("Options...", "iterations=10 count=1 black do=Dilate");
			run("Invert");
			imageCalculator("AND create","sting-b","neurite");
			run("Duplicate...", "title=neurite-1");  close("sting-b");
			//close("Result of sting-b");
		
			run("Set Measurements...", "area mean standard min feret's integrated limit redirect=None decimal=2");
			run("Analyze Particles...", "size=7-infinity pixel circularity=0-1 show=Outlines add in_situ");
			
			selectWindow("Result of sting");
			run("Duplicate...", "title=sting-1");
	        setThreshold(1500,65500); 

	        if(roiManager("Count") > 0){
			run("From ROI Manager");
			}
			roiManager("Measure");
			saveAs("Results", dir1+"data/sting_tuj/"+name+"_STING_Tuj_Neurite_intensity_count.csv"); ///
			close(name+"_STING_Tuj_Neurite_intensity_count.csv"); ///
			run("Clear Results");close("Summary"); close("ROI Manager");
			close("sting-1"); 

	/// Measure number of STING particles in neurites		
			selectWindow("sting-2"); 
	        imageCalculator("Substract create","Result of sting-b","sting-2");
			run("Analyze Particles...", "size=7-infinity pixel circularity=0-1 summarize");
            saveAs("Results", dir1+"data/sting_tuj/"+name+"_STING_Tuj_NEURITE_particle_count.csv"); ///
			close(name+"_STING_Tuj_NEURITE_particle_count.csv"); ///
			run("Clear Results");close("Summary");
			close("ROI Manager"); close("Result of Result of sting-b"); close("Result of sting");
			close("sting-2"); close("neurite-1"); close("Result of sting-b");

			close("Hoescht"); close("proj"); close("H2AX"); close("sting");
			close("TujNucl4"); close("Hoescht-dead-1"); close("body"); 
			close("neurite"); close("neurite_body_full");
			
			run("Close All");

		}
		 
	}

}

run("Close All");

