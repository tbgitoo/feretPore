import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.awt.Scrollbar;

import feretPore.tools.VersionIndicator;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.PlugInFilter;
import ij.process.AutoThresholder;
import feretPore.tools.FeretPoreTools;
import ij.process.ImageProcessor;

// The basic idea of this plugin is that one draws straight
// lines through the image and then quantifies the length
// of intersections with the pore space

public class FeretPore implements PlugInFilter,DialogListener {

	// image to work on
	protected ImagePlus imp_greyscale;
	
	protected ImagePlus imp_rgb=null;

	// image processor at the time of starting the analysis
	protected ImageProcessor ip;

	

	// Possible choices are 
	// 0 => "Summary statistics only"
	// 1 => "Raw list of intersection lengths"
	// 2 => "Histogram"
	public static String[] outputChoices;

	// Among the options above, the one that we should actually do
	public static String chosenOutput=null;

	// Possibility to impose a minimal pore length
	public static double minLength=0;

	// Possibility to impose a minimal pore length as a fraction
	// of mean pore length (self consistant solution)
	public static double fraction_min_length=0.1;

	// Minimal value for a pixel to be considered wall space rather than
	// pore space
	public static double threshold = 127;

	// Possibility to calculate threshold via autothreshold
	public static boolean useAutoThreshold = false;

	// Number of random lines to be drawn across the image
	public static int n_lines=10000;


	// Flag to indicate whether a self consistant solution with 
	// a minimal pore length as a given fraction of the mean should be sought
	public static boolean adjust_min_length_as_fraction_of_mean=false;

	// Holds the list intersections with the pore
	public double[] length_list=new double[0];
	

	// For each intersection, indicates to which line it belongs to
	public int[] line_list=new int[0];

	// This function updates the class variables while the 
	// dialog is shown
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {

		// Intermediate variable to accept numbers before filtering
		double n;

		// Get the first number from the dialog
		n = gd.getNextNumber();
		// Do basic checking, should be a valid number
		if (gd.invalidNumber())
			return false;
		// More checking, should be an integer, and > 0
		if(Math.round(n)>0)
		{
			n_lines=(int) Math.round(n);
		}

		useAutoThreshold = gd.getNextBoolean();

		// Get next number, should be strictly positive integer
		n=gd.getNextNumber();
		if(Math.round(n)>0)
		{
			threshold=(int) Math.round(n);
		}



		// Get next number, minimal lenth can be zero or larger
		n=gd.getNextNumber();
		if(Math.round(n)>=0)
		{
			minLength=Math.round(n);
		}

		// Fix the type of output chosen
		chosenOutput=gd.getNextChoice();


		// Read whether the self consistant solution with a cutoff 
		// being a fraction of the mean should be sought
		adjust_min_length_as_fraction_of_mean=gd.getNextBoolean();

		// Input of the fraction of the mean for the cutoff in 
		// case of the self consistant solution
		n=gd.getNextNumber();
		if(n>=0)
		{
			fraction_min_length=n;
		}
		
		

		return true;
	}

	// Constructor, initialize the choices and the initial choice
	public FeretPore()
	{
		outputChoices = new String[3];
		outputChoices[0] = "Summary statistics only";
		outputChoices[1] = "Raw list of intersection lengths";
		outputChoices[2] = "Histogram";
		if(chosenOutput==null)
		{
			chosenOutput=outputChoices[0];
		}

	}

	// Mandatory function for use in ImageJ, via the interface PlugInFilter
	public int setup(String arg, ImagePlus imp) {

		if (imp.getProcessor().getNChannels()==1) // Greyscale
		{
			this.imp_greyscale = imp;
			this.imp_rgb=null;
		} else { // RGB
			this.imp_greyscale=FeretPoreTools.greyFromMaxRGB(imp);
			this.imp_rgb=imp;
					
		}
		// Store an internal reference to the assigned image
		

		// For now, restricted to 8bit greyscale images
		return DOES_RGB+DOES_8G+NO_CHANGES;
	}

	// Mandatory function for use in ImageJ, via the interface PlugInFilter
	// This function actually runs the plugin
	public void run(ImageProcessor theIp) {

		// store the ImageProcessor internall
		ip = theIp;

		// Register this class to avoid garbage collection
		// usually done in ImageJ plugins, not of necessarily of
		// demonstrated usefulness
		IJ.register(this.getClass()); 

		// Show the dialog for choosing the options
		if(!doDialog())
		{
			return;
		}



		// Get the general results table
		
		ResultsTable	rt = Analyzer.getResultsTable();
		

		// if autothresholding is requested, we need to use the
		// corresponding imageJ plugin to get the autothreshold values
		if(useAutoThreshold)
		{
			threshold=getAutoThreshold(ip);
		}

		// The user wants the self-consistant solution, so
		// run the iterative procedure to find the parameters for
		// this
		if(adjust_min_length_as_fraction_of_mean)
		{
			adjust_for_fraction();
		}



		// Do the actual pore size analysis
		do_pore_size_analysis();

		// Write the results into the results table
		outputResultsToTable(rt);

		
		

		// Make sure small numbers are displayed correctly

		rt.setPrecision(-4);

		// Show the results table
		rt.show("Results");








	}
	
	// Calculates an autothreshold to distinguish background and foreground
	// The algorithm is as follows:
	// on the 80% least bright pixels, estimate a mean and standard deviation
	// on a log-scaled histogram
	// From this, evaluate the boundraries for a symmetric interval of confidence
	// of 99%
	// Compare this to 5% of the maximum intensity
	// Whatever number is higher, take as the threshold
	// The idea is that if there is a strong background, we want to 
	// avoid counting background pixels as foreground by mistake
	// but if the background is very clean, the main danger is to attribute pixels
	// to the foreground, even though their brightness mainly results from 
	// out-of-plane contribution due to the finite size of the 
	// confocal detection volume, particularly in z-direction

	public static int getAutoThreshold(ImageProcessor ip)
	{
		
		
		
		int[] histo_int=ip.getHistogram();
		
		
		double[] histo_double = new double[histo_int.length];
		
		for(int ind=0; ind<histo_double.length; ind++)
		{
			histo_double[ind] = histo_int[ind];
		}

		int cutoff_background=(int)Math.ceil(FeretPoreTools.getQuantile(histo_double, 0.8));

		if(cutoff_background <0 )
		{
			cutoff_background=0;
		}
		if(cutoff_background > 255)
		{
			cutoff_background=255;
		}

		double[] bg_histo = new double[cutoff_background];

		System.arraycopy(histo_double, 0, bg_histo, 0, cutoff_background);

		double[] bg_histo_for_mean = new double[bg_histo.length];

		for(int ind=0; ind<bg_histo.length; ind++)
		{
			bg_histo_for_mean[ind]=bg_histo[ind]*((double) ind);
		}

		double bg_mean_guess = FeretPoreTools.sum(bg_histo_for_mean);

		double[] bg_histo_for_sd = new double[bg_histo.length];

		for(int ind=0; ind<bg_histo.length; ind++)
		{
			bg_histo_for_sd[ind]=bg_histo[ind]*((double) ind-bg_mean_guess)*((double) ind-bg_mean_guess);
		}

		double bg_sd_guess = Math.sqrt(FeretPoreTools.sum(bg_histo_for_sd));

		// For fitting, use 2x the mean
		
		int to_use = (int)Math.min((int)Math.round(bg_mean_guess*2), histo_double.length);
		
		bg_histo = new double[to_use];

		System.arraycopy(histo_double, 0, bg_histo, 0, to_use);
		
		
		

		double[] m_sd= FeretPoreTools.gaussian_fit_histogram_log(bg_histo, bg_mean_guess, bg_sd_guess); 
		
		// 99% interval:
		// Fromt the normal distribution
		double quantile_99p_symmetric=2.575829;
		
		double upper = m_sd[0]+m_sd[1]*quantile_99p_symmetric;
		
		
		int limit_from_background=(int)Math.ceil(upper);
		
		// Reasonably, we should also be above 0.1*the highest intensity
		// This condition is to avoid taking out-of-plane pixels as foreground in the
		// presence of a very weak background
		
		int max=0;
		for(int ind=0; ind<histo_double.length; ind++)
		{
			if(histo_double[ind]>0)
			{
				max=ind;
			}
		}
		
		
		
		return Math.max((int)Math.ceil((double)max*0.1), limit_from_background);

		

	}



	// Helper function: weighted mean of the values in an array
	// length of x and length of weights needs to be identical
	public static double mean(double [] x, double [] weights)
	{
		if(weights==null)
		{
			return mean(x);
		}
		double s = 0;
		double sum_weights=0;
		for(int ind=0; ind<x.length; ind++)
		{
			s = s + x[ind]*weights[ind];
			sum_weights = sum_weights+weights[ind];
		}
		if(sum_weights>0)
		{
			return (s/sum_weights);
		}
		return 0;

	}


	// Helper function: Calculates the arithmetic mean of the values in the list
	public static double  mean(double [] x)
	{
		double s = 0;
		for(int ind=0; ind<x.length; ind++)
		{
			s = s + x[ind];
		}
		if(x.length>0)
		{
			return (s/(double)x.length);
		}
		return 0;

	}

	// Calculates the sum of the values in the list
	public static double  sum(double [] x)
	{
		double s = 0;
		for(int ind=0; ind<x.length; ind++)
		{
			s = s + x[ind];
		}
		if(x.length>0)
		{
			return (s);
		}
		return 0;

	}


	// Calculates a weighted standard deviation; the idea here is that if all the weights
	// are equal to 1, this reduces to the standard N-1 formula
	public static double standard_deviation(double [] x, double [] weights)
	{

		if(weights==null)
		{
			return standard_deviation(x);
		}
		double s = 0;
		double sum_weights=0;

		double m = mean(x,weights);

		for(int ind=0; ind<x.length; ind++)
		{
			s = s + weights[ind]*(x[ind]-m)*(x[ind]-m);
			sum_weights = sum_weights+weights[ind];
		}
		if(x.length>1)
		{
			// If all the weights are 1, or in fact equal, 
			//  we should return the standard N-1 formula
			return Math.sqrt(s/(sum_weights*(double)(x.length-1)/(double)x.length));
		}
		return 0;

	}

	// Calculates the standard of the values in the list, uses
	// the n-1 formula
	public static double standard_deviation(double [] x)
	{
		double s = 0;
		for(int ind=0; ind<x.length; ind++)
		{
			s = s + x[ind]*x[ind];
		}
		if(x.length>1)
		{
			double m=mean(x);
			return Math.sqrt((s-(double)x.length*m*m)/((double)x.length-1));
		}
		return 0;

	}

	// Self consistent solution:
	// Tries to find a cutoff such that minimum length is a 
	// defined fraction (fraction_min_length) of the mean pore size
	// Implementation: With a first guess (entered by user) for the cutoff, evaluate
	// the mean pore size (number weighted), then take the desired fraction of the mean
	// pore size as the new guess, and repeat 3x
	public void adjust_for_fraction()
	{
		for(int ind=1; ind <=3; ind++)
		{
			do_pore_size_analysis();

			double m = mean(length_list);

			minLength = m*fraction_min_length;
		}

	}



	// Selects the lines, and compiles the final length list, stores it internally (in length_list and 
	public void do_pore_size_analysis()
	{

		// Start with empty list of lengths ...
		double[] final_list = new double[0];
		// and empty index list (the list which indicates the line to 
		// which each fragment belongs)
		int[] final_index_list = new int[0];

		// Loop: produce crossecting line, collect the intersecting segments
		for(int indexLine=1; indexLine <= n_lines; indexLine++)
		{
			// random line
			Line theLine=selectRandomLine();


			// get the pore sizes, in pixels
			double[] list=pore_size(theLine.getPixels());

			// get the size of the pixels from the known length of the line and its raw length in pixels
			double pixelSize = theLine.getLength()/theLine.getRawLength();

			// multiply the pore sizes in pixels with the pixelSize 
			// to get the true length
			for(int ind=0; ind<list.length; ind++)
			{
				list[ind]=list[ind]*pixelSize;
			}

			// Go through the length list to keep only the segments that
			// are long enough (larger than minLength)
			for(int ind=0; ind<list.length; ind++)
			{
				if(list[ind]>=minLength)
				{
					double[] new_final_list=new double[final_list.length+1];
					System.arraycopy(final_list, 0, new_final_list, 0, final_list.length);
					final_list = new_final_list;
					final_list[final_list.length-1]=list[ind];

					int[] new_final_index_list=new int[final_index_list.length+1];
					System.arraycopy(final_index_list, 0, new_final_index_list, 0, final_index_list.length);
					final_index_list = new_final_index_list;
					final_index_list[final_index_list.length-1]=indexLine;

				}
			}

		}


		// Update the class variables for the length of the segments and the line identity
		length_list=final_list;

		line_list = final_index_list;





	}







	// Returns a list of pore lengths in pixels by analysis of the pixels of 1 line
	// It's the lengths of the enclosed pore space between wall pixels
	public double[] pore_size(double[] pixels)
	{



		// We start the first pore after having had wall pixels
		boolean poreStarted=false;
		// For counting the pores
		int n=0;



		// Well we'll have to reallocate each time since Java doesn't provide decent arrays
		double[] poreSizeList = new double[n];

		n=0; // no pores identified yet

		double currentPoreLength = 0;

		// By linear interpolation, the fraction of a pixel required to reach the 
		// threshold
		double start_fraction = 0;

		// By linear interpolation, the fraction of a pixel required to reach the threshold at the end of the pore
		double stop_fraction = 0;

		boolean first_wall_reached = false;

		for(int ind =0; ind<pixels.length; ind++)
		{
			double theValue = pixels[ind];

			// in a pore
			if(theValue < threshold)
			{
				// We should only count the pore if we had already a wall before, otherwise, this
				// measures the distance to the edge of the image instead
				if(first_wall_reached)
				{
					if(!poreStarted) // Freshly in a pore, so initiate the length
					{
						// Since we have already passed a wall, it cannot be the first pixel

						double x0 = 0;
						double x1 = 1;
						double y0 = pixels[ind-1];
						double y1 = pixels[ind];

						start_fraction=(x1-x0)*(y1-threshold)/(y1-y0);



						poreStarted=true;
						currentPoreLength=1;
					} else // Already running inside a pore, increase length
					{

						currentPoreLength++;

					}
				}
			} else // in the wall
			{
				// If we had started a true pore before, then
				// reaching teh wall means completion of the pore
				if(poreStarted) // This is a complete pore ...
				{

					double x0 = 0;
					double x1 = 1;
					double y0 = pixels[ind-1];
					double y1 = pixels[ind];

					stop_fraction=(x1-x0)*(y1-threshold)/(y1-y0);

					// increase the length of the 
					// segment list and 
					//add an element	
					double[] newPoreSizeList=new double[n+1];
					System.arraycopy(poreSizeList, 0, newPoreSizeList, 0, poreSizeList.length);
					poreSizeList = newPoreSizeList;
					poreSizeList[n]=currentPoreLength+start_fraction-stop_fraction;


					// reset the pore length
					currentPoreLength=0;
					// increment pore counter
					n++;
					// reset the flag indicated we stared a pore
					poreStarted=false;
				}
				// In any case, running in the wall means that we 
				// can set the flag indicating we reached a wall
				first_wall_reached = true;
			}


		}



		return poreSizeList;

	}


	// Selection of random line (I-Randomness according to 
	// Coleman R., Metallography 6, 103-114 (1973)

	// The idea is to first select a random point in the interior
	// of the image (random x, random y, both uniform on the size of the image)
	// and to place a line with random direction (uniform angle between 0 and pi)

	public Line selectRandomLine()
	{
		// First, determine the anchorage point ("radiator")
		double radiator_x= Math.floor(Math.random()*imp_greyscale.getWidth());
		double radiator_y= Math.floor(Math.random()*imp_greyscale.getHeight());

		// Then, the random angle
		double angle = Math.random()*Math.PI;

		// Directional vector corresponding to the angle
		double dx = Math.cos(angle);
		double dy = Math.sin(angle);

		// Actual start point, needs to be calculated to be on periphery of the image
		int start_x=0;
		int start_y=0;

		// Actual end point, needs to be calculated to be on peripher of the image
		int end_x=0;
		int end_y=0;

		// Flag indicating whether we already got start and end point or whether we need to continue
		boolean done = false;

		// First trivial case: horizontal line
		if(dy==0)
		{
			start_x = 0;
			start_y = (int) radiator_y;

			end_x = imp_greyscale.getWidth()-1;
			end_y = (int) radiator_y;

			done=true;

		}
		// Second trivial case: vertical line
		if(dx==0 && !done)
		{
			start_x = (int) radiator_x;
			start_y = 0;

			end_x = (int) radiator_x;
			end_y = imp_greyscale.getHeight()-1;

			done=true;
		}
		// General case: Some oblique direction
		if(!done)
		{
			// Try: Putative start x coordinate at 0 , end coordinate at image width
			double sx=0;
			double ex=imp_greyscale.getWidth();

			// Calculation of associated y values by extrapolation with the direction vector dx,dy
			double sy=radiator_y+dy/dx*(sx-radiator_x);
			double ey=radiator_y+dy/dx*(ex-radiator_x);

			// Case where sy falls below the 0 line of the image
			if(Math.round(sy)<0)
			{
				// in that case the correct start y is 0
				sy=0;
				// and we need to calculate x with the direction vector
				sx=radiator_x + dx/dy*(sy-radiator_y);

			}

			// other case where adjustment is needed: start y beyond image height
			if(Math.round(sy)>(imp_greyscale.getHeight()-1))
			{
				// in that case, put start y at highest possible value
				sy=imp_greyscale.getHeight()-1;
				// and calculate start x from the direction vector
				sx=radiator_x + dx/dy*(sy-radiator_y);

			}


			// Same analysis for putative end y (ey)
			// Case where it falls below 0
			if(Math.round(ey)<0)
			{
				// Zero is the correct value
				ey=0;
				// and x needs to be estimated
				ex=radiator_x + dx/dy*(ey-radiator_y);

			}

			// sy falls beyond the total image height
			if(Math.round(ey)>(imp_greyscale.getHeight()-1))
			{
				// highest value line is the correct value for ey
				ey=imp_greyscale.getHeight()-1;
				// and again ex needs to be estimated
				ex=radiator_x + dx/dy*(ey-radiator_y);

			}


			// We can now definitely define the start and end point of the line

			start_x=(int)Math.round(sx);

			start_y=(int)Math.round(sy);

			end_x=(int)Math.round(ex);

			end_y=(int)Math.round(ey);


		}





		// Construct the graphical selection (imageJ line)

		Line theRoi = new Line(start_x,start_y,end_x,end_y);

		// To extract the pixel values along the line, we need to 
		// associate the line with the ImagePlus object of this plugin
		theRoi.setImage(imp_greyscale);

		return theRoi;


	}


	public boolean doDialog()
	{

		// Open a dialog to get the use variables
		// As a particular feature of ImageJ, does not open in macro mode but
		// is substituted with macro parameter values instead
		GenericDialog gd = new GenericDialog("Feret Pore diameter plugin");
		
		
		// Add the fields
		gd.addNumericField("Number of lines to be analyzed", n_lines, 0);

		gd.addCheckbox("Use autothreshold", useAutoThreshold);

		gd.addNumericField("Threshold for wall pixels", threshold, 0);



		gd.addNumericField("Minimal length of pore intersection", minLength, 2);


		gd.addChoice("Output", outputChoices, chosenOutput);

		gd.addCheckbox("Find miminal length as fraction of mean", adjust_min_length_as_fraction_of_mean);


		gd.addNumericField("Target min length as fraction of mean pore size", fraction_min_length, 3);

		// We need to follow the dialog to update the class variables
		gd.addDialogListener(this);
		// Show the dialog
		
		
		
		gd.showDialog();                    // input by the user (or macro) happens here
		// Do not proceed when the use pushes cancel
		
		
		
		return (!gd.wasCanceled());



	}


	public void outputResultsToTable(ResultsTable rt)
	{

		// The chosen output is the raw list
		if(chosenOutput.equals(outputChoices[1]) )
		{
			// Run through the raw list of intersections
			for(int ind=0; ind<length_list.length; ind++)
			{

				// Add new line in the results table to hold the 
				// results for this intersection
				rt.incrementCounter();
				// To which line does this interseciton below
				rt.addValue("Line", line_list[ind]);
				// Pore intersection index
				rt.addValue("Pore #", ind);
				// Length of the intersection
				rt.addValue("Section length", length_list[ind]);
				// Threshold value
				rt.addValue("Threshold", threshold);

				rt.addValue("Version FeretPore_.jar", VersionIndicator.versionJar);
				
			}

		}

		// Case where the user only wants the summary statistics
		if(chosenOutput.equals(outputChoices[0])){

			rt.incrementCounter();
			rt.addValue("Image", imp_greyscale.getTitle());

			// Mean pore size: Arithmetic mean of the intersections lengths; since 
			// small intersections are by definition more frequent if involving the
			// same number of pixels, this is "number" weighting
			rt.addValue("Mean pore size (Number weighted)", mean(length_list));

			// Mean pore size: Weighted mean, using the length of the
			// intersections as the weight. This is technically "length weighted"
			// but can also be considered "pixel weighted", as gives each intersection
			// the weight it has due to its numbers of pixels			
			rt.addValue("Mean pore size (Length weighted)", mean(length_list,length_list));

			// Standard deviations, with the two weighting methods
			rt.addValue("Std pore size (Number weighted)",standard_deviation(length_list));
			rt.addValue("Std pore size (Length weighted)",standard_deviation(length_list,length_list));

			// Final minimal pore length, this can be different from the user settings
			// in the case of the self-
			rt.addValue("Minimal pore length", minLength);
			// Number of intersection segments taken into account
			rt.addValue("N", length_list.length);
			rt.addValue("Threshold", threshold);
			
			rt.addValue("Version FeretPore_.jar", VersionIndicator.versionJar);

		}

		// Output a histogram. We need to collect the data first
		if(chosenOutput.equals(outputChoices[2]))
		{
			// The intersection lengths are not usually entire numbers
			// we collect them here to integer numbers in terms of the
			// unit chosen
			int[] rounded = new int[length_list.length];

			// For the histogram, the longest intersection
			int maxVal=0;

			// The shortest intersection
			int minVal=0;

			// Run through the intersections, find longest and shortest one
			for(int ind=0;ind<length_list.length;ind++)
			{
				rounded[ind]=(int) Math.round(length_list[ind]);
				if(ind==0)
				{
					maxVal=rounded[ind];
					minVal=rounded[ind];
				} else
				{
					if(rounded[ind]>maxVal)
					{
						maxVal=rounded[ind];
					}
					if(rounded[ind]<minVal)
					{
						minVal=rounded[ind];
					}
				}
			}


			// Given that we collect round lengths, we need
			// nbin = maxVal-minVal+1
			// We need to store the actual pore sizes
			int [] histogram_pore_size = new int[maxVal-minVal+1];
			// The raw frequency
			double [] histogram_frequency = new double[histogram_pore_size.length];
			// and the length-weighted frequency
			double [] histogram_length_frequency = new double[histogram_pore_size.length];

			// Initialize the histogram variables
			for(int ind=0; ind<histogram_pore_size.length; ind++)
			{
				histogram_pore_size[ind]=ind+minVal;
				histogram_frequency[ind]=0;
				histogram_length_frequency[ind]=0;
			}

			// run the intersection segments and distribute into the
			// histogram bins
			for(int ind=0;ind<rounded.length;ind++)
			{
				histogram_frequency[rounded[ind]-minVal]++;
				histogram_length_frequency[rounded[ind]-minVal]+=rounded[ind];


			}

			// Normalization such that the sum of the frequencies
			// is 1
			// Evaluate the current sums ...
			double hf=sum(histogram_frequency);
			double hfl=sum(histogram_length_frequency);

			// ... and divide the values by it 
			for(int ind=0;ind<histogram_frequency.length;ind++)
			{

				histogram_frequency[ind]=histogram_frequency[ind]/hf;
				histogram_length_frequency[ind]=histogram_length_frequency[ind]/hfl;

			}

			// Output the histogram, 1 line per entry in the results table
			for(int ind=0; ind<histogram_frequency.length; ind++)
			{


				rt.incrementCounter();
				rt.addValue("Pore size", histogram_pore_size[ind]);
				rt.addValue("Frequency number weighted", histogram_frequency[ind]);
				rt.addValue("Frequency length weighted", histogram_length_frequency[ind]);
				rt.addValue("Threshold", threshold);
				rt.addValue("Version FeretPore_.jar", VersionIndicator.versionJar);
			}



		}



	}


}
