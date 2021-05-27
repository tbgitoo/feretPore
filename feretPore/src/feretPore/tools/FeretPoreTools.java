package feretPore.tools;

import org.orangepalantir.leastsquares.fitters.MarquardtFitter;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.Blitter;
import ij.process.ByteProcessor;
import ij.plugin.ChannelSplitter;

public class FeretPoreTools {
	
	public static String[] unique_string_values(String[] list)
	{
		String[] theUniqueList=new String[0];
		String[] newUniqueList=new String[0];
		for(int ind=0; ind<list.length; ind++)
		{
			if(theUniqueList.length>0)
			{
				boolean found=false;
				for(int listIndex=0; listIndex<theUniqueList.length; listIndex++)
				{
					if(theUniqueList[listIndex].equalsIgnoreCase(list[ind]))
					{
						found=true;
					}
				}
				if(!found)
				{
					newUniqueList=new String[theUniqueList.length+1];
					System.arraycopy(theUniqueList, 0, newUniqueList, 0, theUniqueList.length);
					newUniqueList[theUniqueList.length]=list[ind];
					theUniqueList=newUniqueList;
				}
				
			} else {
				theUniqueList=new String[1];
				theUniqueList[0]=list[ind];
			}
		}
		
		return theUniqueList;
		
	}
	
	// Light color combination
	public static String combine_main_color_string(String c1, String c2)
	{
		// Anything that contains white ends up white
		if(c1.equalsIgnoreCase("W")) {return "W"; }
		if(c2.equalsIgnoreCase("W")) {return "W"; }
		// identical arguments
		if(c1.equalsIgnoreCase(c2)) {return c1; }
		
		// Binary mixes with one of the components
		if(c1.equalsIgnoreCase("M") && c2.equalsIgnoreCase("B")) {return "M"; }
		if(c1.equalsIgnoreCase("M") && c2.equalsIgnoreCase("R")) {return "M"; }
		if(c2.equalsIgnoreCase("M") && c1.equalsIgnoreCase("B")) {return "M"; }
		if(c2.equalsIgnoreCase("M") && c1.equalsIgnoreCase("R")) {return "M"; }
		if(c1.equalsIgnoreCase("Y") && c2.equalsIgnoreCase("G")) {return "Y"; }
		if(c1.equalsIgnoreCase("Y") && c2.equalsIgnoreCase("R")) {return "Y"; }
		if(c2.equalsIgnoreCase("Y") && c1.equalsIgnoreCase("G")) {return "Y"; }
		if(c2.equalsIgnoreCase("Y") && c1.equalsIgnoreCase("R")) {return "Y"; }
		if(c1.equalsIgnoreCase("C") && c2.equalsIgnoreCase("G")) {return "C"; }
		if(c1.equalsIgnoreCase("C") && c2.equalsIgnoreCase("B")) {return "C"; }
		if(c2.equalsIgnoreCase("C") && c1.equalsIgnoreCase("G")) {return "C"; }
		if(c2.equalsIgnoreCase("C") && c1.equalsIgnoreCase("B")) {return "C"; }
		
		// Binary mixes with completion to white
		if(c1.equalsIgnoreCase("M") && c2.equalsIgnoreCase("G")) {return "W"; }
		if(c1.equalsIgnoreCase("M") && c2.equalsIgnoreCase("C")) {return "W"; }
		if(c1.equalsIgnoreCase("M") && c2.equalsIgnoreCase("Y")) {return "W"; }
		if(c2.equalsIgnoreCase("M") && c1.equalsIgnoreCase("G")) {return "W"; }
		if(c2.equalsIgnoreCase("M") && c1.equalsIgnoreCase("C")) {return "W"; }
		if(c2.equalsIgnoreCase("M") && c1.equalsIgnoreCase("Y")) {return "W"; }
		if(c1.equalsIgnoreCase("Y") && c2.equalsIgnoreCase("B")) {return "W"; }
		if(c1.equalsIgnoreCase("Y") && c2.equalsIgnoreCase("C")) {return "W"; }
		if(c1.equalsIgnoreCase("Y") && c2.equalsIgnoreCase("M")) {return "W"; }
		if(c2.equalsIgnoreCase("Y") && c1.equalsIgnoreCase("B")) {return "W"; }
		if(c2.equalsIgnoreCase("Y") && c1.equalsIgnoreCase("C")) {return "W"; }
		if(c2.equalsIgnoreCase("Y") && c1.equalsIgnoreCase("M")) {return "W"; }
		if(c1.equalsIgnoreCase("C") && c2.equalsIgnoreCase("R")) {return "W"; }
		if(c1.equalsIgnoreCase("C") && c2.equalsIgnoreCase("Y")) {return "W"; }
		if(c1.equalsIgnoreCase("C") && c2.equalsIgnoreCase("M")) {return "W"; }
		if(c2.equalsIgnoreCase("C") && c1.equalsIgnoreCase("R")) {return "W"; }
		if(c2.equalsIgnoreCase("C") && c1.equalsIgnoreCase("Y")) {return "W"; }
		if(c2.equalsIgnoreCase("C") && c1.equalsIgnoreCase("M")) {return "W"; }
		
		//Binary mixes
		if(c1.equalsIgnoreCase("R") && c2.equalsIgnoreCase("B")) {return "M"; }
		if(c1.equalsIgnoreCase("B") && c2.equalsIgnoreCase("R")) {return "M"; }
		if(c1.equalsIgnoreCase("R") && c2.equalsIgnoreCase("G")) {return "Y"; }
		if(c1.equalsIgnoreCase("G") && c2.equalsIgnoreCase("R")) {return "Y"; }
		if(c1.equalsIgnoreCase("G") && c2.equalsIgnoreCase("B")) {return "C"; }
		if(c1.equalsIgnoreCase("B") && c2.equalsIgnoreCase("G")) {return "C"; }
		
		
		return "combination?";
		
		
		
	}
	
	public static String main_color_string(int[] rgb_value)
	{
		if (rgb_value[0]>rgb_value[1] ) // R > G
		{
			if (rgb_value[2]>rgb_value[0]) // B > R > G
			{
				return "B"; // This is mainly blue
			}
			if (rgb_value[2]== rgb_value[0]) // B=R > G
			{
				return "M"; // This is magenta
			}
			return "R";
		}
		if (rgb_value[0]==rgb_value[1]) // R=G
		{
			if (rgb_value[2]>rgb_value[0])
			{
				return "B"; // This is again mainly blue
			}
			if (rgb_value[2]==rgb_value[0]) // R=G=B
			{
				return "W"; // For white
			}
			return "Y";
		}
		if (rgb_value[0]<rgb_value[1] )
		{
			if (rgb_value[2]>rgb_value[1])
			{
				return "B"; // This is again mainly blue
			}
			if (rgb_value[2]==rgb_value[1]) // R<G=B
			{
				return "C"; // For white
			}
			return "G";
		}
		return null;
	}
	
	public static int[] brightest_pixel_color(ImagePlus imp, double x, double y, double radius)
	{
		int[] result = new int[3];
		result[0]=0;
		result[1]=0;
		result[2]=0;
		
		int current_max = 0;
		
		double min_D=Math.abs(radius)*10;
		
		
		int search_x1 = (int) Math.round(x-Math.abs(radius)); 
		int search_y1 = (int) Math.round(y-Math.abs(radius));
		int search_x2 = (int) Math.round(x+Math.abs(radius)); 
		int search_y2 = (int) Math.round(y+Math.abs(radius));
		
		if(search_x1<0) {search_x1 = 0; }
		if(search_y1<0) {search_y1 = 0; }
		if(search_x2<0) {search_x2 = 0; }
		if(search_y2<0) {search_y2 = 0; }
		
		if(search_x1>=imp.getWidth()) {search_x1 = imp.getWidth()-1; }
		if(search_y1>=imp.getHeight()) {search_y1 = imp.getHeight()-1; }
		if(search_x2>=imp.getWidth()) {search_x2 = imp.getWidth()-1; }
		if(search_y2>=imp.getHeight()) {search_y2 = imp.getHeight()-1; }
		
		for(int xs=search_x1; xs<=search_x2; xs++)
		{
			for(int ys=search_y1; ys<=search_y2; ys++)
			{
				int[] thePixel = imp.getPixel(xs, ys);
				if ((thePixel[0]>0) | (thePixel[1]>0) | (thePixel[2]>0)) 
				{
					if (Math.max(Math.max(thePixel[0], thePixel[1]),thePixel[2])>current_max)
					{
					    current_max = Math.max(Math.max(thePixel[0], thePixel[1]),thePixel[2]);
					    min_D = Math.sqrt((x-xs)*(x-xs)+(y-ys)*(y-ys));
					    result[0]=thePixel[0];
					    result[1]=thePixel[1];
					    result[2]=thePixel[2];
					    
					    
					} else { // same maximum (maybe not same color, though), take the closer
						if (Math.max(Math.max(thePixel[0], thePixel[1]),thePixel[2])==current_max)
						{
							if (Math.sqrt((x-xs)*(x-xs)+(y-ys)*(y-ys))<min_D){
								 current_max = Math.max(Math.max(thePixel[0], thePixel[1]),thePixel[2]);
								 min_D = Math.sqrt((x-xs)*(x-xs)+(y-ys)*(y-ys));
								 result[0]=thePixel[0];
								 result[1]=thePixel[1];
								 result[2]=thePixel[2];
							}
						}
					}
				}
			}
			
		}
		
		
		
		
		
		return result;
	}
	
	public static double getQuantile(double [] hist, double p)
	{

		// Normalize (in case)
		// Also, calculate cumulative sums, with 1 element more than
		// the histogram to have 0 and 1 in it
		double total = 0;
		double[] cumsum = new double[hist.length+1];
		for(int ind=0; ind<hist.length; ind++)
		{
			total=total+hist[ind];
			cumsum[ind]=0;
		}


		for(int ind=0; ind<hist.length; ind++)
		{
			hist[ind]=hist[ind]/total;
			cumsum[ind+1]=cumsum[ind]+hist[ind];
		}
		// To be sure it's really exactly 1 and some close value due to
		// rounding errors
		cumsum[hist.length]=1;

		// limiting cases and non-treatable values
		if(p<=0) { return -1; }
		if(p>=1) { return hist.length; }

		// nominal case





		int current_ind=0;
		while(cumsum[current_ind]<p && current_ind<=hist.length)
		{

			current_ind++;
		}

		// Now we dispose of the element where the cumulative sum is just bigger than the
		// desired p value

		// Do interpolation to get a finer estimate

		double p_upper = cumsum[current_ind];
		double p_lower = cumsum[current_ind-1];

		double q_upper = current_ind;
		double q_lower = current_ind-1;

		// Degenerate case where the is no entry into the histogram here
		if(p_upper==p_lower)
		{
			int index_to_lower = current_ind-1;
			while(index_to_lower>0 && p_lower==p_upper)
			{
				index_to_lower--;
				p_lower = cumsum[index_to_lower];
				q_lower = index_to_lower;
			}

		}

		// This still hasn't helped, return the mean of the associated quantiles
		if(p_lower==p_upper)
		{
			return (q_lower+q_upper)/2;
		}

		double linear_inter_q = q_lower + (q_upper-q_lower)/(p_upper-p_lower)*(p-p_lower);

		return(linear_inter_q);




	}
	
	
	
	public static double[] gaussian_fit_histogram_log(double[] histogram, double mean_guess, double sd_guess)
	{

		// First of all, clean histogram from 0 values for which we cannot calculate a log
		int n_OK=0;
		for(int ind=0; ind<histogram.length; ind++)
		{
			if(histogram[ind]>0)
			{
				n_OK++;
			}
		}

		int[] indexes=new int[n_OK];
		double[] values=new double[n_OK];
		n_OK=0;
		for(int ind=0; ind<histogram.length; ind++)
		{
			if(histogram[ind]>0)
			{
				indexes[n_OK]=ind;
				values[n_OK]=Math.log(histogram[ind]);
				n_OK++;
			}
		}


		gaussianHistoLog theFun = new gaussianHistoLog();

		double[][] X = new double[values.length][1];



		for(int ind=0; ind<values.length; ind++)
		{
			X[ind][0]=indexes[ind];

		}

		double [] Z = values;

		for(int ind=0; ind<values.length; ind++)
		{
			Z[ind]=values[ind];

		}


		MarquardtFitter lft=new MarquardtFitter(theFun);
		lft.setData(X, Z);

		lft.setParameters(new double[]{mean_guess,Math.log(sd_guess),0});

		lft.fitData();

		double[] output = lft.getParameters();

		double[] final_output = new double[2];

		final_output[0]=output[0];
		final_output[1]=Math.exp(output[1]);





		return final_output;


	}
	
	public static ImagePlus greyFromMaxRGB(ImagePlus imp1)
	{
		ImageStack stack2=new ImageStack(imp1.getWidth(),imp1.getHeight());
		
		int nSlices=imp1.getImageStackSize();
		
		if (imp1.getProcessor().getNChannels()==1) // it's already greyscale, just
			// make sure it's returning 8bit
				{
				
				for(int i=1; i<=nSlices; i++) {
					String label = imp1.getStack().getSliceLabel(1);
					ImageProcessor ip = imp1.getStack().getProcessor(i);
					stack2.addSlice(label, ip.convertToByte(false));
				}
				
				
		
			
		} else { // it's rgb get maximum of colors
			for(int i=1; i<=nSlices; i++) {
				String label = imp1.getStack().getSliceLabel(i);
				
				ImageProcessor ip = new ByteProcessor(imp1.getWidth(),imp1.getHeight());
				
				ImagePlus RGB=new ImagePlus("RGB", imp1.getStack().getProcessor(i));
				for(int ic=1; ic<=imp1.getProcessor().getNChannels(); ic++)
				{
				
				  ip.copyBits(ChannelSplitter.getChannel(RGB,ic).getProcessor(1), 0,0,Blitter.MAX);
				
				}
						
						
				stack2.addSlice(label, ip);
			}
		}
		ImagePlus imp2 =new ImagePlus(imp1.getTitle(),stack2);
		return imp2;
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

}
