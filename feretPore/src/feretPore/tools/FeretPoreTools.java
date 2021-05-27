package feretPore.tools;

import org.orangepalantir.leastsquares.fitters.MarquardtFitter;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.Blitter;
import ij.process.ByteProcessor;
import ij.plugin.ChannelSplitter;

public class FeretPoreTools {
	
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

}
