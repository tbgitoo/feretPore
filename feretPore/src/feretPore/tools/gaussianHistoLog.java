package feretPore.tools;

import org.orangepalantir.leastsquares.Function;

// Directly calculates the log of a Gaussian distribution
// parameters A and sd are also given in log

public class gaussianHistoLog implements Function {
    
    @Override
    public double evaluate(double[] values, double[] parameters) {
     // Get the variables
    	
    	double mean=parameters[0];
    	double sd_log=parameters[1];
    	double A_log=parameters[2];
    	
    	double x=values[0];
    	
    	
    	

    	return( A_log-sd_log-Math.log(Math.sqrt(2.0*Math.PI))-(x-mean)*(x-mean)/2/Math.exp(sd_log)/Math.exp(sd_log));
    	
    	
    }

    
    public int getNParameters() {
        return 3;
    }

    
    public int getNInputs() {
        return 1;
    }
}