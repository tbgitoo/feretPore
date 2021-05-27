import feretPore.tools.VersionIndicator;
import ij.IJ;
import ij.plugin.PlugIn;

public class VersionDisplay implements PlugIn {

	
	public void run(String arg) {
		IJ.showMessage("FeretPore","FeretPore_.jar version = "+VersionIndicator.versionJar);
		
	}

}
