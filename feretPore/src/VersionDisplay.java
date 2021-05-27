import feretPore.tools.VersionIndicator;
import ij.IJ;
import ij.plugin.PlugIn;

public class VersionDisplay implements PlugIn {

	
	public void run(String arg) {
		IJ.showMessage("FeretPore","Feret.jar version = "+VersionIndicator.versionJar);
		
	}

}
