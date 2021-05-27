# feretPore
Pore size determination by straight-line segments

This is an ImageJ plugin that determines pore size by placing randomly oriented straight lines across an image. Pore intersections are segments of these straight lines where the
the local pixel value remains below a threshold.

To run this plugin, download FeretPore_.jar from a release and place it in the plugins folder of your Fiji or ImageJ installation. Restart ImageJ/Fiji, and the plugin should appear
in the "Plugins" section of the menu, under "FeretPoreSize".

This plugin depends on the presence of the jama library. This library is usually provided in the jars folder of the Fiji or ImageJ application. It may however be necessary to 
download it and manually place it there ("jama-1.0.3.jar").

Credits are also given to Matthew B. Smith, whose code from https://github.com/odinsbane/least-squares-in-java/ we used here (MIT license).
