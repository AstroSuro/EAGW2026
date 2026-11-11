# ZTF-lightcurves-extraction-from-IRAS
This is a simple example to extract ZTF light curves from IRAS.

This example will download and plot the ZTF light curves (in multiple bands) of ZTF19abanrhr (AGN J124942.3+344929). The light curves will compare with the trigger time of GW190521 (see Graham+ 2020, PRL).

The source position is at RA = 192.42625°, Dec = +34.82472°.

If the download from the notebook does not work due to time out, you can download the data with the following script.

<pre><code>
curl -o lc.csv "https<span>://</span>irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%20192.42625%2034.82472%200.000278&FORMAT=csv"</code></pre>
