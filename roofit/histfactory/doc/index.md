\defgroup HistFactory HistFactory
\brief Factory classes to create RooFit workspaces from histograms
\ingroup Roostats

<p>
This is a package that creates a RooFit probability density function from ROOT histograms 
of expected distributions and histograms that represent the +/- 1 sigma variations 
from systematic effects. The resulting probability density function can then be used
with any of the statistical tools provided within RooStats, such as the profile 
likelihood ratio, Feldman-Cousins, etc.  In this version, the model is directly
fed to a likelihood ratio test, but it needs to be further factorized.</p>

<p>
The user needs to provide histograms (in picobarns per bin) and configure the job
with XML.  The configuration XML is defined in the file `$ROOTSYS/config/HistFactorySchema.dtd`, but essentially
it is organized as follows (see the examples in `${ROOTSYS}/tutorials/histfactory/`)</p>

<ul>
<li> a top level 'Combination' that is composed of:</li>
<ul>
<li> several 'Channels' (eg. ee, emu, mumu), which are composed of:</li>
<ul>
<li> several 'Samples' (eg. signal, bkg1, bkg2, ...), each of which has:</li>
<ul>
<li> a name</li>
<li> if the sample is normalized by theory (eg N = L*sigma) or not (eg. data driven)</li>
<li> a nominal expectation histogram</li>
<li> a named 'Normalization Factor' (which can be fixed or allowed to float in a fit)</li>
<li> several 'Overall Systematics' in normalization with:</li>
<ul>
<li> a name</li>
<li> +/- 1 sigma variations (eg. 1.05 and 0.95 for a 5% uncertainty)</li>
</ul>
<li> several 'Histogram Systematics' in shape with:</li>
<ul>
<li> a name (which can be shared with the OverallSyst if correlated)</li>
<li> +/- 1 sigma variational histograms</li>
</ul>
</ul>
</ul>
<li> several 'Measurements' (corresponding to a full fit of the model) each of which specifies</li>
<ul>
<li> a name for this fit to be used in tables and files</li>
<li> what is the luminosity associated to the measurement in picobarns</li>
<li> which bins of the histogram should be used</li>
<li> what is the relative uncertainty on the luminosity </li>
<li> what is (are) the parameter(s) of interest that will be measured</li>
<li> which parameters should be fixed/floating (eg. nuisance parameters)</li>
</ul>
</ul>
</ul>

For more documentation see [CERN-OPEN-2012-016](https://cds.cern.ch/record/1456844/).

For tutorials see \ref tutorial_histfactory.
