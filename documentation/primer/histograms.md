# Histograms #

Histograms play a fundamental role in any type of physics analysis, not
only to visualise measurements but being a powerful form of data
reduction. ROOT offers many classes that represent histograms, all
inheriting from the `TH1` class. We will focus in this chapter on uni-
and bi- dimensional histograms the bin contents of which are represented by
floating point numbers [^4], the `TH1F` and `TH2F` classes respectively.

## Your First Histogram ##

Let's suppose you want to measure the counts of a Geiger detector located in
proximity of a radioactive source in a given time interval. This would
give you an idea of the activity of your source. The count distribution
in this case is a Poisson distribution. Let's see how operatively you
can fill and draw a histogram with the following example macro.

``` {.cpp .numberLines}
@ROOT_INCLUDE_FILE macros/macro5.C
```

Which gives you the following plot (Figure [5.1](#f51)):

[f51]: figures/poisson.png "f51"
<a name="f51"></a>

![The result of a counting (pseudo) experiment. Only bins corresponding
to integer values are filled given the discrete nature of the poissonian
distribution. \label{f51}][f51]

Using histograms is rather simple. The main differences with respect to
graphs that emerge from the example are:

-   line *5*: The histograms have a name and a title right from the
    start, no predefined number of entries but a number of bins and a
    lower-upper range.

-   line *15*: An entry is stored in the histogram through the
    `TH1F::Fill` method.

-   line *18* and *21*: The histogram can be drawn also normalised, ROOT
    automatically takes cares of the necessary rescaling.

-   line *24* to *30*: This small snippet shows how easy it is to access
    the moments and associated errors of a histogram.

## Add and Divide Histograms ##

Quite a large number of operations can be carried out with histograms.
The most useful are addition and division. In the following macro we
will learn how to manage these procedures within ROOT.

``` {.cpp .numberLines}
@ROOT_INCLUDE_FILE macros/macro6.C
```

The plots that you will obtain are shown in Figures [5.2](#f52) and [5.3](#f53).

[f52]: figures/histo_sum.png "f52"
<a name="f52"></a>

![The sum of two histograms.\label{f52}][f52]

[f53]: figures/histo_ratio.png "f53"
<a name="f53"></a>

![The ratio of two histograms.\label{f53}][f53]

Some lines now need a bit of clarification:

-   line *3*: Cling, as we know, is also able to interpret more than one
    function per file. In this case the function simply sets up some
    parameters to conveniently set the line of histograms.

-   line *19* to *21*: Some `C++` syntax for conditional
    statements is used to fill the histograms with different numbers of
    entries inside the loop.

-   line *30*: The sum of two histograms. A weight, which can be negative, can
    be assigned to the added histogram.

-   line *41*: The division of two histograms is rather straightforward.

-   line *44* to *62*: When you draw two quantities and their ratios, it
    is much better if all the information is condensed in one single
    plot. These lines provide a skeleton to perform this operation.

## Two-dimensional Histograms ##

Two-dimensional histograms are a very useful tool, for example to
inspect correlations between variables. You can exploit the
bi-dimensional histogram classes provided by ROOT in a simple way.
Let's see how in this macro:

``` {.cpp}
@ROOT_INCLUDE_FILE macros/macro7.C
```

Two kinds of plots are provided within the code, the first one
containing three-dimensional representations (Figure [5.4](#f54)) and the second one
projections and profiles (Figure [5.5](#f55)) of the bi-dimensional histogram.

[f54]: figures/th2f.png "f54"
<a name="f54"></a>

![Different ways of representing bi-dimensional
histograms.\label{f54}][f54]

[f55]: figures/proj_and_prof.png "f55"
<a name="f55"></a>

![The projections and profiles of bi-dimensional
histograms.\label{f55}][f55]

When a projection is performed along the x (y) direction, for every bin
along the x (y) axis, all bin contents along the y (x) axis are summed
up (upper the plots of Figure [5.5](#f55)). When a profile is performed along the x (y)
direction, for every bin along the x (y) axis, the average of all the
bin contents along the y (x) is calculated together with their RMS and
displayed as a symbol with error bar (lower two plots of Figure [5.5](#f55)).

Correlations between the variables are quantified by the methods
`Double_t GetCovariance()` and `Double_t GetCorrelationFactor()`.

\newpage

## Multiple histograms ##

The class `THStack` allows to manipulate a set of histograms as a single entity.
It is a collection of `TH1` (or derived) objects. When drawn, the X and Y axis
ranges are automatically computed such as all the histograms will be visible.
Several drawing option are available for both 1D and 2D histograms. The next
macros shows how it looks for 2D histograms:

``` {.cpp .numberLines}
@ROOT_INCLUDE_FILE macros/hstack.C
```

- Line *4*: creates the stack.

- Lines *4-18*: create two histograms to be added in the stack.

- Lines *20-21*: add the histograms in the stack.

- Line *23*: draws the stack as a lego plot. The colour distinguish the two histograms [5.6](#f56).

[f56]: figures/hstack.png "f56"
<a name="f56"></a>

![Two 2D histograms stack on top of each other.\label{f56}][f56]

[^4]: To optimise the memory usage you might go for one byte (TH1C), short (TH1S), integer (TH1I), long64 (TH1L64) or double-precision (TH1D) bin-content.
