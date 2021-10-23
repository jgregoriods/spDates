<h1>spDates</h1>
<h2>Analysis of Spatial Gradients in Radiocarbon Dates</h2>

**Jonas Gregorio de Souza**<br/>
jonas.gregorio@gmail.com<br/>
[![ORCiD](https://img.shields.io/badge/ORCiD-0000--0001--7879--4531-green.svg)](https://orcid.org/0000-0001-6032-4443)<br/>

<p>Tools to perform time-space regressions, an approach often used by archaeologists examining the expansion of cultural phenomena. In essence, one plots the radiocarbon ages of archaeological sites against their distances from a hypothetical origin. If a cultural advance has indeed taken place, the expectation is that radiocarbon dates will be more recent as one moves away from the center of origin. If a significant correlation is found, the intercept of the regression can be used as an estimate of the start date for the dispersal, while the regression slope provides an estimate of the speed of advance. Most applications have been focused on the Neolithic expansion from the Near East to Europe (Ammerman and Cavalli-Sforza 1971; <a href="https://doi.org/10.1017/S0003598X00061330">Gkiasta et al. 2003</a>; <a href="https://doi.org/10.1371/journal.pbio.0030410">Pinhasi et al. 2005</a>), but other case studies include the Paleolithic recolonization of Northern Europe (<a href="https://doi.org/10.1017/S0959774304000046">Fort et al. 2004</a>), the Clovis expansion in North America (<a href="https://doi.org/10.1073/pnas.0704215104">Hamilton and Buchanan 2007</a>), the human colonization of the Americas from Asia (<a href="https://doi.org/10.1371/journal.pone.0012472">Hamilton and Buchanan 2010</a>), the Lapita spread in Austronesia (<a href="https://doi.org/10.1017/S0003598X00092577">Fort 2003</a>), and the Bantu spread in Africa (<a href="https://doi.org/10.1371/journal.pone.0215573">Isern and Fort 2019</a>).</p>

<h3>Installation</h3>

<p>To install from CRAN:</p>

```R
install.packages("spDates")
```

<p>To install from the github repository:</p>

```R
devtools::install_github("jgregoriods/spDates")
```

<h3>Examples</h3>

<p>The package includes data sets with radiocarbon dates of Neolithic sites and potential centers of expansion modified from Pinhasi et al. (<a href="https://doi.org/10.1371/journal.pbio.0030410">2005</a>). The radiocarbon dates have already been filtered to retain only the earliest date per site - since including the more recent dates would affect the results of the regression (we are interested in the time of first arrival of the Neolithic). Let us load the data sets and perform a first regression of the dates versus distances from Jericho - a site that is commonly used as a hypothetical center of origin:</p>

```R
library(spDates)
data(neof)
data(centers)
jericho <- centers[centers$Site=="Jericho",]
model <- modelDates(neof, "C14Age", jericho, method="ols")
plot(model)
```

<img src="https://github.com/jgregoriods/spDates/blob/master/man/figures/model.jpeg" width="300">

<p>Normally, regression is performed on dates versus distances, given the assumption that most of the error will be concentrated on the former (<a href="https://doi.org/10.1371/journal.pbio.0030410">Pinhasi et al. 2005</a>). Nevertheless, distances may also be uncertain, with great-circle distances being only an approximation to the actual path travelled to the site. To account for that, regression on distances versus dates can also be run. In the plot above, the solid line corresponds to the dates-versus-distances regression, while the dashed line shows the distances-versus-dates regression.</p>
<p>To mitigate the uncertainty in radiocarbon dates, the robustness of the regression can be assessed by a bootstrapping procedure (<a href="https://doi.org/10.1017/S0003598X00061330">Gkiasta et al. 2003</a>). Here, the modelDates() function executes 999 regressions, each time sampling a single year from the calibrated age ranges. The lines of each regression are shown in the plot, providing an uncertanty envelope (red for dates-versus-distances, blue for distances-versus-dates). The black lines correspond to the average of each bootstrapping.</p>
<p>One can check the estimates for the expansion start date and speed:</p>

```R
summary(model)
                   Time-versus-distance  Distance-versus-time
Start date:          9782 +/- 11 cal BP   11115 +/- 21 cal BP
Speed of advance: 0.97 +/- 0.0073 km/yr 0.64 +/- 0.0063 km/yr
```

<p>Another method that has been used in time-space regressions is reduced major axis (RMA), which, unlike OLS, assumes a symmetrical distribution of error between both variables and has been shown to be robust to outliers (<a href="https://doi.org/10.1016/j.jas.2010.03.007">Steele 2010</a>; <a href="https://doi.org/10.1371/journal.pone.0087854">Russell et al. 2014</a>):</p>

```R
rmamodel <- modelDates(neof, "C14Age", jericho, method="rma")
plot(rmamodel)
```

<img src="https://github.com/jgregoriods/spDates/blob/master/man/figures/rmamodel.jpeg" width="300">

<p>So far, we have used all of the sites in the analysis. One can also apply a binning procedure to retain only the earliest site per spatial bins - defined by regular distance intervals from the hypothetical origin (<a href="https://doi.org/10.1073/pnas.0704215104">Hamilton and Buchanan 2007</a>; <a href="https://doi.org/10.1016/j.jas.2010.03.007">Steele 2010</a>). Let us apply spatial bins of 500 km (RMA is executed by default):</p>

```R
rmabins <- modelDates(neof, "C14Age", jericho, binWidth=500)
plot(rmabins)
```

<img src="https://github.com/jgregoriods/spDates/blob/master/man/figures/rmabins.jpeg" width="300">

<p>As mentioned above, some level of uncertainty has to be taken into account for the distances as well as for the dates. That is because the exact routes travelled are unknown, and, so far, all distances have been calculated from great circles. It is also possible to incorporate a cost surface in order to calculate least-cost paths. The package includes a cost surface where the coast is easier to travel, but sea and land above 1750 m are barriers:</p>

```R
data(cost)
rmacost <- modelDates(neof, "C14Age", jericho, binWidth=500, cost=cost)
plot(rmacost)
```

<img src="https://github.com/jgregoriods/spDates/blob/master/man/figures/rmacost.jpeg" width="300">

<p>Finally, one can iterate over many sites to test for hypothetical origins, selecting the one with the highest correlation coefficient as the most likely center of origin. Here, we will use 9 sites in the Near East that have been considered as potential Neolithic "cradles" by Pinhasi et al. (<a href="https://doi.org/10.1371/journal.pbio.0030410">2005</a>). We will use spatial bins of 500 km for all cases, but a sequence of widths can also be passed as an argument to test the effect of using different spatial bins:</p>

```R
iter <- iterateSites(neof, "C14Age", centers, "Site", binWidths=500)
iter$results

        r p bin  n              site
1: 0.9640 0 500  9           Aswad**
2: 0.9516 0 500  9          Cayönü**
3: 0.9496 0 500  8     Çatal Höyük**
4: 0.9488 0 500  8 Shillourokambos**
5: 0.9469 0 500  9         Jericho**
6: 0.9451 0 500  9     Qermez Dere**
7: 0.9429 0 500 10        Abu Madi**
8: 0.9316 0 500  9     Abu Hureyra**
9: 0.9188 0 500 11        Ali Kosh**

plot(iter$model)
```

<img src="https://github.com/jgregoriods/spDates/blob/master/man/figures/iter.jpeg" width="300">

<p>One can plot a map with the results of the iteration, showing an interpolated surface with the correlation coefficient of all sites tested as potential origins:</p>

```R
plot(iter$map)
```

<img src="https://github.com/jgregoriods/spDates/blob/master/man/figures/itermap.jpeg" width="400">

<h3>Important</h3>

<p>When preparing your own data, it is crucial that you include two columns named "cal" and "med", containing, respectively, calibrated dates in the form of CalDates objects (from the rcarbon package) and the median of each calibrated date (for display in the plots). These can be created in the following way:</p>

```R
library(rcarbon)
mydates <- read.csv("myfile.csv")
mydates$cal <- calibrate(mydates$C14Age, mydates$C14SD)
mydates$med <- medCal(mydates$cal)
```

<p>In this example we are reading a csv file but your dataset, of course, can already be stored as a SpatialPointsDataFrame.</p>

<h3>References</h3>
<p>Ammerman, A J, and L L Cavalli-Sforza. 1971. “Measuring the Rate of Spread of Early Farming in Europe.” Man 6 (4): 674–88.</p>
<p>Fort, Joaquim. 2003. <a href="https://doi.org/10.1017/S0959774304000046">“Population Expansion in the Western Pacific (Austronesia): A Wave of Advance Model.”</a> Antiquity 77 (297): 520–30.</p>
<p>Fort, Joaquim, Toni Pujol, and Luigi Luca Cavalli-Sforza. 2004. <a href="https://doi.org/10.1017/S0959774304000046">“Palaeolithic Populations and Waves of Advance.”</a> Cambridge Archaeological Journal 14 (1): 53–61.</p>
<p>Gkiasta, Marina, Thembi Russell, Stephen Shennan, and James Steele. 2003. <a href="https://doi.org/10.1017/S0003598X00061330">“Neolithic Transition in Europe: The Radiocarbon Record Revisited.”</a> Antiquity 77 (295): 45–62.</p>
<p>Hamilton, Marcus J, and Briggs Buchanan. 2007. <a href="https://doi.org/10.1073/pnas.0704215104">“Spatial Gradients in Clovis-Age Radiocarbon Dates across North America Suggest Rapid Colonization from the North.”</a> Proceedings of the National Academy of Sciences 104 (40): 15625–30.</p>
<p>Hamilton, Marcus J, and Briggs Buchanan. 2010. <a href="https://doi.org/10.1371/journal.pone.0012472">“Archaeological Support for the Three-Stage Expansion of Modern Humans across Northeastern Eurasia and into the Americas.”</a> PLoS One 5 (8): e12472.</p>
<p>Isern, Neus, and Joaquim Fort. 2019. <a href="https://doi.org/10.1371/journal.pone.0215573">“Assessing the Importance of Cultural Diffusion in the Bantu Spread into Southeastern Africa.”</a> PLOS ONE 14 (5): 1–18.</p>
<p>Pinhasi, Ron, Joaquim Fort, and Albert J Ammerman. 2005. <a href="https://doi.org/10.1371/journal.pbio.0030410">“Tracing the Origin and Spread of Agriculture in Europe.”</a> PLoS Biology 3 (12): e410.</p>
<p>Russell, Thembi, Fabio Silva, and James Steele. 2014. <a href="https://doi.org/10.1371/journal.pone.0087854">“Modelling the Spread of Farming in the Bantu-Speaking Regions of Africa: An Archaeology-Based Phylogeography.”</a> PLoS One 9 (1): e87854.</p>
<p>Steele, James. 2010. <a href="https://doi.org/10.1016/j.jas.2010.03.007">“Radiocarbon Dates as Data: Quantitative Strategies for Estimating Colonization Front Speeds and Event Densities.”</a> Journal of Archaeological Science 37 (8): 2017–30.</p>
