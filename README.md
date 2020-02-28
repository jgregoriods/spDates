<h1>spDates</h1>
<h2>Analysis of Spatial Gradients in Radiocarbon Dates</h2>

<p>This package offers tools to perform time-space regressions, an approach often used by archaeologists examining the expansion of cultural phenomena. In essence, one plots the radiocarbon ages of archaeological sites against their distances from a hypothetical origin. If a cultural advance has indeed taken place, the expectation is that radiocarbon dates will be more recent as one moves away from the center of origin. If a significant correlation is found, the intercept of the regression can be used as an estimate of the start date for the dispersal, while the regression slope provides an estimate of the speed of advance. Most applications have been focused on the Neolithic expansion from the Near East to Europe (Ammerman and Cavalli-Sforza 1971; Gkiasta et al. 2003; Pinhasi et al. 2005), but other case studies include the Paleolithic recolonization of Northern Europe (Fort et al. 2004), the Clovis expansion in North America (Hamilton and Buchanan 2007), the human colonization of the Americas from Asia (Hamilton and Buchanan 2010), the Lapita spread in Austronesia (Fort 2002), and the Bantu spread in Africa (Isern and Fort 2019).</p>

<h3>Installation</h3>

<p>To install from the github repository:</p>

<pre><code>devtools::install_github("jgregoriods/spDates")</pre></code>

<h3>Examples</h3>

<p>The package includes data sets with radiocarbon dates of Neolithic sites and potential centers of expansion modified from Pinhasi et al. (2005). The radiocarbon dates have already been filtered to retain only the earliest date per site - since including the more recent dates would affect the results of the regression (we are interested in the time of first arrival of the Neolithic). Let us load the data sets and perform a first regression of the dates versus distances from Jericho - a site that is commonly used as a hypothetical center of origin:</p>

<pre><code>library(spDates)
data("neof")
data("centers")
jericho <- centers[centers$Site=="Jericho",]
model <- modelDates(neof, "C14Age", jericho, method="ols")
plot(model)</pre></code>

<img src="https://github.com/jgregoriods/spDates/blob/master/model.jpeg" width="300">

<p>Normally, regression is performed on dates versus distances, given the assumption that most of the error will be concentrated on the former (Pinhasi et al. 2005). Nevertheless, distances may also be uncertain, with great-circle distances being only an approximation to the actual path travelled to the site. To account for that, regression on distances versus dates can also be run. In the plot above, the solid line corresponds to the dates-versus-distances regression, while the dashed line shows the distances-versus-dates regression.</p>
<p>To mitigate the uncertainty in radiocarbon dates, the robustness of the regression can be assessed by bootstrapping (Gkiasta et al. 2003). Here, the modelDates() function executes 999 regressions, each time sampling from the 
