# Class07
Qihao Liu (PID:U08901197)

Today we will explore some fundamental machine learning methods
including clustering and dimensionality reduction.

## K-means clustering

To see how this works, Let’s first make up some data to cluster where we
konw what the answer should be, we can use the `rnorm()` function to
help here.

- this gives us n number of data drawn from a normal distribution with
  pre-determined mean and SD.

``` r
rnorm(50)
```

     [1]  0.77300614  0.06353409  2.13129195 -0.45315588  0.86126340 -0.85725897
     [7] -0.27885652 -1.62571691  0.14258987  0.91244834 -1.63751177  0.72265343
    [13]  0.65249173  0.31576313  1.50201575 -0.07916741  0.63396716 -0.73936955
    [19] -0.16199600  0.67704077  0.42585021  2.95327825  1.10026512  0.07654442
    [25]  0.81256381  0.12183260 -0.31715268 -0.71432964  1.25712006  1.90861854
    [31]  0.17567563 -1.17235979 -1.41417270 -1.94515633 -1.10668505  0.21749153
    [37] -0.12347052 -0.82957970 -0.54912145  1.18259119 -0.03146326  0.70603290
    [43] -0.52224864  1.94809784 -0.83733106  0.07056984 -0.22061211  0.26798813
    [49] -0.86767334 -0.20396957

``` r
hist(rnorm(50,15,2))
```

![](class07_files/figure-commonmark/unnamed-chunk-1-1.png)

Now let’s make up a data in which there’s clearly two distinct groups

``` r
x <- c(rnorm(30, mean=-3), rnorm(30,mean=+3))
y <- rev(x) #reverse x

cbind(x,y) #column bind
```

                   x          y
     [1,] -2.6583209  3.7248024
     [2,] -2.4491343  3.1766113
     [3,] -2.8471253  3.9108827
     [4,] -2.0719146  4.0682433
     [5,] -4.3370205  3.4744595
     [6,] -3.2605061  3.6524321
     [7,] -3.8807720  3.6191795
     [8,] -2.0042183  2.8900215
     [9,] -2.7710801  3.5863785
    [10,] -2.7816296  2.3918025
    [11,] -2.8590628  3.4666928
    [12,] -2.0392687  4.0140325
    [13,] -3.6565463  2.3178300
    [14,] -3.3596218  2.6420861
    [15,] -2.9505363  3.8117434
    [16,] -3.6504794  4.3446319
    [17,] -1.7859288  4.2546278
    [18,] -3.0896978  3.2888149
    [19,] -4.6578339  2.2730711
    [20,] -1.3047579  1.2078571
    [21,] -2.4450843  1.9319427
    [22,] -3.2751904  0.9506581
    [23,] -3.7435933  3.5199516
    [24,] -1.9235416  3.4815341
    [25,] -1.7112184  3.2264439
    [26,] -2.8392728  3.2928864
    [27,] -2.7792043  4.0946143
    [28,] -2.4053621  3.6033361
    [29,] -3.0962577  2.6866265
    [30,] -2.0036935  1.6482107
    [31,]  1.6482107 -2.0036935
    [32,]  2.6866265 -3.0962577
    [33,]  3.6033361 -2.4053621
    [34,]  4.0946143 -2.7792043
    [35,]  3.2928864 -2.8392728
    [36,]  3.2264439 -1.7112184
    [37,]  3.4815341 -1.9235416
    [38,]  3.5199516 -3.7435933
    [39,]  0.9506581 -3.2751904
    [40,]  1.9319427 -2.4450843
    [41,]  1.2078571 -1.3047579
    [42,]  2.2730711 -4.6578339
    [43,]  3.2888149 -3.0896978
    [44,]  4.2546278 -1.7859288
    [45,]  4.3446319 -3.6504794
    [46,]  3.8117434 -2.9505363
    [47,]  2.6420861 -3.3596218
    [48,]  2.3178300 -3.6565463
    [49,]  4.0140325 -2.0392687
    [50,]  3.4666928 -2.8590628
    [51,]  2.3918025 -2.7816296
    [52,]  3.5863785 -2.7710801
    [53,]  2.8900215 -2.0042183
    [54,]  3.6191795 -3.8807720
    [55,]  3.6524321 -3.2605061
    [56,]  3.4744595 -4.3370205
    [57,]  4.0682433 -2.0719146
    [58,]  3.9108827 -2.8471253
    [59,]  3.1766113 -2.4491343
    [60,]  3.7248024 -2.6583209

``` r
z <- cbind(x,y)
plot(z)
```

![](class07_files/figure-commonmark/unnamed-chunk-2-1.png)

- Q: What does cbind do?

The function for K-means clustering in “base” R is called `kmeans()` (k
in kmeans means how many cluster you want)

``` r
k <- kmeans(z, centers = 2)
k
```

    K-means clustering with 2 clusters of sizes 30, 30

    Cluster means:
              x         y
    1  3.151747 -2.821262
    2 -2.821262  3.151747

    Clustering vector:
     [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1
    [39] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

    Within cluster sum of squares by cluster:
    [1] 41.05722 41.05722
     (between_SS / total_SS =  92.9 %)

    Available components:

    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

To get the results of the returned list object, we can use the dollar
`$` syntax, for example, to get the size of each cluster:

- Note the components listed fro K-means are essentially different types
  of results

> Q: how many points are in each cluster

``` r
k$size
```

    [1] 30 30

> Q: What ‘component’ gives - Membership/ Cluster Assignment - Cluster
> Center

``` r
k$cluster #Membership
```

     [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1
    [39] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

``` r
k$centers #Cluster Center
```

              x         y
    1  3.151747 -2.821262
    2 -2.821262  3.151747

> Q: Make a result figure of the data colored by cluster membership and
> show cluster centers

``` r
plot(z, col = c("red","darkgreen")) #this colors the dots in alternating colors
```

![](class07_files/figure-commonmark/unnamed-chunk-6-1.png)

``` r
plot(z, col = 2) #color by bumber
```

![](class07_files/figure-commonmark/unnamed-chunk-6-2.png)

``` r
plot(z, col = k$cluster, pch=20) # this basically tells R, use the number indicating each pt's membership to color them.
points(k$centers, col = "lightblue", cex=2, pch=15)
```

![](class07_files/figure-commonmark/unnamed-chunk-6-3.png)

K-means clustering is very popular as it is very fast and relatively
straight forward: it takes numeric data as input and returns the cluster
membership vector etc. But the “issue/feature” of it is that you need to
tell `kmeans()` how many cluster we want!

> Q: Run Kmeans() again and cluster into 4 groups, and plot the results
> So the result is not ideal…

``` r
k4 <- kmeans(z, centers=4)
k4$cluster
```

     [1] 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 2 3 1 1 1 1 1 3
    [39] 2 2 2 3 3 1 3 1 3 3 1 1 2 1 1 3 3 3 1 1 1 1

``` r
k4$centers
```

              x         y
    1  3.640190 -2.406346
    2  1.626094 -2.362071
    3  3.181908 -3.673233
    4 -2.821262  3.151747

``` r
plot(z, col = k4$cluster, pch=20)
points(k$centers, col = "lightblue", cex=2, pch=15)
```

![](class07_files/figure-commonmark/unnamed-chunk-7-1.png)

An alternative way is to make a **Scree Plot**, graphing total variation
within a cluster (measured as Total Within-Cluster Sum of
Squares)against the number of centers, and we can take the “elbow”
center (the centers n that has the greatest decrease in total variation
within a cluster compare to n-1) to be the optimal one.

- Note scree is the place underneath a cliff that the rubbles collects

``` r
k1 <- kmeans(z, centers=1)
k2 <- kmeans(z, centers=2)
k3 <- kmeans(z, centers=3)
k4 <- kmeans(z, centers=4)
k4$tot.withinss
```

    [1] 47.88374

``` r
k5 <- kmeans(z, centers=5)
k6 <- kmeans(z, centers=6)
k7 <- kmeans(z, centers=7)
x <- c(1,2,3,4,5,6,7)
y <- c(k1$tot.withinss, k2$tot.withinss, k3$tot.withinss, k4$tot.withinss, k5$tot.withinss, k6$tot.withinss, k7$tot.withinss)
plot(x,y, type="l", xlab = "Centers", ylab = "Total Within SS")
text(x[2], y[2], labels = "centers = 2", pos = 3, col = "darkred")
```

![](class07_files/figure-commonmark/unnamed-chunk-8-1.png)

or we can use `for()` loop

``` r
n <- NULL
for(i in 1:5) {
  n <- c(n,kmeans(x, centers = i)$tot.withinss)
  plot(n, typ="b")
}
```

![](class07_files/figure-commonmark/unnamed-chunk-9-1.png)

![](class07_files/figure-commonmark/unnamed-chunk-9-2.png)

![](class07_files/figure-commonmark/unnamed-chunk-9-3.png)

![](class07_files/figure-commonmark/unnamed-chunk-9-4.png)

![](class07_files/figure-commonmark/unnamed-chunk-9-5.png)

## Hierarchical Clustering

this meain “base” R function for Hierarchical Clustering is called
`hclust()`

- the only required input of `hclust()` is d, a disimilarity matrix
  generated by the `dist()` function

``` r
hclust()
```

    Error in hclust(): argument "d" is missing, with no default

- `dist()` function calculate a distance matrix for our data, its a
  analysis of all paralle comparision of every point in data (it gives
  us a triangle because it is symmetrical)

``` r
d <- dist(z)
hc <- hclust(d)
hc
```


    Call:
    hclust(d = d)

    Cluster method   : complete 
    Distance         : euclidean 
    Number of objects: 60 

``` r
plot(hc)
abline(h=10, col = "red")
```

![](class07_files/figure-commonmark/unnamed-chunk-11-1.png)

``` r
cutree(hc,h=10)
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

if we look at the dendrogram, we can clearly see that there are two
large clusters, and higher numbers are clustered on one side of the
plot, and lower numbers are clustered on the other side of the plot
(this is related to how we build the data, we get 30 numbers with one
mean, and 30 other numbers with a different mean)

Horizontal line contains the important information. We can see that if
we cut the dendrogram at the dark red line with the `cutree ()`
function, it would yield a membership output

To get our cluster membership output (i.e. our main clustering results),
we can cut the tree at a given height or at a height that yields a diven
“k”

``` r
grps <- cutree(hc,h=10)
#OR
cutree(hc,k=2)
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

> Q: Plot the data with our hclut results coloring

``` r
plot(z, col= grps)
```

![](class07_files/figure-commonmark/unnamed-chunk-13-1.png)

## Dimensionality Reduction: Principle Component Analysis(PCA)

The major goal of PCA is dimensionality reduction and lost as little
information/essence of data as possible

2D-Example: - The first principle component (PC) follows a best fit
through the data points. - The second principle component picks up the
wiggling of the data along the first PC - Combined, these captures most
of the spread of the data and can act as new axis - The data should have
max variation on PC1, PC2 captures the rest of the variants (Generally,
for PC1,2,3, the variation captured decreases as the numbering goes up)

# Principal Component Analysis (PCA)

## PCA of UK Food Data

Import food data from an online CSV file.

``` r
url <- "https://tinyurl.com/UK-foods"
x <- read.csv(url)
head(x)
```

                   X England Wales Scotland N.Ireland
    1         Cheese     105   103      103        66
    2  Carcass_meat      245   227      242       267
    3    Other_meat      685   803      750       586
    4           Fish     147   160      122        93
    5 Fats_and_oils      193   235      184       209
    6         Sugars     156   175      147       139

> Q2. Which approach to solving the ‘row-names problem’ mentioned above
> do you prefer and why? Is one approach more robust than another under
> certain circumstances?

the for the first column, the name is X, and the row names are not food
but numbers To fix this: We can try….

``` r
rownames(x) <- x[,1]
x <- x[,-1]
x
```

                        England Wales Scotland N.Ireland
    Cheese                  105   103      103        66
    Carcass_meat            245   227      242       267
    Other_meat              685   803      750       586
    Fish                    147   160      122        93
    Fats_and_oils           193   235      184       209
    Sugars                  156   175      147       139
    Fresh_potatoes          720   874      566      1033
    Fresh_Veg               253   265      171       143
    Other_Veg               488   570      418       355
    Processed_potatoes      198   203      220       187
    Processed_Veg           360   365      337       334
    Fresh_fruit            1102  1137      957       674
    Cereals                1472  1582     1462      1494
    Beverages                57    73       53        47
    Soft_drinks            1374  1256     1572      1506
    Alcoholic_drinks        375   475      458       135
    Confectionery            54    64       62        41

This is no good because the code is overwriting itself, so each time we
run it, we lose some data.

Instead, we do

``` r
x <- read.csv(url,row.names=1)
x
```

                        England Wales Scotland N.Ireland
    Cheese                  105   103      103        66
    Carcass_meat            245   227      242       267
    Other_meat              685   803      750       586
    Fish                    147   160      122        93
    Fats_and_oils           193   235      184       209
    Sugars                  156   175      147       139
    Fresh_potatoes          720   874      566      1033
    Fresh_Veg               253   265      171       143
    Other_Veg               488   570      418       355
    Processed_potatoes      198   203      220       187
    Processed_Veg           360   365      337       334
    Fresh_fruit            1102  1137      957       674
    Cereals                1472  1582     1462      1494
    Beverages                57    73       53        47
    Soft_drinks            1374  1256     1572      1506
    Alcoholic_drinks        375   475      458       135
    Confectionery            54    64       62        41

Now we can try some base figures:

``` r
barplot(as.matrix(x), beside=T, col=rainbow(nrow(x)))
```

![](class07_files/figure-commonmark/unnamed-chunk-17-1.png)

> Q3: Changing what optional argument in the above barplot() function
> results in the following plot?

``` r
barplot(as.matrix(x), beside=F, col=rainbow(nrow(x)))
```

![](class07_files/figure-commonmark/unnamed-chunk-18-1.png)

These plot are usless

> Q5: We can use the pairs() function to generate all pairwise plots for
> our countries. Can you make sense of the following code and resulting
> figure? What does it mean if a given point lies on the diagonal for a
> given plot?

There is one plot that can be useful FOR SAMLL DATASET

``` r
pairs(x, col=rainbow(10), pch=16)
```

![](class07_files/figure-commonmark/unnamed-chunk-19-1.png)

To read this plot, if we take the first row where it starts with
“England”, we can see England vs Wales, England vs Scotland, England vs
N. Ireland

- If the country eats exactly the same, it will be an diagnal line
- deviation from the diagnal highligths the differences in eating
- from this we can see that N. Ireland is very different \>Main point:
  it can be difficult to spot major trends and patterns even in
  relatively small multivariate dataset (here we only have 17 dimension,
  typically it could be 1000s)

## PCA to the rescue

Now let’s try PCA on this data to see if it is actually useful. The main
functionin “base” R for PCA is called `prcomp()`

- First we transpose the data so food is in the columns. This allow us
  to do PCA on the food with `t()` transpose function.

``` r
t(x)
```

              Cheese Carcass_meat  Other_meat  Fish Fats_and_oils  Sugars
    England      105           245         685  147            193    156
    Wales        103           227         803  160            235    175
    Scotland     103           242         750  122            184    147
    N.Ireland     66           267         586   93            209    139
              Fresh_potatoes  Fresh_Veg  Other_Veg  Processed_potatoes 
    England               720        253        488                 198
    Wales                 874        265        570                 203
    Scotland              566        171        418                 220
    N.Ireland            1033        143        355                 187
              Processed_Veg  Fresh_fruit  Cereals  Beverages Soft_drinks 
    England              360         1102     1472        57         1374
    Wales                365         1137     1582        73         1256
    Scotland             337          957     1462        53         1572
    N.Ireland            334          674     1494        47         1506
              Alcoholic_drinks  Confectionery 
    England                 375             54
    Wales                   475             64
    Scotland                458             62
    N.Ireland               135             41

``` r
PCA <- prcomp(t(x))
PCA
```

    Standard deviations (1, .., p=4):
    [1] 3.241502e+02 2.127478e+02 7.387622e+01 3.175833e-14

    Rotation (n x k) = (17 x 4):
                                 PC1          PC2         PC3          PC4
    Cheese              -0.056955380  0.016012850  0.02394295 -0.694538519
    Carcass_meat         0.047927628  0.013915823  0.06367111  0.489884628
    Other_meat          -0.258916658 -0.015331138 -0.55384854  0.279023718
    Fish                -0.084414983 -0.050754947  0.03906481 -0.008483145
    Fats_and_oils       -0.005193623 -0.095388656 -0.12522257  0.076097502
    Sugars              -0.037620983 -0.043021699 -0.03605745  0.034101334
    Fresh_potatoes       0.401402060 -0.715017078 -0.20668248 -0.090972715
    Fresh_Veg           -0.151849942 -0.144900268  0.21382237 -0.039901917
    Other_Veg           -0.243593729 -0.225450923 -0.05332841  0.016719075
    Processed_potatoes  -0.026886233  0.042850761 -0.07364902  0.030125166
    Processed_Veg       -0.036488269 -0.045451802  0.05289191 -0.013969507
    Fresh_fruit         -0.632640898 -0.177740743  0.40012865  0.184072217
    Cereals             -0.047702858 -0.212599678 -0.35884921  0.191926714
    Beverages           -0.026187756 -0.030560542 -0.04135860  0.004831876
    Soft_drinks          0.232244140  0.555124311 -0.16942648  0.103508492
    Alcoholic_drinks    -0.463968168  0.113536523 -0.49858320 -0.316290619
    Confectionery       -0.029650201  0.005949921 -0.05232164  0.001847469

``` r
summary(PCA)
```

    Importance of components:
                                PC1      PC2      PC3       PC4
    Standard deviation     324.1502 212.7478 73.87622 3.176e-14
    Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00

Note this summary tells us that PC1 captures 67.44% of the variants, PC2
29.05%, and PC3 3.5%. Together, PC1 and 2 captures \>95% of the variance

Now let’s try to plot the results in “base” R and ggplot

- First figure out what we are plotting, we are plotting X, which
  contains the PCs, and since PC1 and PC2 capture most of the variation,
  we are only gonna plot them two

``` r
PCA$x
```

                     PC1         PC2        PC3           PC4
    England   -144.99315   -2.532999 105.768945 -4.894696e-14
    Wales     -240.52915 -224.646925 -56.475555  5.700024e-13
    Scotland   -91.86934  286.081786 -44.415495 -7.460785e-13
    N.Ireland  477.39164  -58.901862  -4.877895  2.321303e-13

``` r
cols <- c("orange","darkgreen","darkblue","darkred")
plot(PCA$x[,1], PCA$x[,2], col=cols,pch=16, xlab = "PC1", ylab = "PC2")
```

![](class07_files/figure-commonmark/unnamed-chunk-21-1.png)

If we are gonna use ggplot \>Q7. Complete the code below to generate a
plot of PC1 vs PC2. The second line adds text labels over the data
points.

``` r
library(ggplot2)
ggplot(PCA$x)+
  aes(PC1,PC2)+
  geom_point(col=cols)
```

![](class07_files/figure-commonmark/unnamed-chunk-22-1.png)

From this graph, we can clearly see that N. Ireland is distinct from the
rest of the countrys from the survey.

But what food contributes to these PCs? We can figure this out, with
`PCA$rotation`, we can see the contribution of each food to the
different PCs, now if we plot it we can visualize what contributes to
the PCs

``` r
PCA$rotation
```

                                 PC1          PC2         PC3          PC4
    Cheese              -0.056955380  0.016012850  0.02394295 -0.694538519
    Carcass_meat         0.047927628  0.013915823  0.06367111  0.489884628
    Other_meat          -0.258916658 -0.015331138 -0.55384854  0.279023718
    Fish                -0.084414983 -0.050754947  0.03906481 -0.008483145
    Fats_and_oils       -0.005193623 -0.095388656 -0.12522257  0.076097502
    Sugars              -0.037620983 -0.043021699 -0.03605745  0.034101334
    Fresh_potatoes       0.401402060 -0.715017078 -0.20668248 -0.090972715
    Fresh_Veg           -0.151849942 -0.144900268  0.21382237 -0.039901917
    Other_Veg           -0.243593729 -0.225450923 -0.05332841  0.016719075
    Processed_potatoes  -0.026886233  0.042850761 -0.07364902  0.030125166
    Processed_Veg       -0.036488269 -0.045451802  0.05289191 -0.013969507
    Fresh_fruit         -0.632640898 -0.177740743  0.40012865  0.184072217
    Cereals             -0.047702858 -0.212599678 -0.35884921  0.191926714
    Beverages           -0.026187756 -0.030560542 -0.04135860  0.004831876
    Soft_drinks          0.232244140  0.555124311 -0.16942648  0.103508492
    Alcoholic_drinks    -0.463968168  0.113536523 -0.49858320 -0.316290619
    Confectionery       -0.029650201  0.005949921 -0.05232164  0.001847469

``` r
ggplot(PCA$rotation)+
  aes(PC1,rownames(PCA$rotation))+
  geom_col()
```

![](class07_files/figure-commonmark/unnamed-chunk-23-1.png)

This graph tells us that Ireland consume more soft drink and fresh
potatoes

PCA looks super useful, and we will come back to describe this further
next class
