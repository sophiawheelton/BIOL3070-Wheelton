Warm-up mini-Report: Mosquito Blood Hosts in Salt Lake City, Utah
================
Sophia Wheelton
2025-10-27

- [Scatterplot with line(s) of best
  fit](#scatterplot-with-lines-of-best-fit)
- [ABSTRACT](#abstract)
- [BACKGROUND](#background)
- [STUDY QUESTION and HYPOTHESIS](#study-question-and-hypothesis)
  - [Questions](#questions)
  - [Hypothesis](#hypothesis)
  - [Prediction](#prediction)
- [METHODS](#methods)
  - [First Analysis and Plot](#first-analysis-and-plot)
  - [Second Analysis](#second-analysis)
- [DISCUSSION](#discussion)
  - [Interpretation of First
    Analysis](#interpretation-of-first-analysis)
  - [Interpretation of Second
    Analysis](#interpretation-of-second-analysis)
- [CONCLUSION](#conclusion)
- [REFERENCES](#references)

``` r
# load data
life_expectancy <- readr::read_csv('https://raw.githubusercontent.com/rfordatascience/tidytuesday/main/data/2023/2023-12-05/life_expectancy.csv')
```

    ## `curl` package not installed, falling back to using `url()`
    ## Rows: 20755 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): Entity, Code
    ## dbl (2): Year, LifeExpectancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# Linear regression years before 1883 (Koch's germ theory)
lm_b1883 <- lm(LifeExpectancy ~ Year, data = life_expectancy[life_expectancy$Year < 1883,])
summary(lm_b1883)
```

    ## 
    ## Call:
    ## lm(formula = LifeExpectancy ~ Year, data = life_expectancy[life_expectancy$Year < 
    ##     1883, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -23.204  -2.009   0.527   2.682  12.073 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -7.294260   5.899927  -1.236    0.217    
    ## Year         0.025738   0.003214   8.008 5.33e-15 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 4.979 on 654 degrees of freedom
    ## Multiple R-squared:  0.08931,    Adjusted R-squared:  0.08791 
    ## F-statistic: 64.13 on 1 and 654 DF,  p-value: 5.329e-15

``` r
# Linear regression years after 1883 (Koch's germ theory)
lm_a1883 <- lm(LifeExpectancy ~ Year, data = life_expectancy[life_expectancy$Year > 1883,])
summary(lm_a1883)
```

    ## 
    ## Call:
    ## lm(formula = LifeExpectancy ~ Year, data = life_expectancy[life_expectancy$Year > 
    ##     1883, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -52.007  -6.987   2.335   7.901  18.287 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -4.781e+02  5.201e+00  -91.93   <2e-16 ***
    ## Year         2.729e-01  2.626e-03  103.92   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 10.05 on 20083 degrees of freedom
    ## Multiple R-squared:  0.3497, Adjusted R-squared:  0.3497 
    ## F-statistic: 1.08e+04 on 1 and 20083 DF,  p-value: < 2.2e-16

# Scatterplot with line(s) of best fit

``` r
plot(life_expectancy$Year,life_expectancy$LifeExpectancy, pch = 19, col = "grey", xlab = "Year of Birth", ylab = "Life Expectancy")
abline(lm_b1883, lt = 3)
abline(lm_a1883, lt = 1)
```

![](Warm-up-mosquitoes-Wheelton_files/figure-gfm/scatterplot-1.png)<!-- -->

# ABSTRACT

West Nile Virus is potentially deadly upon infection. Mosquitoes get
infected upon consuming a blood meal from certain birds and then
transmitting it to their next blood meal host. These hosts could be
other birds, or they could be humans. We sought to determine if there
was a relationship between house finches and West Nile Virus in Salt
Lake City. Specifically, we determined if house finches were important
amplifying hosts of West Nile Virus. Analysis of house finch populations
in areas (hot spots) with and without West Nile Virus visually suggested
a relationship. This relationship was further statistically confirmed
with a Generalized Linear Model (GLM) statistical test. The test showed
that there was a significant difference between house finch blood meal
counts in locations with and without West Nile Virus, implying that the
relationship between house finches and West Nile Virus is predictive,
and that the house finches are acting as amplifying hosts. As the blood
meal counts of house finches in locations with West Nile virus is much
larger than the blood meal counts of other birds in the same location,
we can infer that house finches are especially important for
amplification.

# BACKGROUND

West Nile Virus is typically spread from certain bird species to
mosquitoes and back again. The mosquito collects a blood meal from an
infected bird and then infects their next blood meal. While that next
blood meal is typically a bird, the mosquito can also affect people and
horses accidentally. In people, West Nile Virus is generally manageable.
However, in some cases, it results in death. Research into the
transmittability of West Nile Virus is important to help manage West
Nile outbreaks and to educate the public so they can help protect
themselves.

A lot of th research related to West Nile Virus analyzes blood meal
counts of different hosts of West Nile virus. To do this, DNA must be
extracted and identified. PCR amplifies the extracted DNA through
several cycles of heating and cooling in a Thermocycler. The DNA is
denatured so that a primer can bind to it. DNA polymerase then
replicates the strand. This denaturation and replication happens several
times until there is enough DNA to be sequenced (“Polymerase chain
reaction (PCR) fact sheet”, 2020). Sequencing can then be completed
using numerous methods. MinION sequencing uses a nanopore sequencer to
sequence the DNA, which is compiled using base-calling and aligned using
read counting and a reference gene sequences. The cytochrome oxidase I
(COI) gene is used to compare different species and individuals within a
species. The COI can be identified using the GenBank database.

By analyzing the DNA sequences of the blood meals of mosquitoes, the
hosts of the West Nile Virus-infected mosquitoes can be identified.
Using the viremia duration (Kumar et al., 2003) bar plot, we can
visually determine the importance of house finches in WNV transmission.
This provides the foundation for our hypothesis and prediction that
house finches are important to the amplification of West Nile Virus, and
that locations with more house finches will be predictive of West Nile
Virus hot spots.

The relationship between certain bird species and West Nile Virus
transmission is crucial for prevention efforts. By knowing how the virus
is spread, we can pay more attention to certain bird population
locations and focus our preventative measures in those areas.

``` r
# Manually transcribe duration (mean, lo, hi) from the last table column
duration <- data.frame(
  Bird = c("Canada Goose","Mallard", 
           "American Kestrel","Northern Bobwhite",
           "Japanese Quail","Ring-necked Pheasant",
           "American Coot","Killdeer",
           "Ring-billed Gull","Mourning Dove",
           "Rock Dove","Monk Parakeet",
           "Budgerigar","Great Horned Owl",
           "Northern Flicker","Blue Jay",
           "Black-billed Magpie","American Crow",
           "Fish Crow","American Robin",
           "European Starling","Red-winged Blackbird",
           "Common Grackle","House Finch","House Sparrow"),
  mean = c(4.0,4.0,4.5,4.0,1.3,3.7,4.0,4.5,5.5,3.7,3.2,2.7,1.7,6.0,4.0,
           4.0,5.0,3.8,5.0,4.5,3.2,3.0,3.3,6.0,4.5),
  lo   = c(3,4,4,3,0,3,4,4,4,3,3,1,0,6,3,
           3,5,3,4,4,3,3,3,5,2),
  hi   = c(5,4,5,5,4,4,4,5,7,4,4,4,4,6,5,
           5,5,5,7,5,4,3,4,7,6)
)

# Choose some colors
cols <- c(rainbow(30)[c(10:29,1:5)])  # rainbow colors

# horizontal barplot
par(mar=c(5,12,2,2))  # wider left margin for names
bp <- barplot(duration$mean, horiz=TRUE, names.arg=duration$Bird,
              las=1, col=cols, xlab="Days of detectable viremia", xlim=c(0,7))

# add error bars
arrows(duration$lo, bp, duration$hi, bp,
       angle=90, code=3, length=0.05, col="black", xpd=TRUE)
```

<img src="Warm-up-mosquitoes-Wheelton_files/figure-gfm/viremia-1.png" style="display: block; margin: auto auto auto 0;" />

# STUDY QUESTION and HYPOTHESIS

## Questions

What bird species is acting as WNV amplifying host in Salt Lake City?

## Hypothesis

House finches are acting as important amplifying hosts of WNV in Salt
Lake City.

## Prediction

If house finches are acting as important amplifying hosts, we predict
that trapping locations where mosquitoes feed on house finches will also
have higher rates of confirmed WNV in tested mosquito pools.

# METHODS

Mosquitoes were collected from several locations. We extracted blood
meal DNA from a sample of mosquitoes, both with West Nile Virus and
without. We then pipetted the collected DNA and controls into the wells
of the PCR plate. A primer and DNA polymerase was added and the DNA was
replicated in a Thermocycler. Using MinION sequencing, the mosquito COI
sequence was compiled. We identified the DNA sequences using the NCBI
BLASTn tool in the GenBank database. Finally, we plotted the data and
performed a Generalized Linear Model (GLM) statistical test to analyze
the results.

## First Analysis and Plot

Horizontal plots:

``` r
#This code creates a plot showing the frequency of different bird species at locations with and without West Nile Virus.
## import counts_matrix: data.frame with column 'loc_positives' (0/1) and host columns 'host_*'
counts_matrix <- read.csv("./bloodmeal_plusWNV_for_BIOL3070.csv")

## 1) Identify host columns
host_cols <- grep("^host_", names(counts_matrix), value = TRUE)

if (length(host_cols) == 0) {
  stop("No columns matching '^host_' were found in counts_matrix.")
}

## 2) Ensure loc_positives is present and has both levels 0 and 1 where possible
counts_matrix$loc_positives <- factor(counts_matrix$loc_positives, levels = c(0, 1))

## 3) Aggregate host counts by loc_positives
agg <- stats::aggregate(
  counts_matrix[, host_cols, drop = FALSE],
  by = list(loc_positives = counts_matrix$loc_positives),
  FUN = function(x) sum(as.numeric(x), na.rm = TRUE)
)

## make sure both rows exist; if one is missing, add a zero row
need_levels <- setdiff(levels(counts_matrix$loc_positives), as.character(agg$loc_positives))
if (length(need_levels)) {
  zero_row <- as.list(rep(0, length(host_cols)))
  names(zero_row) <- host_cols
  for (lv in need_levels) {
    agg <- rbind(agg, c(lv, zero_row))
  }
  ## restore proper type
  agg$loc_positives <- factor(agg$loc_positives, levels = c("0","1"))
  ## coerce numeric host cols (they may have become character after rbind)
  for (hc in host_cols) agg[[hc]] <- as.numeric(agg[[hc]])
  agg <- agg[order(agg$loc_positives), , drop = FALSE]
}

## 4) Decide species order (overall abundance, descending)
overall <- colSums(agg[, host_cols, drop = FALSE], na.rm = TRUE)
host_order <- names(sort(overall, decreasing = TRUE))
species_labels <- rev(sub("^host_", "", host_order))  # nicer labels

## 5) Build count vectors for each panel in the SAME order
counts0 <- rev(as.numeric(agg[agg$loc_positives == 0, host_order, drop = TRUE]))
counts1 <- rev(as.numeric(agg[agg$loc_positives == 1, host_order, drop = TRUE]))

## 6) Colors: reuse your existing 'cols' if it exists and is long enough; otherwise generate
if (exists("cols") && length(cols) >= length(host_order)) {
  species_colors <- setNames(cols[seq_along(host_order)], species_labels)
} else {
  species_colors <- setNames(rainbow(length(host_order) + 10)[seq_along(host_order)], species_labels)
}

## 7) Shared x-limit for comparability
xmax <- max(c(counts0, counts1), na.rm = TRUE)
xmax <- if (is.finite(xmax)) xmax else 1
xlim_use <- c(0, xmax * 1.08)

## 8) Plot: two horizontal barplots with identical order and colors
op <- par(mfrow = c(1, 2),
          mar = c(4, 12, 3, 2),  # big left margin for species names
          xaxs = "i")           # a bit tighter axis padding

## Panel A: No WNV detected (loc_positives = 0)
barplot(height = counts0,
        names.arg = species_labels, 
        cex.names = .5,
        cex.axis = .5,
        col = rev(unname(species_colors[species_labels])),
        horiz = TRUE,
        las = 1,
        xlab = "Bloodmeal counts",
        main = "Locations WNV (-)",
        xlim = xlim_use)

## Panel B: WNV detected (loc_positives = 1)
barplot(height = counts1,
        names.arg = species_labels, 
        cex.names = .5,
        cex.axis = .5,
        col = rev(unname(species_colors[species_labels])),
        horiz = TRUE,
        las = 1,
        xlab = "Bloodmeal counts",
        main = "Locations WNV (+)",
        xlim = xlim_use)
```

<img src="Warm-up-mosquitoes-Wheelton_files/figure-gfm/horiz-plot-1.png" style="display: block; margin: auto auto auto 0;" />

``` r
par(op)

## Keep the colors mapping for reuse elsewhere
host_species_colors <- species_colors
```

## Second Analysis

``` r
#This code is for a Generalized Linear Model statistical test. The test determines whether the presence or number of house finch blood meals predicts whether a site had WNV-positive pools (binary) or a higher WNV positivity rate (numeric). This test is a formal evaluation of the relationship suggested in the plot from the first analysis. 

# second-analysis-or-plot, glm with house finch alone against binary +/_
glm1 <- glm(loc_positives ~ host_House_finch,
            data = counts_matrix,
            family = binomial)
summary(glm1)
```

    ## 
    ## Call:
    ## glm(formula = loc_positives ~ host_House_finch, family = binomial, 
    ##     data = counts_matrix)
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)       -0.1709     0.1053  -1.622   0.1047  
    ## host_House_finch   0.3468     0.1586   2.187   0.0287 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 546.67  on 394  degrees of freedom
    ## Residual deviance: 539.69  on 393  degrees of freedom
    ## AIC: 543.69
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
#glm with house-finch alone against positivity rate
glm2 <- glm(loc_rate ~ host_House_finch,
            data = counts_matrix)
summary(glm2)
```

    ## 
    ## Call:
    ## glm(formula = loc_rate ~ host_House_finch, data = counts_matrix)
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)      0.054861   0.006755   8.122 6.07e-15 ***
    ## host_House_finch 0.027479   0.006662   4.125 4.54e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for gaussian family taken to be 0.01689032)
    ## 
    ##     Null deviance: 6.8915  on 392  degrees of freedom
    ## Residual deviance: 6.6041  on 391  degrees of freedom
    ##   (2 observations deleted due to missingness)
    ## AIC: -484.56
    ## 
    ## Number of Fisher Scoring iterations: 2

# DISCUSSION

Overall, we found that there was a significant relationship between
house finch populations and locations for West Nile Virus. This implies
that house finches are important to the spread of West Nile Virus.
However, there might be other confounding variables that might result in
an increased population of both. For example, in suburban areas, there
would be a higher population of house finches and a higher concentration
of West Nile Virus in general. This would weaken the statistical
connection between the two. Additional research should be done to
analyze potential confounding variables. A possible limitation is that
the data is only about Salt Lake City. The relationship should not be
applied outside of the Salt Lake City Valley without further research.

## Interpretation of First Analysis

Visually, there is a significantly greater number of house finches
(blood meal counts) in locations with West Nile Virus than there are in
locations without. This plot shows that there is a relationship between
house finches and West Nile Virus. While this could indicated that house
finches are important in the spread of West Nile Virus, without
performing a statistical test, we cannot claim that the relationship is
not based on chance. The plot also does not tell us whether house
finches can be predictive of West Nile Virus.

## Interpretation of Second Analysis

This GLM statistical test shows us that the relationship between house
finches (blood meal counts) and West Nile Virus is significant.
Specifically, the p-values are both significant at 6.07e-15 and
4.54e-05. This means that more house finches will be associated with and
predictive of more West Nile Virus, and that house finches are an
important amplifying host for West Nile Virus in Salt Lake City.

# CONCLUSION

We found a statistically significant link between house finches
population and hot spots for West Nile Virus when compared to house
finch populations without West Nile Virus. This relationship is
predictive. This shows that finches are acting as important amplifying
hosts of West Nile Virus in Salt Lake City. Additional research should
be done to identify and account for other confounding variables that
might be present, or if the relationship is applicable outside the Salt
Lake City Valley. Knowing what birds relate to increased spreadability
of West Nile Virus can help us with certain preventative measures.
Increased public awareness of the presence of certain bird species
(i.e. house finches) in an area might help make people more aware of the
risk of contracting West Nile Virus. This would further encourage them
to protect themselves against mosquitoes during high-risk seasons for
West Nile Virus.

# REFERENCES

1.  Komar N, Langevin S, Hinten S, Nemeth N, Edwards E, Hettler D, Davis
    B, Bowen R, Bunning M. Experimental infection of North American
    birds with the New York 1999 strain of West Nile virus. Emerg Infect
    Dis. 2003 Mar;9(3):311-22. <https://doi.org/10.3201/eid0903.020628>

2.  ChatGPT. OpenAI, version Jan 2025. Used as a reference for functions
    such as plot() and to correct syntax errors. Accessed 2025-10-27.

3.  Polymerase chain reaction (PCR) fact sheet. Genome.gov. (2020).
    <https://www.genome.gov/about-genomics/fact-sheets/Polymerase-Chain-Reaction-Fact-Sheet>
