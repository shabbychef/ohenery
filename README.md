

# ohenery

[![Build Status](https://travis-ci.org/shabbychef/ohenery.png)](https://travis-ci.org/shabbychef/ohenery)
[![codecov.io](http://codecov.io/github/shabbychef/ohenery/coverage.svg?branch=master)](http://codecov.io/github/shabbychef/ohenery?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/ohenery)](https://cran.r-project.org/package=ohenery)
[![Downloads](http://cranlogs.r-pkg.org/badges/ohenery?color=green)](http://www.r-pkg.org/pkg/ohenery)
[![Total](http://cranlogs.r-pkg.org/badges/grand-total/ohenery?color=green)](http://www.r-pkg.org/pkg/ohenery)
![RCpp](https://img.shields.io/badge/RCpp-inside-blue.svg)

Performs softmax regression for ordered outcomes under the [Harville](http://dx.doi.org/10.1080/01621459.1973.10482425)
and [Henery](http://dx.doi.org/10.1111/j.2517-6161.1981.tb01153.x) models.



-- Steven E. Pav, shabbychef@gmail.com

## Installation

This package can be installed 
from CRAN, 
via [drat](https://github.com/eddelbuettel/drat "drat"), or
from github:


```r
# via CRAN:
install.packages("ohenery")
# via drat:
if (require(drat)) {
    drat:::add("shabbychef")
    install.packages("ohenery")
}
# get snapshot from github (may be buggy)
if (require(devtools)) {
  install_github('shabbychef/ohenery')
}
```

## What is it?

![ohenery](tools/static/ohenery2.png)

The softmax regression generalizes logistic regression, wherein only one of two
possible outcomes is observed, to the case where one of many possible outcomes is observed.
As in logistic regression, one models the log odds of each outcome as a linear
function of some independent variables.
Moreover, a softmax regression allows one to model data from events where the number
of possible outcomes differs.

Some examples where one might apply a softmax regression:

  * Model which among five sprinters takes first place in a race,
    having observed characteristics of the runners over many races.
  * Model which film is awarded the Academy Award for Best Picture,
    based on genre information, and co-nomination information.
    (The number of nominees has varied over the years.)
  * Model which major city in the U.S. experiences the most rain
    in a given calendar year.

Note that in the examples illustrated above, one might be better served
by modeling some continuous outcome, instead of modeling the winner.
For example, one might model the speed of each racer, or the number of
votes each film garnered, or the total amount of rain each city experienced.
Discarding this information in favor of modeling the binary outcome is
likely to cause a loss of statistical power, and is 
[generally discouraged](https://statmodeling.stat.columbia.edu/2014/02/25/basketball-stats-dont-model-probability-win-model-expected-score-differential/).
However, in some cases the continuous outcome is _not_ observed, as
in the case of the Best Picture awards, or in horse racing where finishing
times are often not available.
Softmax regression can be used in these cases.

The softmax regression can be further generalized to model the case where
one observes place information about participants. For example, one might
observe which of first through fifth place each sprinter claims in each race.  
Or one might only observe some limited information, like the 
Gold, Silver and Bronze medal winners in an Olympic event.

There is more than one way to generalize the softmax to deal with ranked
outcomes like these:

  * The Harville model, where the ratio of probabilities of two participants
    taking first place is equal to the ratio of the conditional probabilities that they
    take second place, conditional on neither of them taking first place.
    Effectively in the Harville model once one has observed the first place
    finisher, the probabilities for second place simply rescale.
    I believe that if finishing times for racers are exponentially distributed,
    then their finishing places are distributed under a Harville model.
  * The Henery model generalizes the Harville model to the case where the
    ratio of probabilities of two participants
    taking first place is equal to the ratio of the conditional probabilities that they
    take second place, conditional on neither of them taking first place,
    _raised to some power_, called gamma.  The Harville model is the Henery model with all
    gamma constants equal to one.
    I believe that if finishing times for racers are (log) normally distributed, 
    then their finishing places are nearly distributed like the Henery model.


This package supports fitting softmax regressions under both models.

## Basic Usage

### Best Picture

Here we use softmax regression to model the Best Picture winner odds
as linear in some features: whether the film also was nominated for
Best Director, Best Actor or Actress, or Best Film Editing, as well
as genre designations for Drama, Romance and Comedy.
We only observe the winner of the Best Picture award, and not
runners-up, so we weight the first place finisher with a one,
and all others with a zero.
This seems odd, but it ensures that the regression does not try
to compare model differences between the runners-up.
We find the strongest absolute effect on odds from the Best Director
and Film Editing co-nomination variables.



```r
library(ohenery)
library(dplyr)
library(magrittr)
data(best_picture)
best_picture %<>%
  mutate(place=ifelse(winner,1,2)) %>%
  mutate(weight=ifelse(winner,1,0))

fmla <- place ~ nominated_for_BestDirector + nominated_for_BestActor + nominated_for_BestActress + nominated_for_BestFilmEditing + Drama + Romance + Comedy

osmod <- harsm(fmla,data=best_picture,group=year,weights=weight) 
print(osmod)
```

```
## --------------------------------------------
## Maximum Likelihood estimation
## BFGS maximization, 65 iterations
## Return code 0: successful convergence 
## Log-Likelihood: -548 
## 7  free parameters
## Estimates:
##                                   Estimate Std. error t value Pr(> t)    
## nominated_for_BestDirectorTRUE      3.2361     0.3389    9.55 < 2e-16 ***
## nominated_for_BestActorTRUE         1.0068     0.1292    7.79 6.6e-15 ***
## nominated_for_BestActressTRUE       0.0448     0.1425    0.31  0.7532    
## nominated_for_BestFilmEditingTRUE   1.9634     0.1638   11.99 < 2e-16 ***
## Drama                              -0.6598     0.2475   -2.67  0.0077 ** 
## Romance                             0.5441     0.3983    1.37  0.1719    
## Comedy                              1.3100     0.4143    3.16  0.0016 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## --------------------------------------------
##    R2: 0.77 
## --------------------------------------------
```

### Prediction

Here we use the `predict` method to get back predictions under
the Harville model produced above. 
Three different types of prediction are supported:

1. Predictions of the odds in odds space, the 'eta'.
1. Predictions of the probability of taking first place, the 'mu'.
1. The expected rank under a Harville model.

We do not currently have the ability to compute the
expected rank under the Henery model, as it is too computationally
intensive. (File an issue if this is important.)


```r
prd <- best_picture %>%
	mutate(prd_erank=as.numeric(predict(osmod,newdata=.,group=year,type='erank',na.action=na.pass))) %>%
	mutate(prd_eta=as.numeric(predict(osmod,newdata=.,group=year,type='eta'))) %>%
	mutate(prd_mu=as.numeric(predict(osmod,newdata=.,group=year,type='mu'))) 
```

### Horse Racing

The package is bundled with a dataset of three weeks of
thoroughbred and harness race results from tracks around the world.
First, let us compare the 'consensus odds', as computed from
the win pool size, with the probability of winning a race.
This is usually plotted as the empirical win probability
over horses grouped by their consensus win probability.

We present this plot below.
Clearly visible are the 'longshot bias' and 
'sureshot bias' (or favorite bias).
The longshot bias is the effect where longshots are overbet,
resulting in elevated consensus odds. 
This appears as the point off the line in the lower left. 
One can think of these as points which were moved 
_to the right_ of the plot from a point on the diagonal
by overbetting.
The sureshot bias is the complementary effect where
favorites to win are underbet, effectively shifting the
points to the left.



```r
library(ohenery)
library(dplyr)
library(magrittr)
library(ggplot2)

data(race_data)

ph <- race_data %>%
  group_by(EventId) %>%
    mutate(mu0=WN_pool / sum(WN_pool)) %>%
  ungroup() %>%
  mutate(mubkt=cut(mu0,c(0,10^seq(-2,0,length.out=14)),include.lowest=TRUE)) %>%
  mutate(tookfirst=as.numeric(coalesce(Finish==1,FALSE))) %>%
  group_by(mubkt) %>%
    summarize(winprop=mean(tookfirst),
              medmu=median(mu0),
              nrace=length(unique(EventId)),
              nhorse=n()) %>%
  ungroup() %>%
  ggplot(aes(medmu,winprop)) + 
  geom_point(aes(size=nhorse)) + 
  geom_abline(intercept=0,slope=1) + 
  scale_x_log10(labels=scales::percent) + 
  scale_y_log10(labels=scales::percent) +
  labs(x='consensus odds (from win pool)',
       y='empirical win probability',
       title='efficiency of the win pool')
print(ph)
```

<img src="tools/figure/race_data_efficiency-1.png" title="plot of chunk race_data_efficiency" alt="plot of chunk race_data_efficiency" width="700px" height="600px" />

Here we use softmax regression under the Harville model to capture
the longshot and sureshot biases.
Our model computes consensus odds as the inverse softmax
function of the consensus probabilities, which are constructed
from the win pool.
We define a factor variable that captures longshot, sureshot,
and 'vanilla' bets based on the win pool-implied probabilities.
Here we are only concerned with the win probabilities,
so we adapt a Harville model and weight only the winning
finishers.


```r
# because of na.action, make ties for fourth
df <- race_data %>%
  mutate(Outcome=coalesce(Finish,4L)) 

# create consensus odds eta0
df %<>%
  mutate(weights=coalesce(as.numeric(Finish==1),0) ) %>%
  group_by(EventId) %>%
    mutate(big_field=(n() >= 6)) %>%
    mutate(mu0=WN_pool / sum(WN_pool)) %>%
    mutate(eta0=inv_smax(mu0)) %>%
  ungroup() %>%
  dplyr::filter(big_field) %>%
  dplyr::filter(!is.na(eta0)) %>%
  mutate(bettype=factor(case_when(mu0 < 0.025 ~ 'LONGSHOT',
                                  mu0 > 0.50  ~ 'SURESHOT',
                                  TRUE        ~ 'VANILLA'))) 

# Harville Model with market efficiency
effmod <- harsm(Outcome ~ eta0:bettype,data=df,group=EventId,weights=weights)
print(effmod)
```

```
## --------------------------------------------
## Maximum Likelihood estimation
## BFGS maximization, 31 iterations
## Return code 0: successful convergence 
## Log-Likelihood: -56868 
## 3  free parameters
## Estimates:
##                      Estimate Std. error t value Pr(> t)    
## eta0:bettypeLONGSHOT  1.26078    0.03397    37.1  <2e-16 ***
## eta0:bettypeSURESHOT  1.18968    0.01450    82.1  <2e-16 ***
## eta0:bettypeVANILLA   1.11580    0.00894   124.8  <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## --------------------------------------------
##    R2: 0.53 
## --------------------------------------------
```

We see that in this model the beta coefficient for consensus odds
is around 1.12 for 'vanilla' bets,
slightly higher for sure bets, and much higher for longshots.
The interpretation is that long shots are relatively overbet,
but that otherwise bets are nearly efficient.

One would like to make a correction to inefficiencies in the 
win pool by generalizing this idea to multiple cut points along the consensus
probabilities.
However, the softmax model works by tweaking the consensus odds;
the odds translate into win probabilities in a way that depends on the
other participants, and so there is no one curve.

### Offsets

The print display of the Harville softmax regression includes an 'R-squared' 
and sometimes a 'delta R-squared'. 
The latter is only reported when an offset is used in the model.
The R-squared is the improvement in spread in the predicted ranks from the 
model compared to the null model which assumes all log odds are equal.
When an offset is given, the R-squared includes the offset,
but a 'delta' R-squared is reported which gives the improvement
in spread of predicted ranks over the model based on the offset:


```r
# Harville Model with offset
offmod <- harsm(Outcome ~ offset(eta0) + bettype,data=df,group=EventId,weights=weights)
print(offmod)
```

```
## --------------------------------------------
## Maximum Likelihood estimation
## BFGS maximization, 22 iterations
## Return code 0: successful convergence 
## Log-Likelihood: -56952 
## 2  free parameters
## Estimates:
##                 Estimate Std. error t value Pr(> t)    
## bettypeSURESHOT   0.8165     0.0522   15.65  <2e-16 ***
## bettypeVANILLA    0.4369     0.0438    9.97  <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## --------------------------------------------
##    R2: 0.1 
## delR2: -0.86 
## --------------------------------------------
```

Unfortunately, the R-squareds are hard to interpret, and can sometimes give
negative values. The softmax regression is not guaranteed to improve the
residual sum of squared errors in the rank.

### Henery Model

Here we fit a Henery model to the horse race data.
Note that in the data we only observe win, place, and show finishes, 
so we assign zero weight to all other finishers.  
Again this is so the regression does not attempt to model distinctions between runners-up.
Moreover, because of how `na.action` works on the input data,
the response variable has to be coalesced with tie-breakers.



```r
# because of na.action, make ties for fourth
df <- race_data %>%
  mutate(Outcome=coalesce(Finish,4L)) %>%
  mutate(weights=coalesce(as.numeric(Finish<=3),0) ) %>%  # w/p/s 
  group_by(EventId) %>%
    mutate(big_field=(n() >= 5)) %>%
    mutate(mu0=WN_pool / sum(WN_pool)) %>%
    mutate(eta0=inv_smax(mu0)) %>%
  ungroup() %>%
  dplyr::filter(big_field) %>%
  dplyr::filter(!is.na(eta0)) %>%
	mutate(fac_age=cut(Age,c(0,3,5,7,Inf),include.lowest=TRUE)) %>%
  mutate(bettype=factor(case_when(mu0 < 0.025 ~ 'LONGSHOT',
                                  mu0 > 0.50  ~ 'SURESHOT',
                                  TRUE        ~ 'VANILLA'))) 

# Henery Model with ...
bigmod <- hensm(Outcome ~ eta0:bettype + fac_age,data=df,group=EventId,weights=weights,ngamma=3)
print(bigmod)
```

```
## --------------------------------------------
## Maximum Likelihood estimation
## BFGS maximization, 74 iterations
## Return code 0: successful convergence 
## Log-Likelihood: -59188 
## 8  free parameters
## Estimates:
##                      Estimate Std. error t value Pr(> t)    
## fac_age(3,5]           0.2022     0.0422    4.79 1.7e-06 ***
## fac_age(5,7]           0.1702     0.0481    3.54  0.0004 ***
## fac_age(7,Inf]         0.1529     0.0514    2.98  0.0029 ** 
## eta0:bettypeLONGSHOT   1.5192     0.0388   39.12 < 2e-16 ***
## eta0:bettypeSURESHOT   1.1094     0.0223   49.84 < 2e-16 ***
## eta0:bettypeVANILLA    1.1107     0.0140   79.34 < 2e-16 ***
## gamma2                 0.7053     0.0130   54.21 < 2e-16 ***
## gamma3                 0.5275     0.0116   45.61 < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## --------------------------------------------
```

Note that the gamma coefficients are fit here as
0.71, 0.53. 
Values around 0.8 or smaller are typical.
(In a future release it would probably make sense to
make the gammas relative to the value 1,
which would make it easier to reject the Harville
model from the print summary.)

### Diving

The package is bundled with a dataset of 100 years of Olympic Men's
Platform Diving Records, originally sourced by Randi Griffin
and delivered on 
[kaggle](https://www.kaggle.com/heesoo37/120-years-of-olympic-history-athletes-and-results).

Here we convert the medal records into finishing places of 1, 2, 3 and 4 (no
medal), add weights for the fitting,
make a factor variable for age, factor the NOC (country) of the athlete.
Because Platform Diving is a subjective competition, based on scores from
judges, we investigate whether there is a 'home field advantage'
by creating a Boolean variable indicating whether the athlete is representing
the host nation.

We then fit a Henery model to the data. Note that the gamma terms come
out very close to one, indicating the Harville model would be sufficient.
The home field advantage appears real in this analysis, though modest.



```r
library(forcats)
data(diving)
fitdat <- diving %>%
  mutate(Finish=case_when(grepl('Gold',Medal)   ~ 1,
                          grepl('Silver',Medal) ~ 2,
                          grepl('Bronze',Medal) ~ 3,
                          TRUE ~ 4)) %>%
  mutate(weight=ifelse(Finish <= 3,1,0)) %>%
  mutate(cut_age=cut(coalesce(Age,22.0),c(12,19.5,21.5,22.5,25.5,99),include.lowest=TRUE)) %>%
  mutate(country=forcats::fct_relevel(forcats::fct_lump(factor(NOC),n=5),'Other')) %>%
  mutate(home_advantage=NOC==HOST_NOC)

hensm(Finish ~ cut_age + country + home_advantage,data=fitdat,weights=weight,group=EventId,ngamma=3)
```

```
## --------------------------------------------
## Maximum Likelihood estimation
## BFGS maximization, 74 iterations
## Return code 0: successful convergence 
## Log-Likelihood: -1907 
## 12  free parameters
## Estimates:
##                    Estimate Std. error t value Pr(> t)    
## cut_age(19.5,21.5]   0.0303     0.1402    0.22 0.82864    
## cut_age(21.5,22.5]  -0.7276     0.1758   -4.14 3.5e-05 ***
## cut_age(22.5,25.5]   0.0950     0.1270    0.75 0.45408    
## cut_age(25.5,99]    -0.1838     0.1377   -1.33 0.18197    
## countryGBR          -0.6727     0.2693   -2.50 0.01249 *  
## countryGER           1.0776     0.1662    6.48 8.9e-11 ***
## countryMEX           0.7159     0.1589    4.50 6.6e-06 ***
## countrySWE           0.6205     0.1853    3.35 0.00081 ***
## countryUSA           2.3202     0.1534   15.12 < 2e-16 ***
## home_advantageTRUE   0.5790     0.1377    4.20 2.6e-05 ***
## gamma2               1.0054     0.0956   10.52 < 2e-16 ***
## gamma3               0.9674     0.0993    9.75 < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## --------------------------------------------
```
