# /usr/bin/r
#
# Copyright 2019-2019 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav 
#
# This file is part of ohenery.
#
# ohenery is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ohenery is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with ohenery.  If not, see <http://www.gnu.org/licenses/>.
#
# Created: 2019.09.03
# Copyright: Steven E. Pav, 2019
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

suppressMessages(library(docopt))       # we need docopt (>= 0.3) as on CRAN

doc <- "Usage: gen_best_picture.r [-v] [-L <LOC>] [-O <OUTFILE>]

-L LOC --loc=LOC                 Location of upstream data [default: ~/github/shabbychef.github.io/content/data/best_picture_2.csv]
-O OUTFILE --outfile=OUTFILE     Location of output data [default: best_picture.csv]
-v --verbose                     Be more verbose
-h --help                        show this help text"

opt <- docopt(doc)

suppressMessages({
	library(readr)
	library(dplyr)
	library(tidyr)
	library(usethis)
})

inp <- readr::read_csv(opt$loc)

#best_picture <- inp %>%
  #filter(category=='Best Picture') %>%
  #select(film,year,matches('^nominated_for_(Writing|Best(Director|Actor|Actress|FilmEditing))$'),
         #winner,
         #Action, Adventure, Animation, Biography, Comedy, Crime, Drama,
         #Family, Fantasy, `Film-Noir`, History, Horror, Music, Musical,
         #Mystery, Romance, `Sci-Fi`, Sport, Thriller, War, Western) 

#best_picture %>%
  #select(Action, Adventure, Animation, Biography, Comedy, Crime, Drama,
         #Family, Fantasy, `Film-Noir`, History, Horror, Music, Musical,
         #Mystery, Romance, `Sci-Fi`, Sport, Thriller, War, Western) %>%
  #tidyr::gather(key=series,value=value) %>%
  #group_by(series) %>%
    #summarize(totweight=sum(value)) %>%
  #ungroup() %>%
  #arrange(totweight) %>%
  #mutate(sumtot=cumsum(totweight)) %>%
  #arrange(desc(totweight)) %>%
  #knitr::kable()

#|series    | totweight|  sumtot|
#|:---------|---------:|-------:|
#|Drama     |   240.633| 595.000|
#|Romance   |    77.600| 354.367|
#|Comedy    |    45.150| 276.767|
#|Biography |    33.483| 231.617|
#|War       |    24.717| 198.133|
#|Crime     |    23.350| 173.417|
#|Adventure |    22.850| 150.067|
#|History   |    21.117| 127.217|
#|Thriller  |    18.450| 106.100|
#|Musical   |    12.367|  87.650|
#|Mystery   |    11.733|  75.283|
#|Family    |    10.583|  63.550|
#|Fantasy   |     9.700|  52.967|
#|Action    |     9.500|  43.267|
#|Western   |     8.900|  33.767|
#|Music     |     8.117|  24.867|
#|Sport     |     6.317|  16.750|
#|Sci-Fi    |     4.167|  10.433|
#|Film-Noir |     3.817|   6.267|
#|Animation |     1.450|   2.450|
#|Horror    |     1.000|   1.000|

best_picture <- inp %>%
  filter(category=='Best Picture') %>%
	distinct(film,year,winner,.keep_all=TRUE) %>%
  select(film,year,matches('^nominated_for_(Writing|Best(Director|Actor|Actress|FilmEditing))$'),
         winner,
         Action, Adventure, Animation, Biography, Comedy, Crime, Drama,
         Family, Fantasy, `Film-Noir`, History, Horror, Music, Musical,
         Mystery, Romance, `Sci-Fi`, Sport, Thriller, War, Western) %>%
  mutate(Other= Mystery + Family + Fantasy + Action + Western + Music + Sport + `Sci-Fi` + `Film-Noir` + Animation + Horror) %>%
  select(-Mystery, -Family, -Fantasy, -Action, -Western, -Music, -Sport, -`Sci-Fi`, -`Film-Noir`, -Animation, -Horror) %>%
  arrange(year,desc(winner)) %>%
	select(year,film,winner,everything())

best_picture %>% readr::write_csv(opt$outfile)
usethis::use_data(best_picture,overwrite=TRUE)

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
