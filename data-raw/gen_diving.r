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
# Created: 2019.09.24
# Copyright: Steven E. Pav, 2019
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

suppressMessages(library(docopt))       # we need docopt (>= 0.3) as on CRAN

doc <- "Usage: gen_diving.r [-v] [-L <LOC>] [-O <OUTFILE>]

-L LOC --loc=LOC                 Location of olympics data [default: ~/github/ohenery/data-raw/athlete_events.csv.bz2]
-O OUTFILE --outfile=OUTFILE     Location of output data [default: diving.csv]
-v --verbose                     Be more verbose
-h --help                        show this help text"

opt <- docopt(doc)

suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(magrittr)
})

indat <- readr::read_csv(opt$loc,col_types=cols(ID = col_integer(),
                                                Age = col_double(),
                                                Height = col_double(),
                                                Weight = col_double(),
                                                Year = col_integer(),
                                                .default = col_character()
                                                )) 

lookup <- tibble::tribble(~City,   ~HOST_NOC,
                          "Athina",  "GRE",
                          "London", "GBR",
                          "Stockholm", "SWE",
                          "Antwerpen", "BEL",
                          "Paris", "FRA",
                          "Amsterdam", "NED",
                          "Los Angeles", "USA",
                          "Berlin", "GER",
                          "Helsinki", "FIN",
                          "Melbourne", "AUS",
                          "Roma", "ITA",
                          "Tokyo", "JPN",
                          "Mexico City", "MEX",
                          "Munich", "GER",
                          "Montreal", "CAN",
                          "Moskva", "RUS",
                          "Seoul", "KOR",
                          "Barcelona", "ESP",
                          "Atlanta", "USA",
                          "Sydney", "AUS",
                          "Beijing", "CHN",
                          "Rio de Janeiro", "BRA")

diving <- indat %>%
  filter(grepl('Diving Men.s Platform',Event)) %>%
  arrange(Year) %>%
  mutate(EventId=group_indices(.,Games)) %>%
  mutate(AthleteId=group_indices(.,Name)) %>%
  select(-ID,-Sex,-Sport,-Event,-Season) %>%
  arrange(EventId,Medal) %>%
  left_join(lookup,by='City')

diving %>% readr::write_csv(opt$outfile)
usethis::use_data(diving,overwrite=TRUE)

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
