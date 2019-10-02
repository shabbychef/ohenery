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
# Created: 2019.09.04
# Copyright: Steven E. Pav, 2019
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @title Oscar Award Best Picture Data
#' @description Historical data on the Best Picture nominees and winners
#' from 1934 through 2014. 
#' @usage data(best_picture)
#' @format A \code{data.frame} object with 484 observations and 19 columns. 
#' The columns are defined as follows:
#' \describe{
#'  \item{\code{year}}{The integer year of the nomination. These span from 1934 through 2014. Note that the
#'  number of films nominated per year varies from 5 to 12.}
#'  \item{\code{film}}{The title of the film.}
#'  \item{\code{winner}}{A logical for whether the film won the Oscar for Best
#'  Picture. There is exactly one winning film per year.}
#'  \item{\code{nominated_for_Writing}}{A logical indicating whether the film
#'  was also nominated for a Writing award that year.}
#'  \item{\code{nominated_for_BestDirector}}{A logical indicating whether the film
#'  was also nominated for Best Director award that year.}
#'  \item{\code{nominated_for_BestActress}}{A logical indicating whether the film
#'  was also nominated for at least one Best Actress award that year.}
#'  \item{\code{nominated_for_BestActor}}{A logical indicating whether the film
#'  was also nominated for at least one Best Actor award that year.}
#'  \item{\code{nominated_for_BestFilmEditing}}{A logical indicating whether the film
#'  was also nominated for at least one Best Film Editing award that year.}
#'  \item{\code{Adventure}}{A double computed as a 0/1 indicator of whether
#'  \dQuote{Adventure} was one of the genres tagged for the film in IMDb,
#'  divided by the total count of genres tagged for the film.}
#'  \item{\code{Biography}}{A double computed as a 0/1 indicator of whether
#'  \dQuote{Biography} was one of the genres tagged for the film in IMDb,
#'  divided by the total count of genres tagged for the film.}
#'  \item{\code{Comedy}}{A double computed as a 0/1 indicator of whether
#'  \dQuote{Comedy} was one of the genres tagged for the film in IMDb,
#'  divided by the total count of genres tagged for the film.}
#'  \item{\code{Crime}}{A double computed as a 0/1 indicator of whether
#'  \dQuote{Crime} was one of the genres tagged for the film in IMDb,
#'  divided by the total count of genres tagged for the film.}
#'  \item{\code{Drama}}{A double computed as a 0/1 indicator of whether
#'  \dQuote{Drama} was one of the genres tagged for the film in IMDb,
#'  divided by the total count of genres tagged for the film.}
#'  \item{\code{History}}{A double computed as a 0/1 indicator of whether
#'  \dQuote{History} was one of the genres tagged for the film in IMDb,
#'  divided by the total count of genres tagged for the film.}
#'  \item{\code{Musical}}{A double computed as a 0/1 indicator of whether
#'  \dQuote{Musical} was one of the genres tagged for the film in IMDb,
#'  divided by the total count of genres tagged for the film.}
#'  \item{\code{Romance}}{A double computed as a 0/1 indicator of whether
#'  \dQuote{Romance} was one of the genres tagged for the film in IMDb,
#'  divided by the total count of genres tagged for the film.}
#'  \item{\code{Thriller}}{A double computed as a 0/1 indicator of whether
#'  \dQuote{Thriller} was one of the genres tagged for the film in IMDb,
#'  divided by the total count of genres tagged for the film.}
#'  \item{\code{War}}{A double computed as a 0/1 indicator of whether
#'  \dQuote{War} was one of the genres tagged for the film in IMDb,
#'  divided by the total count of genres tagged for the film.}
#'  \item{\code{Other}}{A double computed as 1 minus the sum of
#'  the other genre indicators. Effectively this is is the
#'  sum of indicators for
#'  \dQuote{Mystery}, \dQuote{Family}, \dQuote{Fantasy}, 
#'  \dQuote{Action}, \dQuote{Western}, \dQuote{Music}, 
#'  \dQuote{Sport}, \dQuote{Sci Fi},  \dQuote{Film-Noir},
#'  \dQuote{Animation}, and \dQuote{Horror} divided
#'  by the total count of genres tagged for the film.}
#' }
#' @source 
#' Awards data were sourced from Wikipedia, while genre data were
#' sourced from IMDb.
#' Any errors in transcription are the fault of the package author. 
#' @template etc
#' @note \dQuote{Oscar} is a copyright property of the Academy of Motion
#' Picture Arts and Sciences. 
#' IMDb is owned by Amazon.
#' @name best_picture
#' @rdname best_picture
#' @docType data
#' @keywords data
#' @examples
#'
#' library(dplyr)
#' data(best_picture)
#'
#' best_picture %>% 
#'   group_by(nominated_for_BestDirector) %>% 
#'   summarize(propwin=mean(winner)) %>% 
#'   ungroup()
#' best_picture %>% 
#'   group_by(nominated_for_BestActor) %>% 
#'   summarize(propwin=mean(winner)) %>% 
#'   ungroup()
#' # hmmmm.
#' best_picture %>% 
#'   group_by(nominated_for_BestActress) %>% 
#'   summarize(propwin=mean(winner)) %>% 
#'  ungroup()
#'
"best_picture"

#' @title Horse Race Data
#' @description Three weeks of horse race data from tracks worldwide.
#' @usage data(race_data)
#' @format A \code{data.frame} object with 36,418 observations and 19 columns. 
#' 
#' The columns are defined as follows:
#' \describe{
#' \item{\code{EventId}}{An integer ID denoting the event (race). These range from 1 to 4486.}
#' \item{\code{TrackId}}{An integer ID number of the the track. There are 64
#' different tracks represented.}
#' \item{\code{Type}}{The type of event, one of \dQuote{Thoroughbred} or \dQuote{Harness}.}
#' \item{\code{RaceNum}}{The integer race number within a group of races at a track on a given date.}
#' \item{\code{CorrectedPostTime}}{The \sQuote{corrected} post time of the race, in the form \code{\%Y-\%m-\%d \%H:\%M:\%S},
#' presumably in the PDT time zone. Has values like \dQuote{2019-03-05 02:30:00}.}
#' \item{\code{Yards}}{The length of the race, in yards.}
#' \item{\code{SurfaceText}}{A string, one of 
#' \dQuote{Turf}, \dQuote{Dirt}, \dQuote{All-Weather} or \code{NA}.}
#' \item{\code{HorseName}}{The string name of the horse.}
#' \item{\code{HorseId}}{A unique integer ID for each horse. As different horses can have the same name, this ID is constructed from the name 
#' of the Horse, the Sire and the Dam.}
#' \item{\code{Age}}{The age of the horse, in integer years, at the time of the event.
#' Typically less than 10.}
#' \item{\code{Sex}}{A single character denoting the sex of the horse. I
#' believe the codes are 
#' \dQuote{M} for \dQuote{Mare} (female four years or older), 
#' \dQuote{G} for \dQuote{Gelding},
#' \dQuote{F} for \dQuote{Filly} (female under four years of age),
#' \dQuote{C} for \dQuote{Colt} (male under four years of age),
#' \dQuote{H} for \dQuote{Horse} (male four years of age and up),
#' \dQuote{R} for \dQuote{Rig} (hard to explain),
#' \dQuote{A} for \dQuote{???}. There are some \code{NA} values as well.}
#' \item{\code{Weight_lbs}}{The weight in integer pounds of the jockey and any equipment. Typically around 120.}
#' \item{\code{PostPosition}}{The integer starting position of the horse.
#' Typically there is a slight advantage to starting at the first or second
#' post position.}
#' \item{\code{Medication}}{One of several codes indicating any medication the horse may
#' be taking at the time of the race. I believe \dQuote{L} stands for
#' \dQuote{Lasix}, a common medication for lung conditions that is thought to give
#' horses a slight boost in speed.}
#' \item{\code{MorningLine}}{A double indicating the \dQuote{morning betting line} for
#' win bets on the horse. It is not clear how to interpret this value, perhaps
#' it is return on a dollar. Values range from 0.40 to 80.}
#' \item{\code{WN_pool}}{The total combined pool in win bets, in dollars, on
#' this horse at post time.}
#' \item{\code{PL_pool}}{The total combined pool in place bets, in dollars, on
#' this horse at post time.}
#' \item{\code{SH_pool}}{The total combined pool in show bets, in dollars, on
#' this horse at post time.}
#' \item{\code{Finish}}{The integer finishing position of the horse. A 1 means first place. We only collect values of 1, 2, and 3, while
#' the remaining finishing places are unknown and left as \code{NA}.}
#' }
#' @source 
#' Data were sourced from the web. Don't ask.
#' @template etc
#' @note The author makes no guarantees regarding correctness of this data.
#' @name race_data
#' @rdname race_data
#' @docType data
#' @keywords data
#' @examples
#'
#' library(dplyr)
#' data(race_data)
#'
#' # compute win bet efficiency
#' efficiency <- race_data %>%
#'   group_by(EventId) %>%
#'     mutate(ImpliedOdds=WN_pool / sum(WN_pool,na.rm=TRUE)) %>%
#'   ungroup() %>%
#'   mutate(OddsBucket=cut(ImpliedOdds,c(0,0.05,seq(0.1,1,by=0.10)),include.lowest=TRUE)) %>%
#'   group_by(OddsBucket) %>%
#'     summarize(PropWin=mean(as.numeric(coalesce(Finish==1,FALSE)),na.rm=TRUE),
#'               MedImpl=median(ImpliedOdds,na.rm=TRUE),
#'               nObs=n()) %>%
#'   ungroup()
#' 
#' \donttest{
#' if (require('ggplot2') && require('scales')) {
#'   efficiency %>%
#'     ggplot(aes(MedImpl,PropWin,size=nObs)) + 
#'     geom_point() + 
#'     scale_x_sqrt(labels=percent) + 
#'     scale_y_sqrt(labels=percent) + 
#'     geom_abline(slope=1,intercept=0,linetype=2,alpha=0.6) + 
#'     labs(title='actual win probability versus implied win probability',
#'          size='# horses',
#'          x='implied win probability',
#'          y='observed win probability')
#'  }
#' }
#'
"race_data"

#' @title Olympic Diving Data
#' @description One hundred years of Men's Olympic Platform Diving records.
#' @usage data(diving)
#' @format A \code{data.frame} object with 695 observations and 13 columns. 
#' 
#' The columns are defined as follows:
#' \describe{
#' \item{\code{Name}}{The participant's name.}
#' \item{\code{Age}}{The age of the participant at the time of the Olympics. Some values missing.}
#' \item{\code{Height}}{The height of the participant at the time of the Olympics, in centimeters. Many values missing.}
#' \item{\code{Weight}}{The height of the participant at the time of the Olympics, in kilograms. Many values missing.}
#' \item{\code{Team}}{The string name of the team (country) which the participant represented.}
#' \item{\code{NOC}}{The string name of the National Olympic Committee which the participant represented. This is a three character code.}
#' \item{\code{Games}}{The string name of the Olympic games, including a year.}
#' \item{\code{Year}}{The integer year of the Olympics. These range from 1906
#' through 2016.}
#' \item{\code{City}}{The string name of the host city.}
#' \item{\code{Medal}}{A string of \dQuote{Gold}, \dQuote{Silver},
#' \dQuote{Bronze} or \code{NA}.}
#' \item{\code{EventId}}{A unique integer ID for each Olympics.}
#' \item{\code{AthleteId}}{A unique integer ID for each participant.}
#' \item{\code{HOST_NOC}}{The string name of the National Olympic Committee of the nation hosting the Olympics. This is a three character code.}
#' }
#'
#' @source 
#' Data were collected by Randi Griffin from the website
#' \dQuote{sports-reference.com}, and staged on Kaggle at
#' \url{https://www.kaggle.com/heesoo37/120-years-of-olympic-history-athletes-and-results}.
#' 
#' @template etc
#' @note The author makes no guarantees regarding correctness of this data.
#' @note Please attribute this data to the upstream harvester.
#' @name diving 
#' @rdname diving 
#' @docType data
#' @keywords data
#' @examples
#'
#' library(dplyr)
#' library(forcats)
#' data(diving)
#'
#' fitdat <- diving %>%
#'   mutate(Finish=case_when(grepl('Gold',Medal)   ~ 1,
#'                           grepl('Silver',Medal) ~ 2,
#'                           grepl('Bronze',Medal) ~ 3,
#'                           TRUE ~ 4)) %>%
#'   mutate(weight=ifelse(Finish <= 3,1,0)) %>%
#'   mutate(cut_age=cut(coalesce(Age,22.0),c(12,19.5,21.5,22.5,25.5,99),include.lowest=TRUE)) %>%
#'   mutate(country=forcats::fct_relevel(forcats::fct_lump(factor(NOC),n=5),'Other')) %>%
#'   mutate(home_advantage=NOC==HOST_NOC)
#' 
#' hensm(Finish ~ cut_age + country + home_advantage,data=fitdat,weights=weight,group=EventId,ngamma=3)
#'
"diving"

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
