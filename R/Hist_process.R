#' Divide sample size to four parts
#'
#' @description Divide sample size to four equal parts for histgram plotting the SV distribution
#'
#' @param hist A numeric vector contains SV number counted per sample (named by sample id).
#'
#' @return A list with four elements (i.e. \code{'section1'}, \code{'section2'}, \code{'section3'} and \code{'section4'} - each represents a SV distribution in one-fourth of samples)
#'
#' @export
histgram_process <- function(hist) {
	
	size = length(hist);
	average_num = floor(size/4);
	yushu = size %% 4;

	section1 = NULL;	section2 = NULL;	section3 = NULL;	section4 = NULL; 
	if ( yushu == 0 ) {
		if ( average_num > 0 ) {
			section1 = hist[1:average_num]; section2 = hist[(average_num+1):(average_num*2)];
			section3 = hist[(average_num*2+1):(average_num*3)]; section4 = hist[(average_num*3+1):(average_num*4)];
		}
	} else if ( yushu == 1 ) {
		if ( average_num == 0 ) {
			section1 = hist[1];
		} else {
			section1 = hist[1:(average_num+1)]; section2 = hist[(average_num+2):(average_num*2+1)];
			section3 = hist[(average_num*2+2):(average_num*3+1)]; section4 = hist[(average_num*3+2):(average_num*4+1)];
		}
	} else if ( yushu == 2 ) {
		if ( average_num == 0 ) {
			section1 = hist[1]; section2 = hist[2];
		} else {
			section1 = hist[1:(average_num+1)]; section2 = hist[(average_num+2):(average_num*2+2)];
			section3 = hist[(average_num*2+3):(average_num*3+2)]; section4 = hist[(average_num*3+3):(average_num*4+2)];
		}
	} else if ( yushu == 3 ) {
		if ( average_num == 0 ) {
			section1 = hist[1]; section2 = hist[2]; section3 = hist[3];
		} else {
			section1 = hist[1:(average_num+1)]; section2 = hist[(average_num+2):(average_num*2+2)];
			section3 = hist[(average_num*2+3):(average_num*3+3)]; section4 = hist[(average_num*3+4):(average_num*4+3)];
		}
	}
	tmp = list(section1=section1, section2=section2, section3=section3, section4=section4)
	return(tmp)
}

