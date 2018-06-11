

#' group_slice
#' 
#' get N groups.
#' #' 
#' @param dt - a grouped tbl
#' @param v - a vector for slicing
#' @examples e.g. - an example
#' this 
#' @return returns 
#' @export 

group_slice<-function(dt,v){
  stopifnot('data.frame'%in%class(dt))
  stopifnot(v>0,(v%%1)==0)
  if(n_groups(dt)==1) warning('Only 1 group to slice')
  #get the grouping variables
  groups=
    unique(dt%>%select())%>%
    group_by()%>%
    slice(v)
  out=inner_join(dt,groups)
  out
}