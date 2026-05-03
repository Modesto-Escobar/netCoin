## timeCoin----
timeCoin <- function(periods, name = "name", start = "start", end = "end", group = NULL,
                     text = NULL, main = NULL, note = NULL, info = NULL, 
                     events = NULL, eventNames = "name", eventPeriod = "period", eventTime = "date",
                     eventColor = NULL, eventShape = NULL,
                     cex = 1, language = c("en","es","ca"), dir = NULL){
  time <- timeline_rd3(periods, name = name, start = start, end = end,
                       group = group, text = text, main = main, note = note, info = info,
                       events = events, eventNames = eventNames, eventPeriod = eventPeriod,
                       eventTime = eventTime, eventColor = eventColor, eventShape = eventShape,
                       cex = cex, language = language, dir = dir)
  class(time) <- c("timeCoin",class(time))
  return(time)
}

summary.timeCoin <- function(object, ...){
  cat(dim(object$periods)[1], "categories.\n")
  cat(object$options$start, "'s distribution:","\n",sep="")
  print(summary(object$periods[[object$options$start]]))
  cat(object$options$end, "'s distribution:","\n",sep="")
  print(summary(object$periods[[object$options$end]]))
}

## mobileEdges ----
mobileEdges<-function(data, name=1, number=2, difference=0) {
  if(!is.numeric(data[[number]])) data[[number]]<-as.numeric(paste(data[[number]]))
  DC<-matrix(NA,nrow=nrow(data),ncol=nrow(data))
  colnames(DC)<-rownames(DC)<-data[[name]]
  for(i in 1:nrow(data))DC[i,]=ifelse(abs(data[[number]][i]-t(data[[number]]))<=difference,(1+difference-abs(data[[number]][i]-t(data[[number]]))),0)
  diag(DC)<-0
  DCLinks<-edgeList(DC,"shape",min=1)
  colnames(DCLinks)[3]<-"sim."
  DCLinks$dist.<-(1+difference-DCLinks$sim.)
  return(DCLinks)
}

## incTime ----
incTime<-function(data, name="name", beginning="birth", end="death") {
  if(!is.integer(data[[beginning]])) data[[beginning]]<-as.integer(paste(data[[beginning]]))
  if(!is.integer(data[[end]])) data[[end]]<-as.integer(paste(data[[end]]))
  anos<-min(na.omit(data[[beginning]])):max(na.omit(data[[end]]))
  E<-matrix(NA,nrow=nrow(data),ncol=length(anos))
  colnames(E)<-anos
  for(i in 1:nrow(data)) E[i,]<-ifelse(anos>=data[[beginning]][i] & (anos<=data[[end]][i] | is.na(data[[end]][i])),1,0)
  Datos<-as.data.frame(t(E))
  colnames(Datos)<-data[[name]]
  return(Datos)
}
