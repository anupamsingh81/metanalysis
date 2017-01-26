study = studynu
study = study[1:95,]
str(study$culture)
str(study)
# replacing "" with value character http://stackoverflow.com/questions/3357743/replacing-character-values-with-na-in-a-data-frame
str(study$culture)
# change factor to character then only blank or NA replacement done,important
study$culture = as.character(study$culture)
study$culture[is.na(study$culture) ] = "Negative"
# Reconvert to factor
study$culture = as.factor(study$culture)
summary(study$culture)
plot(study$culture)

labelname = colnames(study)
labelname[13]

main2 = paste("Pie Chart of ", labelname[13],"sensitivity with sample sizes")
# http://stackoverflow.com/questions/24901061/in-r-how-do-i-compute-factors-percentage-given-on-different-variable
# http://www.statmethods.net/graphs/pie.html
# Pie Chart from data frame with Appended Sample Sizes
mytable <- table(study$amikacin)
mytable
pct =prop.table(mytable)
pct = as.data.frame(pct)
pct1 = round(100*pct[,2],2)

lbls <- paste(names(mytable)," ",pct1,"%","\n", mytable, sep="")
pie(mytable, labels = lbls,col=rainbow(length(lbls)),
    #main="Pie Chart of Culture\n (with sample sizes)"
    main = main2)

fpie = function(x , max =68) {
  
  
  mytable<- table(study[,x])
  pct =prop.table(mytable)
  pct = as.data.frame(pct)
  pct1 = round(100*pct[,2],2)
  labelname = colnames(study)
  main2 = paste("Pie Chart of ", labelname[x],"sensitivity with sample sizes")
  lbls <- paste(names(mytable)," ",pct1,"%","\n", mytable, sep="")
  main3 = paste("'",x,".jpg'")
  jpeg(main3)
  p =  pie(mytable, labels = lbls,col=rainbow(length(lbls)),
      #main="Pie Chart of Culture\n (with sample sizes)"
      main = main2)
  dev.off()
  
}

fpie(15)


fpie(14:19)


lapply(study,fpie)


# ggplot

library(ggplot2)
ggplot(study , aes(x=factor(1), fill= study[,14]))+
  geom_bar(width = 1)+
  coord_polar("y")




# 