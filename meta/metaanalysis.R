#import dat
# add study number variable
meta$studynumber = c(1:6)
# Rename study column
# http://www.cookbook-r.com/Manipulating_data/Renaming_columns_in_a_data_frame/
names(meta)[names(meta)=="study"] <- "Author"
# successfully renamed
library(metafor)
# calculation of proportion meta analysis by proportional log odds method
# http://www.metafor-project.org/doku.php/analyses:stijnen2010
res <- rma(measure="PLO", xi= event, ni= total, data=meta)

res <- rma(measure="PLO", xi= meta$no..of.positive.cultures, ni= meta$sample.size, data=meta)
res2 <- rma(measure="PLO", xi= meta$no..of.gram.positives, ni= meta$no..of.positive.cultures, data=meta)
res3 <- rma(measure="PLO", xi= meta$no..of.gram.negatives, ni= meta$no..of.positive.cultures, data=meta)
print(respft, digits=3)
print(res2,digits = 3)
print(res3,digits=3)

respR = rma(measure="PR", xi= meta$no..of.positive.cultures, ni= meta$sample.size, data=meta)
respft = rma(measure="PFT", xi= meta$no..of.positive.cultures, ni= meta$sample.size, data=meta, method ="DL")

# Back transform log odds ratio

predict(res, transf=transf.ilogit, digits=3)

predict(respft, transf = transf.pft, digits =3)



# ci.lb,ci.ub(fixed effect), cr.lb and cr.ub(random effects)

# studying structure of res, to decide which plots to apply
str(res)
# Thus rma.uni model
metafor::plot.rma.uni(res)

# simple forest plot
forest(res)

# transformed effect , we will use forest aruments

forest(res,transf = transf.logit())

# employing Binomial normal model

res1 <- rma.glmm(measure="PLO", xi= event, ni= total, data=meta)
print(res1, digits=3)

predict(res1, transf=transf.ilogit, digits=3)

# escalc
resescalc = escalc(measure="PLO", xi= event, ni= total, data=meta)


PLO.ma<-rma.uni(yi,vi,method="REML",data=resescalc)
summary(PLO.ma)

# forest plots
forest.rma(PLO.ma)

# https://ecologyforacrowdedplanet.wordpress.com/2013/05/10/using-metafor-and-ggplot-togetherpart-1/

# using ggplot2

library(ggplot2)
ROM.ma =PLO.ma

theme_set(theme_bw(base_size=10))
forrest_data<-rbind(data.frame(ES=ROM.ma$yi,SE=sqrt(ROM.ma$vi),Type="Study",Study= meta$Author),data.frame(ES=ROM.ma$b,SE=ROM.ma$se,Type="Summary",Study="Summary"))
forrest_data$Study2<-factor(forrest_data$Study, levels=rev(levels(forrest_data$Study)) )
levels(forrest_data$Study2)
plot1<-ggplot(data=forrest_data,aes(x= Study2,y= ES,ymax=ES+(1.96*SE),ymin=ES-(1.96*SE),size=factor(Type),colour=factor(Type)))+geom_pointrange()
plot2<-plot1+coord_flip()+geom_hline(aes(x=0), lty=2,size=1)+scale_size_manual(values=c(0.5,1))
plot3<-plot2+xlab("Study")+ylab("log response ratio")+scale_colour_manual(values=c("grey","black"))
plot3+theme(legend.position="none")
plot3
plot1
plot2






