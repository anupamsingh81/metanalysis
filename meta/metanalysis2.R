meta = ER.SHEET
# metanalysis of proportions

# http://www.metafor-project.org/doku.php/analyses:miller1978


meta <- escalc(measure="PFT", xi= meta$no..of.positive.cultures, ni= meta$sample.size, data= meta, add=0)
res <- rma(yi, vi, method="DL", data= meta)
predictRE = predict(res, transf=transf.ipft.hm, targs=list(ni= meta$sample.size))

# Forestplot
dat =meta
pred= predictRE
dat.back <- summary(dat, transf=transf.ipft, ni=dat$sample.size)


forest(dat.back$yi, ci.lb=dat.back$ci.lb, ci.ub=dat.back$ci.ub, psize=1,slab= paste(dat.back$study, sep=", "),
       refline=NA, digits=3, xlab="Proportion")
addpoly(pred$pred, ci.lb=pred$ci.lb, ci.ub=pred$ci.ub, row=- 0.5, digits=3, mlab="FE Model", efac=1.3)

text(-0.5, 9, "Study",               pos=4)
text( 1.4, 9, "Proportion [95% CI]", pos=2)




resfe = rma(yi, vi, method="FE", data= meta)
predictFE =predict(resfe, transf=transf.ipft.hm, targs=list(ni= meta$sample.size))
predictFE

predictfinal = c(predictFE,predictRE)
predictfinal