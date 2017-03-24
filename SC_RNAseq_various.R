#############################PLOT AND FACET FUNCTION#############################################
plotexptsne_facet<-function(scobject,genestring,genetitle,facetvar,logsc=FALSE){
    require(dplyr)
    require(scales)
    require(ggplot2)
    plotdat<-scobject@tsne
    plotdat$Var<-facetvar
    gi<-genestring %in% rownames(scobject@ndata)
    if(sum(gi)<length(gi)){print(paste(sum(gi),"out of",length(gi),"genes are match the expression dataset rownames",sep=" "))}
    genestring<-genestring[gi]
    l<-apply(scobject@ndata[genestring,]-.1,2,sum)+.1
    if (logsc) {
              f <- l == 0
              l <- log2(l)
              l[f] <- NA
            }
    plotdat$label<-l
    ggplot(plotdat %>% arrange(label), aes(x = V1, y = V2, color = label))+geom_point(size = 2)+scale_colour_continuous(low = "steelblue3", high ="darkorange", space = "Lab", na.value = "grey50", 
      guide = "colourbar",name=ifelse(logsc==FALSE,"Counts","Log2Counts"))+facet_wrap(~Var,ncol=2)+xlab("Dim1")+ylab("Dim2")+theme(axis.text=element_text(size=14),axis.title=element_text(size=16),strip.text=element_text(size=14))+ggtitle(genetitle)
    }

########################################################################################################
