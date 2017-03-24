#############################PLOT AND FACET FUNCTION#############################################
#Extends functionality of plotexptsne from RaceID(2) package to include faceting by a categorical variable of interest
#scobject is SCseq object as from RaceID(2)
#genestring is character string matching rows of ndata slot in scobject to plot gene expression for
#genetitle is character string to use as plot tile
#facetvar is a character string with a categorical variable describing the cells, in the same order as rownames of scobject@tsne
#logsc is a logical string to switch between direct count plotting and log2 count plotting
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
