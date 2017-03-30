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


######################ANNOTATE MONOCLE DIFFERENTIALLY EXPRESSED GENE TABLE WITH CELL COUNTS #############
#NOT TESTED
#will annotate a differentially expressed gene table from monocle with cell counts and log2FC in marker-positive cells for a grouping variable in pData of monocle object, for which the DEG table was generated
#assumes two levels of grouping variable (i.e. will take the first two levels)
#cells with mincount transcripts per gene will be used to calculate group means and log2FC (mean group 1: mean group 2)
#pData of the monocle object should have a column "SampleID", the entries of which match column names of expression data
#returns annotated data frame

annotate_DEG_list<-function(monocleObject,monocleDEGtable,groupingVar,mincount){

    require(monocle)
    exp.dat.filt<-as.matrix(monocleObject)
    tabz<-data.frame(matrix(nrow=nrow(monocleDEGtable),ncol=9))
    gVu<-unique(as.character(pData(monocleObject)[,groupingVar]))
    colnames(tabz)<-c(paste0("Expr.",gVu[1]),paste0("Expr.",gVu[2]),paste0("NCells2.",gVu[1]),paste0("NCells2.",gVu[2]),"log2FC","GeneSym",paste0("Expr.",gVu[1],".iC"),paste0("Expr.",gVu[2],".iC","log2FC.iC"))
    exp.A<-exp.dat.filt[,colnames(exp.dat.filt) %in% as.character(pData(monocleObject)$SampleID)[as.character(pData(monocleObject)[,groupingVar]) %in% gVu[1]]]
    exp.B<-exp.dat.filt[,colnames(exp.dat.filt) %in% as.character(pData(monocleObject)$SampleID)[as.character(pData(monocleObject)[,groupingVar]) %in% gVu[2]]]
    for( k in seq_along(1:nrow(monocleDEGtable))){
        tabz[k,paste0("Expr.",gVu[1])]<-mean(exp.A[rownames(exp.A)==rownames(monocleDEGtable)[k],])
        tabz[k,paste0("Expr.",gVu[2])]<-mean(exp.B[rownames(exp.B)==rownames(monocleDEGtable)[k],])    
        tabz[k,paste0("NCells2.",gVu[1])]<-sum(exp.A[rownames(exp.A)==rownames(monocleDEGtable)[k],]>=mincount)
        tabz[k,paste0("NCells2.",gVu[2])]<-sum(exp.B[rownames(exp.B)==rownames(monocleDEGtable)[k],]>=mincount)    
        tmp.A<-exp.A[rownames(exp.A)==rownames(monocleDEGtable)[k],]
        tmp.B<-exp.B[rownames(exp.B)==rownames(monocleDEGtable)[k],]    
        tabz[k,paste0("Expr.",gVu[1],".iC")]<-mean(tmp.A[tmp.A>=mincount])
        tabz[k,paste0("Expr.",gVu[2],".iC")]<-mean(tmp.B[tmp.B>=mincount])
        print(paste0(k,"_processed"))
        }
    tabz[is.na(tabz)]<-0.01
    tabz$log2FC<-log2(tabz[k,paste0("Expr.",gVu[1])]/tabz[k,paste0("Expr.",gVu[2])])
    tabz$log2FC.iC<-log2(tabz[k,paste0("Expr.",gVu[1],".iC")]/tabz[k,paste0("Expr.",gVu[2],".iC")])
    tabz.an<-data.frame(v,tabz,stringsAsFactors=FALSE)
    rownames(tabz.an)<-rownames(monocleDEGtable)
    return(tabz.an)
}

################################################################################################################################

