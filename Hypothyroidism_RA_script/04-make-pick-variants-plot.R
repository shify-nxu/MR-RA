library(ggplot2)
library(reshape2)
library(aplot)

sei_pick_variants <- read.table("Merge-SNP-singleanalysis-Sei_pick_vairiants.txt",header=T,sep='\t',stringsAsFactors=F,check.names=FALSE)
sei_pick_variants_use <- melt(sei_pick_variants)
colnames(sei_pick_variants_use) <- c('id','Features','value')
sei_pick_variants_use_score <- sei_pick_variants_use[sei_pick_variants_use$Features != "b_score",]
b_score <- sei_pick_variants_use[sei_pick_variants_use$Features == "b_score",]

pp <- ggplot(data = sei_pick_variants_use_score,aes(x=id,y=Features,color=value))+
    geom_tile(color="grey70",fill="white",size=0.6)+
    geom_point(shape=19)+
    geom_point(aes(size=abs(value)),show.legend = F)+
    labs(x = NULL,y = NULL,color="Sei score") + 
    scale_color_viridis_c(option = "plasma")+
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0),position="right")+
    theme(text=element_text(family="Arial"),
          axis.text.x=element_text(color="black",angle=45,hjust=1),
          axis.text.y=element_text(color="black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y=element_blank(),
          panel.border = element_rect(fill=NA,color="grey80",size=1, linetype="solid"))+
    scale_size(range=c(1,5),guide=NULL)+
    guides(color=guide_colorbar(direction = "vertical",
    reverse = F,barwidth = unit(.6, "cm"),
    barheight = unit(10,"cm")))
                         

b_score_plot <- ggplot(data = b_score, aes(x=id, y=Features ,fill=value))+
    geom_tile(color = "white",lwd = 1.5,linetype = 1)+
    scale_fill_gradient(low = "blue", high = "red") +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
	  panel.background = element_blank(),# Make the background white
	  panel.border = element_blank())+# Make the border with none
          labs(fill= "b score")

pp1 <- pp %>% insert_top(b_score_plot,height=0.03)

ggsave(pp1, file="sei_pick_variants_use.png", width = 8, height=7, dpi=900)
