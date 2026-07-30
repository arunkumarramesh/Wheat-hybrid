library(data.table)
library(ggplot2)


sample_cols <- c(CS="#0072B2",CSxP="#E69F00",P="#CC79A7")

summarize_cds_methylation <- function(infile,remove_third_col=FALSE) {
  dt <- fread(infile)
  if(remove_third_col) dt <- dt[,-3,with=FALSE]
  dt <- dt[cov_CS>10 & cov_CSxP>10 & cov_P>10]
  dt[,subgenome:=sub("^chr[0-9]+([ABD])$","\\1",chr)]
  dt <- dt[subgenome %in% c("A","B","D")]
  meth_long <- rbind(dt[,.(subgenome,sample="CS",pct=pct_CS)],dt[,.(subgenome,sample="CSxP",pct=pct_CSxP)],dt[,.(subgenome,sample="P",pct=pct_P)])
  meth_summary <- meth_long[,.(N=sum(!is.na(pct)),mean_pct=mean(pct,na.rm=TRUE),sd_pct=sd(pct,na.rm=TRUE)),by=.(subgenome,sample)]
  meth_summary[,se:=sd_pct/sqrt(N)]
  meth_summary[,ci:=qt(0.975,df=N-1)*se]
  meth_summary[,ci_low:=mean_pct-ci]
  meth_summary[,ci_high:=mean_pct+ci]
  meth_summary
}

cds_meth_summary <- rbindlist(list(CG=summarize_cds_methylation("merged_CG_symmetric_CDS.txt.gz"),CHG=summarize_cds_methylation("merged_CHG_symmetric_CDS.txt.gz"),CHH=summarize_cds_methylation("merged_CHH_all_CDS.txt.gz",remove_third_col=TRUE)),idcol="context")

genome_meth_summary <- rbindlist(list(CG=fread("met_classify_cg.tsv"),CHG=fread("met_classify_chg.tsv"),CHH=fread("met_classify_chh.tsv")),idcol="context")

meth_summary_combined <- rbindlist(list("Genome-wide"=genome_meth_summary,"CDS"=cds_meth_summary),idcol="region",fill=TRUE)
meth_summary_combined[,context:=factor(context,levels=c("CG","CHG","CHH"))]
meth_summary_combined[,sample:=factor(sample,levels=c("CS","CSxP","P"))]
meth_summary_combined[,subgenome:=factor(subgenome,levels=c("A","B","D"))]
meth_summary_combined[,region:=factor(region,levels=c("Genome-wide","CDS"))]
meth_summary_combined[,x_group:=factor(paste(region,subgenome,sep="."),levels=c("Genome-wide.A","Genome-wide.B","Genome-wide.D","CDS.A","CDS.B","CDS.D"))]

p_methylation <- ggplot(meth_summary_combined,aes(x=x_group,y=mean_pct,fill=sample,group=sample)) +
  geom_vline(xintercept=3.5,linetype="dashed",linewidth=0.3,colour="grey50") +
  geom_col(position=position_dodge(width=0.75),width=0.65,colour="black",linewidth=0.25) +
  geom_errorbar(aes(ymin=ci_low,ymax=ci_high),position=position_dodge(width=0.75),width=0.2,colour="black",linewidth=0.4) +
  scale_fill_manual(values=sample_cols) +
  scale_x_discrete(labels=c("Genome-wide.A"="A","Genome-wide.B"="B","Genome-wide.D"="D","CDS.A"="A","CDS.B"="B","CDS.D"="D")) +
  facet_wrap(~context,nrow=1,scales="free_y") +
  labs(x=NULL,y="Mean methylation (%)",fill=NULL) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(linewidth=0.2,colour="grey90"),strip.background=element_rect(fill="white",colour="black"),axis.title=element_text(face="bold"),axis.text.x=element_text(size=10),legend.position="right",plot.margin=margin(5.5,5.5,28,5.5)) +
  annotate("text",x=2,y=-Inf,label="Genome-wide",vjust=3.2,size=3.5,fontface="bold") +
  annotate("text",x=5,y=-Inf,label="CDS",vjust=3.2,size=3.5,fontface="bold") +
  coord_cartesian(clip="off")

ggsave("percent_meth_genomewide_CDS.pdf",p_methylation,width=10,height=3)

coverage_summary <- rbindlist(list(CG=fread("mean_coverage_cg.tsv"),CHG=fread("mean_coverage_chg.tsv"),CHH=fread("mean_coverage_chh.tsv")),idcol="context")

coverage_summary[,context:=factor(context,levels=c("CG","CHG","CHH"))]
coverage_summary[,sample:=factor(sample,levels=c("CS","CSxP","P"))]
coverage_summary[,subgenome:=factor(subgenome,levels=c("A","B","D"))]

p_coverage <- ggplot(coverage_summary,aes(x=subgenome,y=mean_cov,fill=sample,group=sample)) +
  geom_col(position=position_dodge(width=0.7),width=0.65,colour="black",linewidth=0.25) +
  geom_errorbar(aes(ymin=ci_low,ymax=ci_high),position=position_dodge(width=0.7),width=0.2,colour="black",linewidth=0.4) +
  scale_fill_manual(values=sample_cols) +
  facet_wrap(~context,nrow=1) +
  labs(x=NULL,y="Mean coverage",fill=NULL) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(linewidth=0.2,colour="grey90"),strip.background=element_rect(fill="white",colour="black"),axis.title=element_text(face="bold"),axis.text.x=element_text(angle=0,hjust=0.5),legend.position="right")

ggsave("coverage_meth.pdf",p_coverage,width=4.5,height=2.5)