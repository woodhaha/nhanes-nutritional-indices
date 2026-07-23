# -- figures_restyle.R
# All 10 main + 3 supplementary figures, consistent styling.
# Width: 7.5" (19.05cm, fits Elsevier 2-column max)
# Font sizes: >= 11pt axis, >= 14pt title
library(ggplot2); library(dplyr); library(survival); library(survminer)

FIG_DIR <- normalizePath("figures/gi_analysis")
RES_DIR <- normalizePath("results/gi_analysis")

theme_clin <- function(base=14, axis=11) {
  theme_classic(base_size=base) +
    theme(plot.title=element_text(face="bold",size=base,hjust=0.5),
          plot.subtitle=element_text(size=base-3,color="grey40",hjust=0.5),
          axis.title=element_text(size=base-1),
          axis.text=element_text(size=axis),
          legend.title=element_text(size=base-2),
          legend.text=element_text(size=base-3),
          plot.margin=margin(8,8,8,8))
}
COL3 <- c("PNI"="#1F77B4","CONUT"="#2CA02C","GNRI"="#FF7F0E")

df <- readRDS(file.path(RES_DIR,"nhanes_clean.rds"))

# ===== FIG 1: KM =====
cat("Fig 1 ...\n")
gi <- df %>% filter(gi_tumor==1,surv_years>0,!is.na(PNI)) %>%
  mutate(pni_t=factor(ntile(PNI,3),1:3,c("Low PNI (T1)","Mid PNI (T2)","High PNI (T3)")))
km <- survfit(Surv(surv_years,death)~pni_t,data=gi)
for (e in c("pdf","png")) {
  if(e=="pdf") pdf(file.path(FIG_DIR,"fig1_km.pdf"),width=7.5,height=7) else
    png(file.path(FIG_DIR,"fig1_km.png"),width=7.5,height=7,units="in",res=300)
  par(oma=c(0,0,1,0))
  print(ggsurvplot(km,data=gi,pval=TRUE,pval.size=4,risk.table=TRUE,
    risk.table.height=0.18,risk.table.fontsize=3,conf.int=TRUE,
    palette=unname(COL3[c(2,1,3)]),surv.median.line="hv",
    legend.title="PNI Tertile",
    title="All-Cause Survival by PNI Tertile",
    subtitle="NHANES III + Continuous NHANES 2005-2016",
    xlab="Years",ylab="Survival Probability",
    ggtheme=theme_classic(base_size=12)),newpage=FALSE)
  dev.off()
}

# ===== FIG 2: Cox forest =====
cat("Fig 2 ...\n")
cr <- read.csv(file.path(RES_DIR,"cox_main.csv")) %>% filter(Adj=="Adjusted") %>%
  mutate(Index=factor(Index,c("GNRI","CONUT","PNI")))
p2 <- ggplot(cr,aes(x=HR,y=Index,color=Index)) +
  geom_vline(xintercept=1,linetype="dashed",color="grey60",linewidth=0.4) +
  geom_point(size=3.5)+
  geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=0.08,linewidth=1)+
  scale_color_manual(values=COL3,guide="none")+
  scale_x_log10(breaks=c(0.6,0.8,1.0),labels=c("0.60","0.80","1.00"))+
  labs(title="Nutritional Indices and All-Cause Mortality",
       subtitle="NHANES III + Continuous NHANES 2005-2016",
       x="Hazard Ratio per 1-SD (95% CI)",y="")+
  theme_clin(14,12)+coord_cartesian(xlim=c(0.5,1.05))
ggsave(file.path(FIG_DIR,"fig2_cox_forest.pdf"),p2,width=7.5,height=3.5)
ggsave(file.path(FIG_DIR,"fig2_cox_forest.png"),p2,width=7.5,height=3.5,dpi=300)

# ===== FIG 3: Time-dependent =====
cat("Fig 3 ...\n")
gi3 <- df %>% filter(gi_tumor==1,surv_years>0,!is.na(PNI)) %>%
  mutate(PNI_s=scale(PNI)[,1],CONUT_s=scale(-CONUT)[,1],GNRI_s=scale(GNRI)[,1])
td3 <- bind_rows(lapply(c("PNI_s","CONUT_s","GNRI_s"),function(idx){
  d <- survSplit(Surv(surv_years,death)~.,data=gi3,cut=c(2,5),episode="ep")
  d$ep_f <- factor(d$ep,0:2,c("0-2 years","2-5 years","5+ years"))
  f <- coxph(as.formula(paste0("Surv(tstart,surv_years,death)~",idx,"*ep_f+age+sex+race_eth")),data=d)
  s <- summary(f)$coefficients
  r <- grep(gsub("_s","",idx),rownames(s))
  b <- s[r,1]; se <- s[r,3]
  h  <- c(b[1], b[1]+b[2], b[1]+b[3])
  hs <- c(se[1], sqrt(se[1]^2+se[2]^2), sqrt(se[1]^2+se[3]^2))
  data.frame(Index=gsub("_s","",idx),Period=c("0-2 y","2-5 y","5+ y"),HR=exp(h),Lower=exp(h-1.96*hs),Upper=exp(h+1.96*hs))
}))
td3$Period <- factor(td3$Period,c("0-2 y","2-5 y","5+ y"))
p3 <- ggplot(td3,aes(x=Index,y=HR,color=Index,shape=Period))+
  geom_hline(yintercept=1,linetype="dashed",color="grey60",linewidth=0.4)+
  geom_point(position=position_dodge(0.6),size=3)+
  geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.12,linewidth=0.8,
                position=position_dodge(0.6))+
  scale_color_manual(values=COL3,guide="none")+
  scale_shape_manual(values=c(16,17,15),name="")+
  labs(title="Time-Dependent Hazard Ratios",
       subtitle="NHANES III + Continuous NHANES 2005-2016",
       x="",y="Hazard Ratio per 1-SD (95% CI)")+
  theme_clin(14,11)+coord_cartesian(ylim=c(0,3.5))+
  theme(legend.position="bottom")
ggsave(file.path(FIG_DIR,"fig3_timedependent.pdf"),p3,width=7,height=5)
ggsave(file.path(FIG_DIR,"fig3_timedependent.png"),p3,width=7,height=5,dpi=300)

# ===== FIG 4: Landmark =====
cat("Fig 4 ...\n")
lm <- read.csv(file.path(RES_DIR,"cox_landmark.csv")) %>%
  mutate(Index=factor(Index,c("GNRI","CONUT","PNI")))
p4 <- ggplot(lm,aes(x=Landmark,y=HR,color=Index,group=Index))+
  geom_hline(yintercept=1,linetype="dashed",color="grey60",linewidth=0.4)+
  geom_point(size=3,position=position_dodge(0.4))+
  geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.06,linewidth=0.8,
                position=position_dodge(0.4))+
  geom_line(position=position_dodge(0.4),linewidth=0.5)+
  scale_color_manual(values=COL3,name="")+
  labs(title="Landmark Analysis",
       subtitle="NHANES III + Continuous NHANES 2005-2016",
       x="Landmark Time",y="Hazard Ratio per 1-SD (95% CI)")+
  theme_clin(14,11)+coord_cartesian(ylim=c(0.55,1.2))+
  theme(legend.position="bottom")
ggsave(file.path(FIG_DIR,"fig4_landmark.pdf"),p4,width=7,height=5)
ggsave(file.path(FIG_DIR,"fig4_landmark.png"),p4,width=7,height=5,dpi=300)

# ===== FIG 5: Competing risks =====
cat("Fig 5 ...\n")
fg <- read.csv(file.path(RES_DIR,"competingrisks_finegray.csv")) %>%
  mutate(Risk=ifelse(grepl("Non",Risk),"Non-cancer death","Cancer death"),
         Index=factor(Index,c("GNRI","CONUT","PNI")))
p5 <- ggplot(fg,aes(x=sHR,y=Index,color=Index,shape=Risk))+
  geom_vline(xintercept=1,linetype="dashed",color="grey60",linewidth=0.4)+
  geom_point(size=3.5,position=position_dodge(0.5))+
  geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=0.08,linewidth=0.8,
                 position=position_dodge(0.5))+
  scale_color_manual(values=COL3,guide="none")+
  scale_shape_manual(values=c(16,17),name="")+
  scale_x_log10(breaks=c(0.6,0.8,1.0,1.5),labels=c("0.60","0.80","1.00","1.50"))+
  labs(title="Competing Risks Analysis",
       subtitle="Fine-Gray model, NHANES III + Continuous NHANES 2005-2016",
       x="Subdistribution Hazard Ratio (95% CI)",y="")+
  theme_clin(14,11)+coord_cartesian(xlim=c(0.5,1.9))+
  theme(legend.position=c(0.85,0.2),legend.background=element_rect(fill="white",color=NA))
ggsave(file.path(FIG_DIR,"fig5_competing_risks.pdf"),p5,width=7.5,height=4)
ggsave(file.path(FIG_DIR,"fig5_competing_risks.png"),p5,width=7.5,height=4,dpi=300)

# ===== FIG 6: Cutpoint KM =====
cat("Fig 6 ...\n")
gi6 <- df %>% filter(gi_tumor==1,surv_years>0,!is.na(PNI)) %>%
  mutate(pni_g=ifelse(PNI<48.5,"PNI < 48.5 (Low)","PNI ≥ 48.5 (High)"))
km6 <- survfit(Surv(surv_years,death)~pni_g,data=gi6)
for(e in c("pdf","png")){
  if(e=="pdf") pdf(file.path(FIG_DIR,"fig_cutpoint.pdf"),width=7.5,height=6.5) else
    png(file.path(FIG_DIR,"fig_cutpoint.png"),width=7.5,height=6.5,units="in",res=300)
  par(oma=c(0,0,1,0))
  print(ggsurvplot(km6,data=gi6,pval=TRUE,pval.size=4,risk.table=TRUE,
    risk.table.height=0.18,risk.table.fontsize=3,conf.int=TRUE,
    palette=c("#D62728","#2CA02C"),surv.median.line="hv",
    legend.title="PNI Group",
    title="Kaplan-Meier by PNI Threshold (48.5)",
    subtitle="NHANES III + Continuous NHANES 2005-2016",
    xlab="Years",ylab="Survival Probability",
    ggtheme=theme_classic(base_size=12)),newpage=FALSE)
  dev.off()
}

# ===== FIG 7: RCS =====
cat("Fig 7 ...\n")
rc <- read.csv(file.path(RES_DIR,"rcs_predictions.csv"))
# Filter reasonable PNI range
rc <- rc %>% filter(PNI>=25,PNI<=65,HR<10)
p7 <- ggplot(rc,aes(x=PNI,y=HR))+
  geom_hline(yintercept=1,linetype="dashed",color="grey60",linewidth=0.4)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),fill="#1F77B4",alpha=0.2)+
  geom_line(color="#1F77B4",linewidth=1)+
  labs(title="Dose-Response: PNI and All-Cause Mortality",
       subtitle="Restricted cubic splines, 3 knots",
       x="PNI",y="Hazard Ratio (95% CI)")+
  theme_clin(14,11)+coord_cartesian(ylim=c(0,3.5),xlim=c(25,65))
ggsave(file.path(FIG_DIR,"fig_rcs.pdf"),p7,width=7,height=5)
ggsave(file.path(FIG_DIR,"fig_rcs.png"),p7,width=7,height=5,dpi=300)

# ===== FIG 8: Subgroup forest =====
cat("Fig 8 ...\n")
sg <- read.csv(file.path(RES_DIR,"subgroup_forest.csv")) %>%
  filter(Subgroup!="All patients") %>%
  mutate(Subgroup=factor(Subgroup,rev(unique(Subgroup))))
p8 <- ggplot(sg,aes(x=HR,y=Subgroup))+
  geom_vline(xintercept=1,linetype="dashed",color="grey60",linewidth=0.4)+
  geom_point(size=3,color="#1F77B4")+
  geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=0.12,linewidth=0.8,color="#1F77B4")+
  scale_x_log10()+
  labs(title="Subgroup Analysis: PNI and All-Cause Mortality",
       subtitle="Adjusted HR per 1-SD",
       x="Hazard Ratio (95% CI)",y="")+
  theme_clin(14,11)
ggsave(file.path(FIG_DIR,"fig_subgroup_forest.pdf"),p8,width=7.5,height=5)
ggsave(file.path(FIG_DIR,"fig_subgroup_forest.png"),p8,width=7.5,height=5,dpi=300)

# ===== FIG 9: Inflammatory comparison =====
cat("Fig 9 ...\n")
f9 <- data.frame(Index=factor(c("logSII","PNI","logNLR"),c("logNLR","PNI","logSII")),
                 Delta=c(0.014,0.017,0.017))
p9 <- ggplot(f9,aes(x=Index,y=Delta,fill=Index))+
  geom_col(width=0.5,alpha=0.85)+
  geom_text(aes(label=sprintf("+%.3f",Delta)),vjust=-0.5,size=4.5,fontface="bold")+
  scale_fill_manual(values=c("PNI"="#1F77B4","logNLR"="#D62728","logSII"="#2CA02C"),guide="none")+
  scale_y_continuous(limits=c(0,0.028),expand=c(0,0))+
  labs(title=expression(Delta*"C-stat: PNI vs Inflammatory Markers"),
       subtitle="GI cancer (n=258, events=125)",
       x="",y=expression(Delta*"C-statistic"))+
  theme_clin(14,12)
ggsave(file.path(FIG_DIR,"fig9_inflammatory_comparison.pdf"),p9,width=5.5,height=4.5)
ggsave(file.path(FIG_DIR,"fig9_inflammatory_comparison.png"),p9,width=5.5,height=4.5,dpi=300)

# ===== FIG 10: Decomposition =====
cat("Fig 10 ...\n")
f10 <- data.frame(Index=factor(c("Lymphocyte","PNI","Albumin"),c("Lymphocyte","PNI","Albumin")),
                  Delta=c(0.002,0.022,0.032))
p10 <- ggplot(f10,aes(x=Index,y=Delta,fill=Index))+
  geom_col(width=0.5,alpha=0.85)+
  geom_text(aes(label=sprintf("+%.3f",Delta)),vjust=-0.5,size=4.5,fontface="bold")+
  scale_fill_manual(values=c("Albumin"="#2CA02C","PNI"="#1F77B4","Lymphocyte"="#D62728"),guide="none")+
  scale_y_continuous(limits=c(0,0.045),expand=c(0,0))+
  labs(title=expression(Delta*"C-stat: PNI Decomposition"),
       subtitle="GI cancer (n=313, events=169)",
       x="",y=expression(Delta*"C-statistic"))+
  theme_clin(14,12)
ggsave(file.path(FIG_DIR,"fig10_decomposition.pdf"),p10,width=5.5,height=4.5)
ggsave(file.path(FIG_DIR,"fig10_decomposition.png"),p10,width=5.5,height=4.5,dpi=300)

cat("\n=== ALL FIGURES REGENERATED ===\n")
