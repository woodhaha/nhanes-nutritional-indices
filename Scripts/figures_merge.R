# -- figures_merge.R
# Merge 10 figures → 6 main panels + 3 supplementary panels
# Layout:
#   Fig 1: KM PNI tertile (top) | Cutpoint 48.5 (bottom)  — taller panel
#   Fig 2: Cox forest (left) | Competing risks sHR (right) — side-by-side
#   Fig 3: Time-dependent HR (top) | Landmark (bottom) — stacked
#   Fig 4: RCS dose-response (standalone)
#   Fig 5: Inflam ΔC-stat (left) | Decomposition ΔC-stat (right) — side-by-side
#   Fig 6: Subgroup forest (standalone)
#
# Supp: S1 CALLY KM | S2 Joint PNI+NLR | S3 RMST

library(ggplot2); library(dplyr); library(survival); library(survminer); library(patchwork)

FIG_DIR <- normalizePath("figures/gi_analysis")
RES_DIR <- normalizePath("results/gi_analysis")

# -- Theme --
theme_clin <- function(base=12, axis=10) {
  theme_classic(base_size=base) +
    theme(plot.title=element_text(face="bold",size=base,hjust=0.5,margin=margin(b=4)),
          plot.subtitle=element_text(size=base-2,color="grey40",hjust=0.5,margin=margin(b=4)),
          axis.title=element_text(size=base-1),
          axis.text=element_text(size=axis),
          legend.title=element_text(size=base-2),
          legend.text=element_text(size=base-3),
          plot.margin=margin(4,4,4,4),
          panel.spacing=unit(4,"pt"))
}
COL3 <- c("PNI"="#1F77B4","CONUT"="#2CA02C","GNRI"="#FF7F0E")

df <- readRDS(file.path(RES_DIR,"nhanes_clean.rds"))

# ======================================================================
# FIGURE 1: KM PNI Tertile + Cutpoint (two-panel, row stack)
# ======================================================================
cat("Fig 1: KM panels ...\n")
gi <- df %>% filter(gi_tumor==1,surv_years>0,!is.na(PNI)) %>%
  mutate(pni_t=factor(ntile(PNI,3),1:3,c("Low PNI (T1)","Mid PNI (T2)","High PNI (T3)")))

# Panel A: KM by tertile (no risk table, saved as grob)
km1 <- survfit(Surv(surv_years,death)~pni_t,data=gi)
g1a <- ggsurvplot(km1,data=gi,pval=TRUE,pval.size=3,conf.int=TRUE,
  palette=unname(COL3[c(2,1,3)]),surv.median.line="hv",
  legend.title="PNI Tertile",legend="none",
  title="All-Cause Survival by PNI Tertile",
  xlab="",ylab="Survival Probability",
  ggtheme=theme_classic(base_size=10))
p1a <- g1a$plot + theme(plot.margin=margin(4,4,0,4))

# Panel B: Cutpoint KM (no risk table)
gi6 <- gi %>% mutate(pni_g=ifelse(PNI<48.5,"PNI < 48.5 (Low)","PNI >= 48.5 (High)"))
km6 <- survfit(Surv(surv_years,death)~pni_g,data=gi6)
g1b <- ggsurvplot(km6,data=gi6,pval=TRUE,pval.size=3,conf.int=TRUE,
  palette=c("#D62728","#2CA02C"),surv.median.line="hv",
  legend.title="PNI Group",legend="none",
  title="Kaplan-Meier by PNI Threshold (48.5)",
  xlab="Years",ylab="Survival Probability",
  ggtheme=theme_classic(base_size=10))
p1b <- g1b$plot + theme(plot.margin=margin(0,4,4,4))

fig1 <- (p1a / p1b) + plot_annotation(tag_levels="A",
  theme=theme(plot.margin=margin(4,4,4,4)))
ggsave(file.path(FIG_DIR,"fig1_main_km.pdf"),fig1,width=7.5,height=7)
ggsave(file.path(FIG_DIR,"fig1_main_km.png"),fig1,width=7.5,height=7,dpi=300)

# ======================================================================
# FIGURE 2: Cox forest (left) + Competing risks (right)
# ======================================================================
cat("Fig 2: Forest + Competing risks ...\n")

# Panel A: Cox forest
cr <- read.csv(file.path(RES_DIR,"cox_main.csv")) %>% filter(Adj=="Adjusted") %>%
  mutate(Index=factor(Index,c("GNRI","CONUT","PNI")))
p2a <- ggplot(cr,aes(x=HR,y=Index,color=Index)) +
  geom_vline(xintercept=1,linetype="dashed",color="grey70",linewidth=0.3)+
  geom_point(size=2.5)+
  geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=0.06,linewidth=0.7)+
  scale_color_manual(values=COL3,guide="none")+
  scale_x_log10(breaks=c(0.6,0.8,1.0),labels=c("0.60","0.80","1.00"))+
  labs(title="All-cause mortality",x="HR per 1-SD (95% CI)",y="")+
  theme_clin(11,9)+coord_cartesian(xlim=c(0.5,1.05))+
  theme(plot.title=element_text(size=10,face="bold",hjust=0.5))

# Panel B: Competing risks sHR
fg <- read.csv(file.path(RES_DIR,"competingrisks_finegray.csv")) %>%
  mutate(Risk=ifelse(grepl("Non",Risk),"Non-cancer","Cancer"),
         Index=factor(Index,c("GNRI","CONUT","PNI")))
p2b <- ggplot(fg,aes(x=sHR,y=Index,color=Index,shape=Risk))+
  geom_vline(xintercept=1,linetype="dashed",color="grey70",linewidth=0.3)+
  geom_point(size=2.5,position=position_dodge(0.4))+
  geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=0.06,linewidth=0.7,
                 position=position_dodge(0.4))+
  scale_color_manual(values=COL3,guide="none")+
  scale_shape_manual(values=c(16,17),name="")+
  scale_x_log10(breaks=c(0.6,0.8,1.0,1.5),labels=c("0.60","0.80","1.00","1.50"))+
  labs(title="Competing risks (Fine-Gray)",x="sHR (95% CI)",y="")+
  theme_clin(11,9)+coord_cartesian(xlim=c(0.5,1.9))+
  theme(plot.title=element_text(size=10,face="bold",hjust=0.5),
        legend.position=c(0.75,0.15),legend.background=element_blank(),
        legend.text=element_text(size=7))

fig2 <- (p2a | p2b) + plot_annotation(tag_levels="A",
  theme=theme(plot.margin=margin(4,4,4,4)))
ggsave(file.path(FIG_DIR,"fig2_main_forest.pdf"),fig2,width=7.5,height=3.5)
ggsave(file.path(FIG_DIR,"fig2_main_forest.png"),fig2,width=7.5,height=3.5,dpi=300)

# ======================================================================
# FIGURE 3: Time-dependent HR (top) + Landmark (bottom)
# ======================================================================
cat("Fig 3: Time-dependent + Landmark ...\n")

# Panel A: Time-dependent
gi3 <- df %>% filter(gi_tumor==1,surv_years>0,!is.na(PNI)) %>%
  mutate(PNI_s=scale(PNI)[,1],CONUT_s=scale(-CONUT)[,1],GNRI_s=scale(GNRI)[,1])
td3 <- bind_rows(lapply(c("PNI_s","CONUT_s","GNRI_s"),function(idx){
  d <- survSplit(Surv(surv_years,death)~.,data=gi3,cut=c(2,5),episode="ep")
  d$ep_f <- factor(d$ep,0:2)
  f <- coxph(as.formula(paste0("Surv(tstart,surv_years,death)~",idx,"*ep_f+age+sex+race_eth")),data=d)
  s <- summary(f)$coefficients
  r <- grep(gsub("_s","",idx),rownames(s))
  b <- s[r,1]; se <- s[r,3]
  h <- c(b[1],b[1]+b[2],b[1]+b[3])
  hs <- c(se[1],sqrt(se[1]^2+se[2]^2),sqrt(se[1]^2+se[3]^2))
  data.frame(Index=gsub("_s","",idx),Period=factor(c("0-2 y","2-5 y","5+ y"),c("0-2 y","2-5 y","5+ y")),
             HR=exp(h),Lower=exp(h-1.96*hs),Upper=exp(h+1.96*hs))
})); td3$Index <- factor(td3$Index,c("PNI","CONUT","GNRI"))
p3a <- ggplot(td3,aes(x=Index,y=HR,color=Index,shape=Period))+
  geom_hline(yintercept=1,linetype="dashed",color="grey70",linewidth=0.3)+
  geom_point(position=position_dodge(0.6),size=2.5)+
  geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.1,linewidth=0.6,
                position=position_dodge(0.6))+
  scale_color_manual(values=COL3,guide="none")+
  scale_shape_manual(values=c(16,17,15),name="")+
  labs(title="Time-dependent hazard ratios",x="",y="HR per 1-SD (95% CI)")+
  theme_clin(11,9)+coord_cartesian(ylim=c(0,2.8))+
  theme(plot.title=element_text(size=10,face="bold",hjust=0.5),
        legend.position="bottom",legend.text=element_text(size=8))

# Panel B: Landmark
lm <- read.csv(file.path(RES_DIR,"cox_landmark.csv")) %>%
  mutate(Landmark_clean=case_when(grepl("1",Landmark)~"1 yr",grepl("3",Landmark)~"3 yr",TRUE~"5 yr"),
         Index=factor(Index,c("GNRI","CONUT","PNI")))
p3b <- ggplot(lm,aes(x=Landmark_clean,y=HR,color=Index,group=Index))+
  geom_hline(yintercept=1,linetype="dashed",color="grey70",linewidth=0.3)+
  geom_point(size=2.5,position=position_dodge(0.4))+
  geom_errorbar(aes(ymin=Lower,ymax=Upper),width=0.05,linewidth=0.6,
                position=position_dodge(0.4))+
  geom_line(position=position_dodge(0.4),linewidth=0.4)+
  scale_color_manual(values=COL3,name="")+
  labs(title="Landmark analysis (conditional survival)",x="",y="HR per 1-SD (95% CI)")+
  theme_clin(11,9)+coord_cartesian(ylim=c(0.55,1.2))+
  theme(plot.title=element_text(size=10,face="bold",hjust=0.5),
        legend.position="bottom",legend.text=element_text(size=8))

fig3 <- (p3a / p3b) + plot_annotation(tag_levels="A",
  theme=theme(plot.margin=margin(4,4,4,4))) + plot_layout(heights=c(1.2,1))
ggsave(file.path(FIG_DIR,"fig3_main_timedep.pdf"),fig3,width=7,height=7)
ggsave(file.path(FIG_DIR,"fig3_main_timedep.png"),fig3,width=7,height=7,dpi=300)

# ======================================================================
# FIGURE 4: RCS dose-response (standalone)
# ======================================================================
cat("Fig 4: RCS ...\n")
rc <- read.csv(file.path(RES_DIR,"rcs_predictions.csv")) %>% filter(PNI>=25,PNI<=65,HR<10)
p4 <- ggplot(rc,aes(x=PNI,y=HR))+
  geom_hline(yintercept=1,linetype="dashed",color="grey70",linewidth=0.4)+
  geom_ribbon(aes(ymin=Lower,ymax=Upper),fill="#1F77B4",alpha=0.2)+
  geom_line(color="#1F77B4",linewidth=1)+
  labs(title="Dose-Response: PNI and All-Cause Mortality",
       subtitle="Restricted cubic splines (3 knots)",
       x="Prognostic Nutritional Index (PNI)",y="Hazard Ratio (95% CI)")+
  theme_clin(12,10)+coord_cartesian(ylim=c(0,4),xlim=c(25,65))+
  theme(plot.margin=margin(8,8,8,8))
ggsave(file.path(FIG_DIR,"fig4_main_rcs.pdf"),p4,width=6,height=5)
ggsave(file.path(FIG_DIR,"fig4_main_rcs.png"),p4,width=6,height=5,dpi=300)

# ======================================================================
# FIGURE 5: Inflam ΔC-stat (left) + Decomposition ΔC-stat (right)
# ======================================================================
cat("Fig 5: C-stat panels ...\n")

# Panel A: Inflammatory comparison
f9 <- data.frame(Index=factor(c("logSII","PNI","logNLR"),c("logNLR","PNI","logSII")),
                 Delta=c(0.014,0.017,0.017))
p5a <- ggplot(f9,aes(x=Index,y=Delta,fill=Index))+
  geom_col(width=0.5,alpha=0.85)+
  geom_text(aes(label=sprintf("+%.3f",Delta)),vjust=-0.5,size=3.5,fontface="bold")+
  scale_fill_manual(values=c("PNI"="#1F77B4","logNLR"="#D62728","logSII"="#2CA02C"),guide="none")+
  scale_y_continuous(limits=c(0,0.028),expand=c(0,0))+
  labs(title=expression(Delta*"C-stat: PNI vs Inflammation"),
       subtitle="GI cancer (n=258, events=125)",
       x="",y=expression(Delta*"C-statistic"))+
  theme_clin(10,9)

# Panel B: Decomposition
f10 <- data.frame(Index=factor(c("Lymphocyte","PNI","Albumin"),c("Lymphocyte","PNI","Albumin")),
                  Delta=c(0.002,0.022,0.032))
p5b <- ggplot(f10,aes(x=Index,y=Delta,fill=Index))+
  geom_col(width=0.5,alpha=0.85)+
  geom_text(aes(label=sprintf("+%.3f",Delta)),vjust=-0.5,size=3.5,fontface="bold")+
  scale_fill_manual(values=c("Albumin"="#2CA02C","PNI"="#1F77B4","Lymphocyte"="#D62728"),guide="none")+
  scale_y_continuous(limits=c(0,0.045),expand=c(0,0))+
  labs(title=expression(Delta*"C-stat: PNI Decomposition"),
       subtitle="GI cancer (n=313, events=169)",
       x="",y=expression(Delta*"C-statistic"))+
  theme_clin(10,9)

fig5 <- (p5a | p5b) + plot_annotation(tag_levels="A",
  theme=theme(plot.margin=margin(4,4,4,4)))
ggsave(file.path(FIG_DIR,"fig5_main_cstat.pdf"),fig5,width=7.5,height=3.5)
ggsave(file.path(FIG_DIR,"fig5_main_cstat.png"),fig5,width=7.5,height=3.5,dpi=300)

# ======================================================================
# FIGURE 6: Subgroup forest (standalone)
# ======================================================================
cat("Fig 6: Subgroup forest ...\n")
sg <- read.csv(file.path(RES_DIR,"subgroup_forest.csv")) %>%
  filter(Subgroup!="All patients") %>%
  mutate(Subgroup=factor(Subgroup,rev(unique(Subgroup))))
p6 <- ggplot(sg,aes(x=HR,y=Subgroup))+
  geom_vline(xintercept=1,linetype="dashed",color="grey70",linewidth=0.4)+
  geom_point(size=3,color="#1F77B4")+
  geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=0.1,linewidth=0.7,color="#1F77B4")+
  scale_x_log10()+
  labs(title="Subgroup Analysis: PNI and All-Cause Mortality",
       subtitle="Adjusted HR per 1-SD",
       x="Hazard Ratio (95% CI)",y="")+
  theme_clin(12,10)+theme(plot.margin=margin(8,8,8,8))
ggsave(file.path(FIG_DIR,"fig6_main_subgroup.pdf"),p6,width=7,height=4)
ggsave(file.path(FIG_DIR,"fig6_main_subgroup.png"),p6,width=7,height=4,dpi=300)

cat("\n=== ALL MERGED FIGURES DONE ===\n")
