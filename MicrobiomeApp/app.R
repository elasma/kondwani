#setwd("~/Documents/Consult/Milk/Analysis/Microbiome/ubiome2")
source("Packages.R")
#source("libs.R")
source("ANCOM_updated_code.R")
source("permanova.R")

ui <- fluidPage(
  
    # Application title
    titlePanel("Key Results for Microbiome data"),

        mainPanel(  "Upload Data",
                      fileInput("phyObj", "Upload Phyloseq object as .rds or RData"),
                      uiOutput("varselect"),
                      uiOutput("grpselect"),
            div(navlistPanel("Overview",
                        tabPanel("Relative Abundance", 
                                 plotOutput("abPlot",width="600px"),selectInput("var", 
                        label = "Choose taxanomic level",
                        choices = c("Phylum", 
                        "Order",
                        "Class", 
                        "Family",
                        "Genus"),
                        selected = "Phylum"),
                        radioButtons("gp", label = "Choose graph type",
                                    choices = list("Box plot" = "geom_box",
                                                   "Bar plot" = "geom_bar"),
                                              selected = "geom_box"),
                        conditionalPanel(
                                     condition = "input.gp == 'geom_bar'",
                                     radioButtons("stat", 
                                        label="Summary", 
                                        choices = list("Median" = "median", 
                                                        "Mean" = "mean"), 
                                        selected = "median"))),
                        "Ordination",
                tabPanel("Pattern Exploration", plotOutput("pattPlot",width="600px"),
                        selectInput("ordi",label ="Ordination method",
                        choices = c("DCA", "CCA", "RDA", "MDS","NMDS", "PCoA"),
                        selected = "PCoA"),
                helpText("Select a grouping variable under 'Groupby'")),
                        #uiOutput("grpselect")),                
                  "Alpha Diversity", 
                       
             tabPanel("Summary Statistics",
                      #uiOutput("grpselect"),
                      helpText("Select a grouping variable under 'Groupby'"),
                      tableOutput("adivTable"),
                      tableOutput("adivOver")),
             tabPanel("Visual Exploration",
                      plotOutput("adivPlot",width = "900px"),
                      #uiOutput("grpselect"),
                      helpText("Select a grouping variable under 'Groupby'")),
             tabPanel("Association Analysis",
                      selectInput("divm","Select diversity measure outcome",
                              choices = list(
                                "Shannon Index"="Shannon",
                                "Chao"="Chao1",
                                "Inverse Simpson"="InvSimpson"
                              )),
                      checkboxInput("alphaout","Treat diversity measure as a covariate"),
                      conditionalPanel(
                          condition = "input.alphaout==true",
                          #uiOutput("grpselect"),
                          helpText("Select an outcome variable under 'Groupby'"),
                          selectInput("outype",label ="Outcome type",
                          choices = list("Continuous"="gaussian","Binary"="binomial",
                                          "Poisson"="poisson","Categorical"="mult"),
                          selected = "binomial")),
                      helpText("Choose covariates to be adjusted",
                               "for from dataset variables"),
                      tableOutput("regress")),
             "Beta diversity",
             tabPanel("PERMANOVA",tableOutput("perm"),
                      #uiOutput("grpselect"),
                      helpText("Select a grouping variable under 'Groupby'")),
             tabPanel("Dirichlet Model",
                      #uiOutput("grpselect"),
                      helpText("Select a grouping variable under 'Groupby'"),
                      tableOutput("DM")),
             "Differential Taxa",
             tabPanel("Kruskal Wallis",  radioButtons("tax", 
                                  label = "Choose taxanomic level",
                                  choices = c("Phylum", 
                                              "Order",
                                              "Class", 
                                              "Family",
                                              "Genus",
                                              "Species"),
                                  selected = "Genus"),
                      radioButtons("dis",label = "Display results as",
                                   choices = list("Table"="tab",
                                                  "Graphs"="gph")),
                      selectInput("padj", 
                                  label = "Adjust pvalue?",
                                  choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                     "fdr", "none"),
                                  selected = "fdr"),
                      helpText("Choose method to adjust pvalue for multiple comparisons"),
                      helpText("Select a grouping variable under 'Groupby'"),
                      sliderInput("num","choose value of alpha", value = 0.01,
                                   min=0,max=1,step = 0.001),
                      tableOutput("kwt"),
                      plotOutput("kwp",width="600px")
                      ),
             # tabPanel("Wilcoxon rank",  selectInput("tax", 
             #         label = "Choose taxanomic level",
             #                   choices = c("Phylum", 
             #                               "Order",
             #                               "Class", 
             #                               "Family",
             #                               "Genus",
             #                               "Species"),
             #                   selected = "Genus"),
             #         selectInput("dispy",label = "Display results as",
             #               choices = list("Table"="tab",
             #                              "Graphs"="gph")),
             #         #uiOutput("grpselect"),
             #  helpText("Select a grouping variable under 'Groupby'"),
             #  sliderInput("num","choose value of alpha", value = 0.1,
             #               min=0,max=1,step = 0.001),
             #  tableOutput("kwt.wx"),
             #  plotOutput("kwp.wx",width="600px")
             # ),
         tabPanel("Negative Binomial",  radioButtons("tax", 
                    label = "Choose taxanomic level",
                    choices = c("Phylum", 
                                "Order",
                                "Class", 
                                "Family",
                                "Genus",
                                "Species"),
                    selected = "Genus"),
                  radioButtons("dispy",label = "Display results as",
                              choices = list("Table"="tab",
                                             "Graphs"="gph")),
                  #uiOutput("grpselect"),
                  helpText("Select a grouping variable under 'Groupby'"),
                  sliderInput("num","choose value of alpha", value = 0.1,
                               min=0,max=1,step=0.001),
                  tableOutput("nb.tb"),
                  plotOutput("nb.pt",width="600px")
         ),
         tabPanel("ANCOM", 
                  radioButtons("tax", 
                                    label = "Choose taxanomic level",
                                    choices = c("Phylum", 
                                                "Order",
                                                "Class", 
                                                "Family",
                                                "Genus",
                                                "Species"),
                                    selected = "Genus"),
                  radioButtons("dispy",label = "Display results as",
                              choices = list("Table"="tab",
                                             "Graphs"="gph")),
                  #uiOutput("grpselect"),
                  helpText("Select a grouping variable under 'Groupby'"),
                  sliderInput("ancAlpha","choose value of alpha", value = 0.1,
                              min=0,max=1,step=0.001),
                  tableOutput("anc.tb"),
                  plotOutput("anc.pt",width="600px")
         ) ,id="nav")
             
                       ))
    
)

# Define server 
server <- function(input, output,session) {
  session$onSessionEnded(stopApp)
  Phy=reactive({
    req(input$phyObj)
    readRDS(input$phyObj$datapath)
  })

  phy.list=reactive({
    phy.1=Phy()
    ln=length(rank_names(Phy()))
    tax.rank=rank_names(Phy())
    tax.rank=tax.rank[ln]
    physeq2=transform_sample_counts(phy.1, function(x) x / sum(x) )
    phy.2.R=tax_glom(physeq2, taxrank =tax.rank)
    phy.2=tax_glom(phy.1, taxrank =tax.rank)
    phy.meta=meta(phy.2)
    nm.meta=names(phy.meta)
    phy.list=list(physeq2=physeq2,phy.2.R=phy.2.R,phy.2=phy.2,
                  phy.meta=phy.meta,nm.meta=nm.meta)
    phy.list
  })
  output$varselect=renderUI({
    varSelectInput("datavars", "Dataset variables", data=phy.list()[[4]],
                   multiple=TRUE,selected =NULL )
  })
  output$grpselect=renderUI({
    selectInput("tblvar", "Groupby", choices=phy.list()[[5]],
                selected = NULL)
  })
  
    output$abPlot <- renderPlot({
        if(input$gp=="geom_bar")
        {
        ff1=function(cc,st){
            cc1=enquo(cc)
            st1=paste(st)
            physeq3 <- tax_glom(phy.list()[[2]], taxrank=cc)
            dat.phy=data.table(psmelt(physeq3))
            # create dataframe from phyloseq object
            dat.phy1=dat.phy %>% select(Abundance,noquote(eval(cc))) %>% 
                group_by_(eval(cc)) %>% filter(Abundance>0) %>%
                summarise_at("Abundance",eval(substitute(st1))) 
            dat.phy1=setorder(dat.phy1,-Abundance)
            dat.phy1=dat.phy1[1:5,]
            ggplot(dat.phy1,aes(x=get(cc),y=Abundance)) +
                geom_bar(stat = "identity",position="dodge")  + 
                labs(x=paste(cc), y=paste(st1, "relative abundance"))
        }
        ff1(input$var,input$stat)
        }
else 
{physeq3 <- tax_glom(phy.list()[[2]], taxrank =input$var)
    # create dataframe from phyloseq object
    dat.phy=data.table(psmelt(physeq3))
    #calculate median rel. abundance
    dat.phy[, median := median(Abundance, na.rm = TRUE), 
            by =eval(input$var)]
    ff=function(cc){
        ggplot(dat.phy[median >=0.01],
               aes(x=get(cc),
                   y=Abundance)) + geom_boxplot()  +  
            labs(x=paste(input$var), y="Relative abundance")  +
            coord_flip()
    }
    ff(input$var)}
        
    })
    output$pattPlot <- renderPlot({
        ordi = ordinate(phy.list()[[2]], method=input$ordi, distance="bray")
        pp=function(cc){
        plot_ordination(phy.list()[[2]], ordi, "samples", color=cc) 
        }
       pp(input$tblvar)

    })
    div_df=reactive({
        div=estimate_richness(phy.list()[[3]],measures=c("Observed", "Chao1","Shannon","InvSimpson"))
        meta_div=cbind(meta(phy.list()[[3]]),div)
        meta_div=as.data.frame(meta_div)
        meta_div
    })
   #Summary tables     
        adivtables=reactive({
   t1= div_df() %>% select(c("Chao1","Shannon","InvSimpson", input$tblvar)) %>% group_by_at(input$tblvar) %>% 
    summarise_at("Chao1",c("median","mean","min","max")) %>% 
            mutate(Measure="Chao") %>% select(Measure,everything())
   t2= div_df() %>% select(c("Chao1","Shannon","InvSimpson", input$tblvar)) %>% group_by_at(input$tblvar) %>% 
       summarise_at("Shannon",c("median","mean","min","max")) %>% 
       mutate(Measure="Shannon Index") %>% select(Measure,everything())
   t3= div_df() %>% select(c("Chao1","Shannon","InvSimpson", input$tblvar)) %>% group_by_at(input$tblvar) %>% 
       summarise_at("InvSimpson",c("median","mean","min","max")) %>% 
       mutate(Measure="Inverse Simpson") %>% select(Measure,everything())
   ttt=bind_rows(t1,t2,t3)
   return(ttt)
        })
    adivOverR=reactive({
        t4= div_df() %>% select(c("Chao1","Shannon","InvSimpson"))%>% 
            summarise_all(c("median","mean","min","max"))
        
        time.vars=list(c("Chao1_median", "Shannon_median", "InvSimpson_median"),
                       c("Chao1_mean", "Shannon_mean", "InvSimpson_mean"),
                       c("Chao1_min", "Shannon_min", "InvSimpson_min"),
                       c("Chao1_max", "Shannon_max", "InvSimpson_max"))
        long.vars=c("median","mean","min","max")
        t4long=reshape(as.data.frame(t4), varying = time.vars,v.names =long.vars,
                       times=c("Chao","Shannon","InvSimpson"),direction ="long")
        t4long$Measure=t4long$time
        t4long$ungrouped="overall"
        t4l=t4long %>% select(Measure,ungrouped,median, mean,min,max)
        rownames(t4l)=NULL
        t4l     
    })
    output$adivOver=renderTable({
      adivOverR()
    })
    
   output$adivTable=renderTable({ 
       adivtables()
   })
#Summary plots
 {
   output$adivPlot=renderPlot({
     divgrp=div_df()%>% select(SampleID,input$tblvar)
     divgrp$grpvar=as.factor( divgrp[,2])
     meta_div=div_df()
     meta_div$grpvar=as.character(divgrp$grpvar)
      
       ch1=ggplot( meta_div,aes(x=grpvar,y=Chao1)) + geom_boxplot()  +  
           labs(x="",y="Chao") 
       sh1=ggplot( meta_div,aes(x=grpvar,y=Shannon)) + geom_boxplot()  +  
           labs(x="", y="Shannon") 
       ch2=   ggplot( meta_div,aes(x=grpvar,y=InvSimpson)) + geom_boxplot()  +  
           labs(x="", y="Inverse Simpson")  
       ggarrange(ch1, sh1,ch2 ,nrow = 1, ncol = 3)  
       

   })
  }
    
    reg.res=reactive({
        
        if(input$alphaout==FALSE){
            x=input$datavars
            if (is.null(x))
                x <- character(0)
            shdf=div_df()%>% select(input$divm,!!!input$datavars)
            shdf$out=shdf[,1]
            shdf=shdf%>% select(out,!!!input$datavars)
            shreg=lm(out~ . ,data=shdf)
            summ_sh=summary(shreg)
            coef_sh=summ_sh$coefficients
          as.data.table(coef_sh,keep.rownames=TRUE)
        }
        else{
        if (input$outype=="mult"){
          x=input$datavars
          if (is.null(x))
            x <- character(0)
          shdf=div_df()%>% select(input$tblvar,input$divm,!!!input$datavars)
          shdf$out=as.factor(shdf[,1])
          if (length(levels(shdf$out))<3){stop("use the binary option")}
          else{
          n=length(levels(as.factor(shdf$out)))-1
          shdf=shdf%>% select(out,input$divm,!!!input$datavars)
          mod1 <- multinom(out~ ., data = shdf)
          
          #Obtain zscores
          summ=summary(mod1)
          cf=summ$coefficients
          cz=colnames(cf)
          czn=length(cz)
          rz=rownames(cf)
          rzn=length(rz)
          cz=rep(cz,each=rzn)
          rz=rep(rz,czn)
          parameter=paste(cz,rz,sep=":")
          cf=round(as.vector(summ$coefficients),3)
          z <- round(as.vector(summ$coefficients/summ$standard.errors),3)
          #obtain Pvalues
          p <- round(as.vector((1 - pnorm(abs(z), 0, 1)) * 2),3)
          se=round(as.vector(summ$standard.errors),3)
          all_res=cbind( parameter,cf,se,z,p) 
          colnames(all_res)=c("Parameter", "Estimate","Std.error","Z","P-value")
          as.data.table(all_res)}}
            else{
              x=input$datavars
              if (is.null(x))
                x <- character(0)
              shdf=div_df()%>% select(input$tblvar,input$divm,!!!input$datavars)
              shdf$out=shdf[,1]
              shdf=shdf%>% select(out,input$divm,!!!input$datavars)
              shglm=glm(out~.,data=shdf,family =input$outype)
              summ_sh=summary(shglm)
              coef_sh=summ_sh$coefficients
              as.data.table(coef_sh,keep.rownames=TRUE)
            }
        }
    })
    output$regress=renderTable({
        reg.res()
    })
    #PERMANOVA
    perm.r=reactive({
      adf=phy.list()[[4]]%>% select(SampleID,input$tblvar)
      adf$grpvar=as.factor(adf[,2])
      adf$grpvar=as.character(adf$grpvar)
      phy.sub=phy.list()[[3]]
      sample_data(phy.sub)$grpvar=adf$grpvar
      physeq.m1=phyloseq::subset_samples(phy.sub, !is.na(grpvar))
      sub.meta=meta(physeq.m1)
      adf=sub.meta%>% select(SampleID,grpvar)
      adf$grpvar=as.factor(adf[,2])
      bray_dist_tp <- vegdist(otu_table(physeq.m1))
      if (length(levels(adf$grpvar))<3)
      {
        adm=adonis(bray_dist_tp ~ grpvar, data =adf)
        ad=summary(adm)
        ad=adm$aov.tab
        head(ad)
        as.data.table(cbind(ad[1],ad[2],ad[3],ad[4]),keep.rownames = TRUE)
      }
      else{
        adm=adonis(bray_dist_tp ~ grpvar, data =adf)
        ad=summary(adm)
        ad=adm$aov.tab
        head(ad)
        as.data.table(cbind(ad[1],ad[2],ad[3],ad[4]),keep.rownames = TRUE)
        pairwise.adonis(otu_table(physeq.m1), adf$grpvar, sim.method = "bray")}
    })
    output$perm=renderTable({
      perm.r()
    })
  DMr=reactive({
    #Dirichlet-multinomial Model
    adf=phy.list()[[4]]%>% select(SampleID,input$tblvar)
    adf$grpvar=as.factor(adf[,2])
    lv=levels(adf$grpvar)
    adf$grpvar=as.character(adf$grpvar)
    phy.sub=phy.list()[[3]]
    sample_data(phy.sub$grpvar)=adf$grpvar
    gg=function(i) {
      keepTaxa3 =adf[which(adf$grpvar==i),]
      miss_samples1=rownames(keepTaxa3)
     c.1= prune_samples(sample_names(phy.list()[[3]]) %in% miss_samples1,phy.list()[[3]])
     sub.otu=as(otu_table(c.1),"matrix")
     return(sub.otu)
    }
    lvv=as.list(lv)
    sub.phy.list=sapply(lvv, gg,USE.NAMES=TRUE, simplify = TRUE)
    cmb=combn(sub.phy.list, 2, FUN = Xdc.sevsample, simplify = TRUE)
    combs=combn(lv,2)
    bb=apply(combs,2,function(x) paste(x,collapse =" vs "))
    colnames(cmb)=bb
    rownames(cmb)=c("Test statistic","P-value")
    as.data.table(cmb,keep.rownames=TRUE)
  })
  
  output$DM=renderTable({
    DMr()
  })
 kwtr=reactive({
   phy.3.R <- tax_glom(phy, taxrank=input$tax)
   p3.meta= meta(phy.3.R)
   p3.taxn=colnames(tax_table(phy.3.R))
   indx=match(input$tax,p3.taxn)
   if (indx!=length(p3.taxn)){
     p3.tax=tax_table(phy.3.R)[,1:indx]
     p3.otu=otu_table(phy.3.R)
     phy.3.R=phyloseq(otu_table(p3.otu,taxa_are_rows =FALSE),sample_data(p3.meta),tax_table(p3.tax))
   }
   #print(meta(phy.3.R)[1:5,1:3])
   pmeta=meta(phy.3.R)
   adf=pmeta%>% select(SampleID,input$tblvar)
   adf$grpvar=as.factor(adf[,2])
   lv=levels(adf$grpvar)
   sample_data(phy.3.R)$grpvar=adf$grpvar  
   
   combs=combn(lv,2) 
   nlv=c(1:(ncol(combs)))
   pmeta=meta(phy.3.R)
   fwfn=function(i){
     keepTaxa3 =pmeta[which(pmeta$grpvar %in% combs[,i]),]
     miss_samples1=rownames(keepTaxa3)
     phy_sub= prune_samples(sample_names(phy.3.R) %in% miss_samples1,phy.3.R)
     otu.df=as(otu_table(phy_sub),"matrix")
     Sample.ID=rownames(otu.df)
     rownames(otu.df)=NULL
     Group=as.factor(sample_data(phy_sub)$grpvar)
     vardata=cbind(Sample.ID,Group)
     
     wx=function(i){pairwise.wilcox.test(i, vardata[,2], p.adjust.method = input$padj)}
     wx.list=apply(otu.df, 2,wx)
     wx.list2=lapply(wx.list, function(i){i[3]})
     wx.list3=lapply(wx.list2, function(i){i<input$num})
     wx.list3=as.matrix(wx.list3)
     wx.list4=unlist(wx.list2)
     names(wx.list4)=NULL
     wx.res=as.data.table(cbind(wx.list3,wx.list4),keep.rownames=TRUE)
     names(wx.res)=c("otu","detected","p.adjusted")
     wx.res.sub=wx.res%>%filter(detected==TRUE)%>%select(otu,p.adjusted)
     if (length(rownames(wx.res.sub))==0){
       anc2= cbind("No differential otus","NA",paste(combs[,i],collapse = " vs "))
       colnames(anc2)=c(input$tax,"p.adjusted","compare")
       anc2}
     else{
       colnames(ancW)=NULL
       ancW=as.vector(wx.res.sub$otu)
       detc.mat=as(tax_table(phy_sub)[ancW,], "matrix")
       otu.names=as.data.table(detc.mat,keep.rownames=TRUE)
     wx.res.sub$compare=paste(combs[,i],collapse = " vs ")
     wx.res.sub=cbind(otu.names,wx.res.sub)
     wx.res.sub=wx.res.sub%>%select(input$tax,p.adjusted,compare)
     wx.res.sub}
   }
 gg=do.call(rbind,lapply(as.list(nlv), fwfn))   
gg
    })

  kwpr=reactive({
    gg=kwtr()
    gen.vec=t(unique(gg[,1]))
    gen.vec=gen.vec[gen.vec!="No differential otus"]
    if (length(gen.vec)==0){paste("No differential otus to plot")}
    else{
    gen.vec=trimws(gen.vec)
    phy.3.R <- tax_glom(phy.list()[[2]], taxrank=input$tax)
    phy.3.meta=data.table(psmelt(phy.3.R))
    phy.3.meta$tax1=phy.3.meta%>% select(input$tax)
 
    phy.3.meta$gpvar=phy.3.meta%>% select(input$tblvar)
    phy.3.meta$gpvar=as.character(phy.3.meta$gpvar)
    phy.3.meta$tax1=trimws(phy.3.meta$tax1)
    
    phy_sub=dplyr::filter(phy.3.meta,tax1 %in% gen.vec)
  
    pp=phy_sub %>% select(Abundance,tax1,gpvar)%>% 
      group_by(tax1,gpvar) %>%
      summarise_at("Abundance",mean) 
    is.null(length(pp$tax1))

      ggplot(pp["gpvar"!="NA"],aes(x=tax1,y=Abundance,fill=gpvar)) +
      geom_col(position="dodge")  + 
      labs(x=paste(input$tax), y=" Mean Relative abundance") +
      facet_wrap(~tax1,scale="free")
    }
  })


  output$kwt=renderTable({
    if (input$dis=="tab")
    kwtr()
  })

  output$kwp=renderPlot({
    if (input$dis=="gph")
    kwpr()
  })

  # kwtr.wx=reactive({
  #   phy.3.R <- tax_glom(phy.list()[[2]], taxrank=input$tax)
  #   p3.meta= meta(phy.3.R)
  #   p3.taxn=colnames(tax_table(phy.3.R))
  #   indx=match(input$tax,p3.taxn)
  #   if (indx!=length(p3.taxn)){
  #     p3.tax=tax_table(phy.3.R)[,1:indx]
  #     p3.otu=otu_table(phy.3.R)
  #     phy.3.R=phyloseq(otu_table(p3.otu,taxa_are_rows =FALSE),sample_data(p3.meta),tax_table(p3.tax))
  #   }
  #   adf=p3.meta%>% select(SampleID,input$tblvar)
  #   adf$grpvar=as.factor(adf[,2])
  #   lv=levels(adf$grpvar)
  #   sample_data(phy.3.R)$grpvar=adf$grpvar  
  #   combs=combn(lv,2) 
  #   nlv=c(1:(ncol(combs)))
  #   fwfn=function(i){
  #     wx=function(i){pairwise.wilcox.test(i, adf[,3], p.adjust.method = "fdr")}
  #     wx.list=apply(otu.df, 2,wx)
  #     wx.list2=lapply(wx.list, function(i){i[3]})
  #     wx.list3=lapply(wx.list2, function(i){i<input$alpha})
  #     wx.list3=as.matrix(wx.list3)
  #     wx.list4=unlist(wx.list2)
  #     names(wx.list4)=NULL
  #     wx.res=as.data.table(cbind(wx.list3,wx.list4),keep.rownames=TRUE)
  #     names(wx.res)=c("otu","detected","p.adjusted")
  #     wx.res.sub=wx.res%>%filter(detected==TRUE)%>%select(otu,p.adjusted)
  #     
  #     kw=test_differential_abundance_Wilcoxon(
  #       phy.3.R, group = "grpvar", 
  #       compare=combs[,i],p.adjust.method = "fdr")  
  #     kwt=kw$table 
  #     kwt$compare=paste(combs[,i],collapse = " vs ")
  #     kwt=kwt%>%select(compare,Taxon,padj,input$tax)%>%filter(padj <input$num)
  #     
  #   }
  #   ggw=do.call(rbind,lapply(as.list(nlv), fwfn))   
  #   ggw
  # })
  # 
  # kwpr.wx=reactive({
  #   ggw=kwtr.wx()
  #   gen.vec=str_squish(unique(ggw[,4]))
  #   if (length(gen.vec)==0){
  #     dlg_message("No differential taxa")$res
  #   }
  #   else{
  #   phy.3.R <- tax_glom(phy.list()[[2]], taxrank=input$tax)
  #   phy.3.meta=data.table(psmelt(phy.3.R))
  #   phy.3.meta$tax1=phy.3.meta%>% select(input$tax)
  #   
  #   phy.3.meta$gpvar=phy.3.meta%>% select(input$tblvar)
  #   phy.3.meta$gpvar=as.character(phy.3.meta$gpvar)
  #   phy.3.meta$tax1=str_squish(phy.3.meta$tax1)
  #   
  #   phy_sub=dplyr::filter(phy.3.meta,tax1 %in% gen.vec)
  #   
  #   pp=phy_sub %>% select(Abundance,tax1,gpvar)%>% 
  #     group_by(tax1,gpvar) %>%
  #     summarise_at("Abundance",mean) 
  #   
  #   ggplot(pp["gpvar"!="NA"],aes(x=tax1,y=Abundance,fill=gpvar)) +
  #     geom_col(position="dodge")  + 
  #     labs(x=paste(input$tax), y=" Mean Relative abundance") +
  #     facet_wrap(~tax1,scale="free")
  # }})
  # 
  # output$kwt.wx=renderTable({
  #   if (input$dispy=="tab")
  #     kwtr.wx()
  # })
  # 
  # output$kwp.wx=renderPlot({
  #   if (input$dispy=="gph")
  #     kwpr.wx()
  # })

  nbtb.r=reactive({
    phy.3.R <- tax_glom(phy.list()[[3]], taxrank=input$tax)
    p3.meta= meta(phy.3.R)
    p3.taxn=colnames(tax_table(phy.3.R))
    indx=match(input$tax,p3.taxn)
    if (indx!=length(p3.taxn)){
      p3.tax=tax_table(phy.3.R)[,1:indx]
      p3.otu=otu_table(phy.3.R)
      phy.3.R=phyloseq(otu_table(p3.otu,taxa_are_rows =FALSE),sample_data(p3.meta),tax_table(p3.tax))
    }
    #print(meta(phy.3.R)[1:5,1:3])
    pmeta=meta(phy.3.R)
    adf=pmeta%>% select(SampleID,input$tblvar)
    adf$grpvar=as.factor(adf[,2])
    lv=levels(adf$grpvar)
    sample_data(phy.3.R)$grpvar=adf$grpvar  
 
    combs=combn(lv,2) 
    nlv=c(1:(ncol(combs)))
    pmeta=meta(phy.3.R)
    fwfn=function(i){
      keepTaxa3 =pmeta[which(pmeta$grpvar %in% combs[,i]),]
      miss_samples1=rownames(keepTaxa3)
      phy_sub= prune_samples(sample_names(phy.3.R) %in% miss_samples1,phy.3.R)
      nb12 = phyloseq_to_deseq2(phy_sub, ~ grpvar)
      nb12 = DESeq(nb12, test="Wald", fitType="parametric",sfType="poscounts")
      res12 = results(nb12, cooksCutoff = FALSE)
      res12$compare=paste(combs[,i],collapse = " vs ")
      res12=res12[which(res12$padj< input$num),]
      sigtab = cbind(as(res12, "data.frame"), as(tax_table(phy_sub)[rownames(res12), ], "matrix"))
      sigtab=sigtab%>%select(compare,log2FoldChange,padj,input$tax)
    }
    ggw=do.call(rbind,lapply(as.list(nlv), fwfn))   
    ggw
  })
  
  nbpt.r=reactive({
    ggw=nbtb.r()
    gen.vec=str_squish(unique(ggw[,4]))
    if (length(gen.vec)==0){
      dlg_message("No differential taxa")$res
    }
    else{
    phy.3.R <- tax_glom(phy.list()[[2]], taxrank=input$tax)
    phy.3.meta=data.table(psmelt(phy.3.R))
    phy.3.meta$tax1=phy.3.meta%>% select(input$tax)
    
    phy.3.meta$gpvar=phy.3.meta%>% select(input$tblvar)
    phy.3.meta$gpvar=as.character(phy.3.meta$gpvar)
    phy.3.meta$tax1=str_squish(phy.3.meta$tax1)
    
    phy_sub=dplyr::filter(phy.3.meta,tax1 %in% gen.vec)
    
    pp=phy_sub %>% select(Abundance,tax1,gpvar)%>% 
      group_by(tax1,gpvar) %>%
      summarise_at("Abundance",mean) 
    
    ggplot(pp["gpvar"!="NA"],aes(x=tax1,y=Abundance,fill=gpvar)) +
      geom_col(position="dodge")  + 
      labs(x=paste(input$tax), y=" Mean Relative abundance") +
      facet_wrap(~tax1,scale="free")
  }})
  
  output$nb.tb=renderTable({
    if (input$dispy=="tab")
      nbtb.r()
  })
  
  output$nb.pt=renderPlot({
    if (input$dispy=="gph")
      nbpt.r()
  })
  
  anctb.r=reactive({
    phy.3.R <- tax_glom(phy.list()[[3]], taxrank=input$tax)
    p3.meta= meta(phy.3.R)
    p3.taxn=colnames(tax_table(phy.3.R))
    indx=match(input$tax,p3.taxn)
    if (indx!=length(p3.taxn)){
      p3.tax=tax_table(phy.3.R)[,1:indx]
      p3.otu=otu_table(phy.3.R)
      phy.3.R=phyloseq(otu_table(p3.otu,taxa_are_rows =FALSE),sample_data(p3.meta),tax_table(p3.tax))
    }
    #print(meta(phy.3.R)[1:5,1:3])
    pmeta=meta(phy.3.R)
    adf=pmeta%>% select(SampleID,input$tblvar)
    adf$grpvar=as.factor(adf[,2])
    lv=levels(adf$grpvar)
    sample_data(phy.3.R)$grpvar=adf$grpvar  
    
    combs=combn(lv,2) 
    nlv=c(1:(ncol(combs)))
    pmeta=meta(phy.3.R)
    fwfn=function(i){
      keepTaxa3 =pmeta[which(pmeta$grpvar %in% combs[,i]),]
      miss_samples1=rownames(keepTaxa3)
      phy_sub= prune_samples(sample_names(phy.3.R) %in% miss_samples1,phy.3.R)
      otu.df=as(otu_table(phy_sub),"matrix")
      Sample.ID=rownames(otu.df)
      rownames(otu.df)=NULL
      otu.df=cbind(Sample.ID,otu.df)
      otu.df=as.data.frame(otu.df)
      Group=as.factor(sample_data(phy_sub)$grpvar)
      vardata=cbind(Sample.ID,Group)
      anc=ANCOM.main(OTUdat=otu.df,
                 Vardat=vardata,
                 adjusted=FALSE,
                 repeated=F,
                 main.var="Group",
                 adj.formula=NULL,
                 repeat.var=NULL,
                 longitudinal=FALSE,
                 random.formula=NULL,
                 multcorr=2,
                 sig=input$ancAlpha,
                 prev.cut=0.97)
      ancW=anc$W.taxa[,c(1,6)]%>%filter(detected_0.6==TRUE)%>% select(otu.names)
      if (length(rownames(ancW))==0){
        anc2= cbind("No differential otus",paste(combs[,i],collapse = " vs "))
        colnames(anc2)=c(input$tax,"compare")
        anc2}
      else{
      colnames(ancW)=NULL
      ancW=as.vector(t(ancW))
      detc.mat=as(tax_table(phy_sub)[ancW,], "matrix")
      detc.mat.df=as.data.table(detc.mat,keep.rownames=TRUE)
      detc.mat.df$compare=paste(combs[,i],collapse = " vs ")
      sigtab=detc.mat.df%>%select(input$tax,compare)
      sigtab
    }}
    ggw=do.call(rbind,lapply(as.list(nlv), fwfn))   
    ggw
 })
  
  ancpt.r=reactive({
    ggw=anctb.r()
    gen.vec=t(unique(ggw[,1]))
    gen.vec=gen.vec[gen.vec!="No differential otus"]
    if (length(gen.vec)==0){paste("No differential otus to plot")}
    else{
    phy.3.R <- tax_glom(phy.list()[[2]], taxrank=input$tax)
    phy.3.meta=data.table(psmelt(phy.3.R))
    phy.3.meta$tax1=phy.3.meta%>% select(input$tax)
    phy.3.meta$gpvar=phy.3.meta%>% select(input$tblvar)
    phy.3.meta$gpvar=as.character(phy.3.meta$gpvar)
    phy.3.meta$tax1=trimws(phy.3.meta$tax1)
    gen.vec=trimws(gen.vec)
    phy_sub=dplyr::filter(phy.3.meta,tax1 %in% gen.vec)
    pp=phy_sub %>% select(Abundance,tax1,gpvar)%>% 
      group_by(tax1,gpvar) %>%
      summarise_at("Abundance",mean) 
    ggplot(pp["gpvar"!="NA"],aes(x=tax1,y=Abundance,fill=gpvar)) +
      geom_col(position="dodge")  + 
      labs(x=paste(input$tax), y=" Mean Relative abundance") +
      facet_wrap(~tax1,scale="free")
    }
  })


  output$anc.tb=renderTable({
    if (input$dispy=="tab")
      anctb.r()
  })
  
  output$anc.pt=renderPlot({
    if (input$dispy=="gph")
      ancpt.r()
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

