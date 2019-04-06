library(shiny)
library(ape)
library(data.table)
library(adegenet)
library(ggplot2)
library(dplyr)
library(inbreedR)
library(BiocManager)
BiocManager::install(c("SNPRelate", "qvalue"))
options(shiny.maxRequestSize=10*1024^2) 


shinyServer(
  function(input, output) {
    
    observeEvent(input$SNPdata, {
    #This function is repsonsible for loading in the selected file
      infile <- input$SNPdata
      if (is.null(infile)) {
        # User has not uploaded a file yet
        return(NULL)
      }
     # snps <- read.PLINK(input$SNPdata$datapath)
     snps <- read.PLINK("data/ORYX_500K_ld.raw")
      
      snps <- fread(input$SNPdata$datapath)
      snps <- fread("data/pine.txt")
      
      genin <- reactive({
        
        locus <- snps[, -c(1, 2, 3, 4, 17:3086)]    
        colnames(locus) <- gsub("\\.", "_", colnames(locus)) # locus names can't have "."
        ind <- as.character(snps$tree_id) # labels of the individuals
        population <- as.character(snps$state) # labels of the populations
        snps1 <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep = "")
        snps1
        
      })
      # summary table - this should be different
      # output$sum_stats <- renderTable({
      #   summary <- data.frame(c("Number of loci", "Number of individuals"),
      #                         c(snps@n.loc, length(snps@ind.names)))
      #   })
      
      dat <- list(toto=c(1,1,0,0), titi=c(NA,1,1,0), tata=c(NA,0,3, NA))
      x <- new("genlight", dat)
      x.mat <- as.matrix(x) # x is a genlight object
      x.mat[x.mat == 0] <- "1/1" # homozygote reference
      x.mat[x.mat == 1] <- "1/2" # heterozygote
      x.mat[x.mat == 2] <- "2/2" # homozygote alternate
      x.gid <- df2genind(x.mat, sep = "/", ploidy = 2)
      
      library(hierfstat)
      
      # calculate fst across all loci
      
      matFst <- pairwise.fst(x.gid)
      
      
      # df to genind
      df <- data.frame(locusA=c("11","11","12","32"),
                       locusB=c(NA,"34","55","15"),locusC=c("22","22","21","22"))
      row.names(df) <- .genlab("genotype",4)
      df
      
      obj <- df2genind(df, ploidy=2, ncode=1)
      obj
      tab(obj)
      
      
      df <- reactive({
        
        # Get maf
        
        x <- as.matrix(snps)
        SNP_ID <- colnames(x)
        x <- as.matrix(x)
        x <- t(x)
        
        ## calc_n
        n0 <- apply(x==0,1,sum,na.rm=T)
        n1 <- apply(x==1,1,sum,na.rm=T)
        n2 <- apply(x==2,1,sum,na.rm=T)
        
        n <- n0 + n1 + n2
        
        ## calculate allele frequencies
        p <- ((2*n0)+n1)/(2*n)
        q <- 1 - p
        MAF <- pmin(p, q)
        
        MAF <- data.frame(SNP_ID) %>%
          mutate(MAF = MAF)
        
        # Get geno rate
        df <- as.data.frame(x)
        df <- df[c(1:1000),]
        
        #test
        obj <- df2genind(t(df), ploidy=2, ncode=1)
        obj
        pop <- c("A", "A", "A", "A", "A",
                 "B", "B", "B", "B", "B",
                 "C", "C", "C", "C", "C",
                 "D", "D", "D", "D", "D")
        
        pop(obj) <- pop
        
        obj@pop
        matFst <- pairwise.fst(obj)
        
        
        
        df[df=="2"]<-NA
        nind <- ncol(df)
        ind_per_snp <- data.frame(n = ncol(x) - rowSums(is.na(df)))
        df$ind_per_snp <- ind_per_snp$n
        df$SNP_ID <- rownames(df)
        
        df <- df %>%
          mutate(geno_rate = (ind_per_snp/ncol(.))*100) %>%
          left_join(MAF)
        
        #x[c(1:10),(1:21)]
        
        df <- filter(df, geno_rate >= input$percent_geno & MAF > input$maf_thresh) # and maf
      #  df <- filter(df, geno_rate >= 90 & MAF > 0) # and maf
        
        
      })
      
      output$mytable = DT::renderDataTable({
        df()
      })
      
      output$sum_stats <- renderTable({
        summary <- data.frame(c("Number of individuals", "Number of loci"),
                              c((ncol(df())-4), nrow(df())))
        })
      
      output$sum_stats2 <- renderTable({
        summary <- data.frame(c("Number of individuals", "Number of loci"),
                              c((ncol(df())-4), nrow(df())))
      })
      
      output$geno <- renderPlot({
        ggplot(df(), aes(geno_rate)) +
          geom_histogram() +
          xlim(c(0,100))
       #   xlim(c(0,nind+1))
        
      })
      
      
      output$maf <- renderPlot({
        ggplot(df(), aes(MAF)) +
          geom_histogram() +
          xlim(c(0,0.5))
      })
      
      
    })
    

    
  }
  
)
      