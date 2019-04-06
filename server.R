library(shiny)
library(ape)
library(data.table)
library(adegenet)
library(ggplot2)
library(dplyr)
library(inbreedR)
library(hierfstat)
library(pegas)

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
     input_data <- read.PLINK("data/ORYX_500K_ld.raw")
      
     # snps <- fread(input$SNPdata$datapath)
     # snps <- fread("data/pine.txt")
      
      # dataframe to genind
      
      
      # genind <- reactive({
      #   
      #   locus <- snps[, -c(1, 2, 3, 4, 17:3086)]    
      #   colnames(locus) <- gsub("\\.", "_", colnames(locus)) # locus names can't have "."
      #   ind <- as.character(snps$tree_id) # labels of the individuals
      #   population <- as.character(snps$state) # labels of the populations
      #   snps1 <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep = "")
      #   snps1
      #   
      # })
      # summary table - this should be different
      # output$sum_stats <- renderTable({
      #   summary <- data.frame(c("Number of loci", "Number of individuals"),
      #                         c(snps@n.loc, length(snps@ind.names)))
      #   })

     
      
      # calculate fst across all loci
      
      # matFst <- pairwise.fst(x.gid)
      
      
      # df to genind
      # df <- data.frame(locusA=c("11","11","12","32"),
      #                  locusB=c(NA,"34","55","15"),locusC=c("22","22","21","22"))
      # row.names(df) <- .genlab("genotype",4)
      # df
      # 
      # obj <- df2genind(df, ploidy=2, ncode=1)
      # obj
      # tab(obj)
      
     
     x <- as.matrix(input_data)
     df <- as.data.frame(t(x))
     
     df <- df[c(1:100),]
     SNP_ID <- rownames(df)
     #df[df=="2"]<-NA
     
     maf <- df %>% mutate(n0 = rowSums(. == 0, na.rm = T),
                          n1 = rowSums(. == 1, na.rm = T), # calc_n
                          n2 = rowSums(. == 2, na.rm = T),
                          n = n0 + n1 + n2,
                          p =  ((2*n0)+n1)/(2*n), # calc allele f
                          q = (1 - p),
                          MAF = pmin(p, q)) %>%
       mutate(SNP_ID = SNP_ID) %>%
       select(SNP_ID, MAF)
     
     snp_geno <- df %>%
       mutate(ind_per_snp = ncol(.) - rowSums(is.na(.)),
              geno_rate = (ind_per_snp / ncol(.))*100) %>%
       mutate(SNP_ID = SNP_ID)
     
     maf_geno <- left_join(maf, snp_geno, by = "SNP_ID")
     
     
      df <- reactive({
        
        df <- maf_geno %>%
          select(SNP_ID, MAF, geno_rate) %>%
         # filter(geno_rate > 0 & MAF > 0.2)
           filter(geno_rate >= input$percent_geno & 
                   MAF > input$maf_thresh)

      })
      
      genind <- reactive({ # needs to be the same as df
        
        SNP_ID <- df()$SNP_ID
        
        genind <- df() %>%
          select(-c(MAF, geno_rate)) %>%
          left_join(maf_geno, by = "SNP_ID") %>%
          select(-c(SNP_ID, MAF, SNP_ID, ind_per_snp, geno_rate)) %>%
          t(.)
        
        genind[genind == 0] <- "1/1" # homozygote reference
        genind[genind == 1] <- "1/2" # heterozygote
        genind[genind == 2] <- "2/2" # homozygote alternate
        
        colnames(genind) <- SNP_ID
        inds <- rownames(genind)
        genind <- df2genind(genind, ploidy=2, sep = "/", NA.char = "NA")
        indNames(genind)
        genind <- genind

      })

      
      
      #~~ Summary table
      output$sum_stats <- renderTable({
        
        summary <- data.frame(c("Number of individuals", "Number of loci",
                                "Missing data", "He", "Ho"),
                              c(ncol(maf_geno)-4, nrow(df()), 
                                summary(genind())$NA.perc,
                                mean(summary(genind())$Hexp),
                                mean(summary(genind())$Hobs)))

        # summary <- data.frame(c("Number of individuals", "Number of loci",
        #                         "Missing data", "He", "Ho"),
        #                       c(ncol(maf_geno)-4, nrow(df), ade_sum$NA.perc,
        #                       mean(ade_sum$Hexp), mean(ade_sum$Hobs)))
        # colnames(summary) <- c("Summary statistics","")
        
      })


      #hwe <- hw.test(genind(), B=0)
      # test <- as.matrix(genind)
      # test[c(1:10),c(1:10)]
      
      output$hwe <- DT::renderDataTable({
        hwe <- hw.test(genind(), B=0)

      })
      
      
      #
      # 
      # pop <- c("A", "A", "A", "A", "A",
      #          "B", "B", "B", "B", "B",
      #          "C", "C", "C", "C", "C",
      #          "D", "D", "D", "D", "D")
      # 
      # pop(genind) <- pop
      # 
      # genind@pop
      # matFst <- pairwise.fst(genind)
      # 
      # 
      output$mytable = DT::renderDataTable({
        df()
      })
      

      
      output$sum_stats2 <- renderTable({
        summary <- data.frame(c("Number of individuals", "Number of loci"),
                              c((ncol(df())-4), nrow(df())))
      })
      
      output$geno <- renderPlot({
        ggplot(df(), aes(geno_rate)) +
          geom_histogram() +
          xlim(c(0,110))
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
      