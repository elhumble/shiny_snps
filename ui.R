library(shiny)
shinyUI(pageWithSidebar(

    titlePanel(
    title = div(img(src="shiny_snps_logo.jpg", height = 150, width = 400), align = "center"), windowTitle = "ShinySNPs"
    #title = h1("wildlife_detectR", align = "center"))
    #titlePanel((title=div("wildlife detectR",(img(src="heading.jpg", height = 150, width = 150)), align = "center"))
    ),

    # 
    # headerPanel(
    #   HTML('<img src="shiny_snps_logo.jpg" height="50%" width="50%" />'),
    #   windowTitle="ShinySNPs"
    # ),
  
    sidebarPanel(
    
    ## conditionalPanel() functions for selected tab
    conditionalPanel(condition="input.tabselected==1",h3("Upload data"),
                     fileInput("SNPdata", "SNP dataset (PLINK raw format)"),
                     sliderInput("percent_geno", "Select genotyping rate", min = 0, max = 100, value = 0, step = 1),
                     sliderInput("maf_thresh", "Minor allele frequency threshold", min = 0, max = 0.5, value = 0, step = 0.01)
                     ),
    
    conditionalPanel(condition="input.tabselected==2", h3("Hardy-Weinberg equilibrium"),
                     h4("Chi-Square values (chi^2), number of degrees of freedom (df) and the associated p-values (Pr(chi^2)) for each SNP")
                     ),
    
    conditionalPanel(condition="input.tabselected==3", h3("Population structure"),
                     h4("Pairwise FST values calculated for each population comparison"),
                     ("How do these values change as you filter SNPs for MAF and genotyping rate?")
                     ),
    conditionalPanel(condition="input.tabselected==4", h3("Inbreeding"),
                     h4("Distribution of multi-locus heterozygosity across individuals"),
                     ("How does this distribution vary as you filter the dataset for MAF and genotyping rate?")
    ),
    
    conditionalPanel(condition="input.tabselected==5", h3("Data"),
                     h4("List of SNPs and their minor allele frequencies and individual genotyping rates")
    )
      ),
  
  
  mainPanel(
    # recommend review the syntax for tabsetPanel() & tabPanel() for better understanding
    # id argument is important in the tabsetPanel()
    # value argument is important in the tabPanle()
    tabsetPanel(
      tabPanel("Summary", value=1,
              # helpText("Upload a dataset of SNPs in PLINK raw format"),
              # hr(),
               tableOutput("sum_stats"),
               #plotOutput("maf", width = "40%"),
               #plotOutput("geno", width = "40%"),
               splitLayout(
                 cellWidths = "50%",
                 cellArgs = list(style = "padding: 6px"),
                 plotOutput("maf"),
                 plotOutput("geno")
               )),

      tabPanel("Population structure", value=3,
               DT::dataTableOutput("fst"),
               renderPlot("pca")),
      tabPanel("Inbreeding", value =4,
               plotOutput("MLH", width = "50%")),
      #tabPanel("Relatedness", value = 5),
      tabPanel("Hardy-Weinberg equilibrium", value=2,
               DT::dataTableOutput("hwe")
      ),
      tabPanel("Data", value = 5,
               DT::dataTableOutput("mytable")),
      id = "tabselected"
    )
  )
))