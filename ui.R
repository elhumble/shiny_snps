library(shiny)
shinyUI(pageWithSidebar(
  titlePanel(
    (title = div(img(src="heading_temp.jpg", height = 150, width = 400), align = "center"))
    #title = h1("wildlife_detectR", align = "center"))
    #titlePanel((title=div("wildlife detectR",(img(src="heading.jpg", height = 150, width = 150)), align = "center"))
    ),
  
    
    sidebarPanel(
    
    ## conditionalPanel() functions for selected tab
    conditionalPanel(condition="input.tabselected==1",h4("Upload data"),
                     fileInput("SNPdata", "SNP dataset"),
                     sliderInput("percent_geno", "Select genotyping rate", min = 0, max = 100, value = 0, step = 1),
                     sliderInput("maf_thresh", "Minor allele frequency threshold", min = 0, max = 1, value = 0, step = 0.01)    ),
    
    conditionalPanel(condition="input.tabselected==2", h3("Hardy-Weinberg equilibrium"),
                     h5("gkdgnkfdng fkdgnkd fndkslnf  fsdjklfj fl jklfjdskl fsfkl")
                     ),
    
    conditionalPanel(condition="input.tabselected==3", h3("Filtering"),
                     radioButtons("method2","Choose a method", choices=c("random match probability" = 1, "likelihood ratio approach" = 2)),
                     textInput("IndividualizationTheta", "Theta", ""),
                     fileInput("UnknownProfile", "Unknown profile file"),
                     fileInput("ReferenceProfiles", "Reference database file"),
                     selectInput("ProfileInput", "File type", choices=c("Genetix", "Genalex", "Genepop", "DArT"), selected = "Genepop")                     ),
    
    
    conditionalPanel(condition="input.tabselected==6", h3("Data")
    )
      ),
  
  
  mainPanel(
    # recommend review the syntax for tabsetPanel() & tabPanel() for better understanding
    # id argument is important in the tabsetPanel()
    # value argument is important in the tabPanle()
    tabsetPanel(
      tabPanel("About", value=1,
               tableOutput("sum_stats"),
               #plotOutput("maf", width = "40%"),
               #plotOutput("geno", width = "40%"),
               splitLayout(
                 cellWidths = "50%",
                 cellArgs = list(style = "padding: 6px"),
                 plotOutput("maf"),
                 plotOutput("geno")
               )),
      tabPanel("Hardy Weinberg Equilibrium", value=2,
               DT::dataTableOutput("hwe")
      ),
      tabPanel("Population structure", value=3),
      tabPanel("Inbreeding", value =4),
      tabPanel("Relatedness", value = 5),
      tabPanel("Data", value = 6,
               DT::dataTableOutput("mytable")),
      id = "tabselected"
    )
  )
))