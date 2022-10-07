library(shiny)
library(shinyjs)
library(shinydashboard)
library(ggplot2)
library(ggridges)
library(DT)
library(data.table)
library(shinycssloaders)
library(ComplexHeatmap)
library(RColorBrewer)
high.NE.col<-"#5e3c99";low.NE.col<-"#e66101"
example.dat<-readRDS('data/example_expr.rds')
sig<-readRDS('data/SCLC_sig.rds')
sig.mouse<-readRDS('data/SCLC_mmu_sig.rds')
comp.score<-function(dat, sig, logged) {
  ind <- match(sig$Symbol, rownames(dat))
  ind <- ind[!is.na(ind)]
  sig <- sig[sig$Symbol %in% rownames(dat),]
  if (!logged) {
    dat[dat < 1] <- 1
    dat <- log2(dat)
  }
  score1 <- apply(dat[ind, ], 2, function(d) tryCatch({cor.test(sig[, 2], d)$estimate},error=function(e) NA))
  score2 <- apply(dat[ind, ], 2, function(d) tryCatch({cor.test(sig[, 3], d)$estimate},error=function(e) NA))
  score <- (score1 - score2)/2
  score
}
ui <- dashboardPage(
  dashboardHeader(title = "NE score calculator"),
  dashboardSidebar(collapsed = T),
  dashboardBody(
    fluidRow(
      box(width = 4,
          title = "Calculate NE scores with SCLC NE signature", status = "primary", solidHeader = TRUE,collapsible = TRUE,
          downloadButton('step1-dl1', label = "Download example input"),
          h6('Prepare transcriptomic data with samples in rows and gene name in columns. You may also download the example data (neuroblastoma cell line microarray data) to test use our tools.'),
          fileInput("expr", "Upload expression data",multiple = TRUE,
                    accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
          textOutput(outputId = 'step1-hint'),
          useShinyjs(),
          actionButton('step1-btn1',label = "Calculate"),
          plotOutput(outputId = "step1-plot", height = 100),
          DT::dataTableOutput(outputId = "step1-tbl"),
          downloadButton('step1-dl2', label = "Download NE scores")
      ),
      uiOutput("box2"),
      uiOutput("box3")
    )
  )
)

server <- function(input, output) {
  options(shiny.maxRequestSize=300*1024^2) 
  source(file.path("Server", "Step1.R"),  local = TRUE)$value
  source(file.path("Server", "Step2.R"),  local = TRUE)$value
  source(file.path("Server", "Step3.R"),  local = TRUE)$value
}

shinyApp(ui, server)
