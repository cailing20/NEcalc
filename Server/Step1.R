shinyjs::hide('step1-btn1')
shinyjs::hide('step1-dl2')
output$`step1-dl1` <- downloadHandler(
  filename=function(){'example_input.csv'},content = function(file){write.csv(example.dat,file,row.names = T)}
)
user.dat<-reactive({
  req(input$expr)
  tmp<-fread(input$expr$datapath);rnames<-tmp[[1]];tmp<-as.matrix(tmp[,-1,with=F]);rownames(tmp)<-rnames
  tmp
})
output$box2<-renderUI({
  box(width = 4,collapsed = T,
      title = "Generate a study-specific NE signature", status = "warning", solidHeader = TRUE,collapsible = TRUE,
      h6("This step requires the SCLC NE scores from step 1 as input.")
  )
})
output$box3<-renderUI({
  box(width = 4,collapsed = T,
      title = "Calculate NE scores with study-specific NE signature", status = "warning", solidHeader = TRUE,collapsible = TRUE,
      h6("This step requires the study-specific NE signature from step 2 as input.")
  )
})
observeEvent(input$expr,{
  shinyjs::hide('step1-dl2')
  shinyjs::hide('step1-btn1')
  output$`step1-plot`<-NULL
  output$`step1-tbl`<-NULL
  output$box2<-renderUI({
    box(width = 4,collapsed = T,
        title = "Generate a study-specific NE signature", status = "warning", solidHeader = TRUE,collapsible = TRUE,
        h6("This step requires the SCLC NE scores from step 1 as input.")
    )
  })
  output$box3<-renderUI({
    box(width = 4,collapsed = T,
        title = "Calculate NE scores with study-specific NE signature", status = "warning", solidHeader = TRUE,collapsible = TRUE,
        h6("This step requires the study-specific NE signature from step 2 as input.")
    )
  })
  human.gs<-reactive({length(intersect(colnames(user.dat()),sig$Symbol))})
  mouse.gs<-reactive({length(intersect(colnames(user.dat()),sig.mouse$Symbol))})
  row.gs<-reactive({length(intersect(rownames(user.dat()),c(sig.mouse$Symbol,sig$Symbol)))})
  if(human.gs()>0){
    if(human.gs()<25){
      output$`step1-hint`<-renderText(paste(human.gs(),'human NE signature genes detected in input data. More genes are needed.'))
    }else{
      output$`step1-hint`<-renderText(paste(human.gs(),'human NE signature genes detected in input data.'))
      shinyjs::show('step1-btn1')
    }
  }else if(mouse.gs()>25){
    if(mouse.gs()<25){
      output$`step1-hint`<-renderText(paste(human.gs(),'mouse NE signature genes detected in input data. More genes are needed.'))
    }else{
      output$`step1-hint`<-renderText(paste(human.gs(),'mouse NE signature genes detected in input data.'))
      shinyjs::show('step1-btn1')
    }
  }else if(row.gs()>0){
    if(row.gs()<25){
      output$`step1-hint`<-renderText(paste(row.gs(),'NE signature genes detected in input data. More genes are needed.'))
    }else{
      output$`step1-hint`<-renderText(paste(row.gs(),'NE signature genes detected in input data. Please transpose input matrix and try upload again'))
    }
  }else{
    output$`step1-hint`<-renderText('no NE signature genes detected in input data. Please reference our example data and prepare your data in the required format.')
  }
})

SCLC.NE<-reactive({
  comp.score(dat=t(user.dat()),sig = switch(as.numeric(identical(colnames(user.dat())[1:4],toupper(colnames(user.dat())[1:4])))+1,sig.mouse,sig),logged = max(user.dat()[1,],na.rm = T)<100)
})
NE.df<-reactive({data.table(sample=rownames(user.dat()),`NE score`=round(SCLC.NE(),2))})

small.NE.diff<-reactive({diff(range(SCLC.NE(),na.rm = T))<.65})
observeEvent(input$`step1-btn1`,{
  output$`step1-tbl`<-DT::renderDataTable(DT::datatable(data = isolate(NE.df()),options = list(lengthMenu = c(5,10,15))))
  gg.NE<-reactive({
    ggplot(data.table(study='data',NE.df()),aes(x=`NE score`,y=study,fill=..x..))+geom_density_ridges_gradient(panel_scaling = F)+
      scale_fill_gradientn("NE score",colours = c(low.NE.col,'white',high.NE.col),breaks = c(-1,0,1),limits=c(-1,1))+
      theme_bw()+theme(panel.grid = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.key.size = unit(0.5, 'cm'))+
      xlab('NE score')+ylab('density')+ggtitle('NE score distribution (by SCLC NE signature)')+xlim(-1,1)
  })
  output$`step1-plot`<-renderPlot({isolate(gg.NE())},height = 100)
  shinyjs::show('step1-dl2')
  output$box2<-renderUI({
    box(width = 4,collapsed = T,
        title = "Generate a study-specific NE signature", status = ifelse(small.NE.diff(),'warning',"primary"), solidHeader = TRUE,collapsible = F,
        h6("Compute correlation between NE score and input expression data"),
        useShinyjs(),
        actionButton('step2-btn1',label = "Run"),
        
        hidden(div(id = 'step2-loader', hr(),withSpinner(DT::dataTableOutput("step2-tbl1")),hr())),
        textOutput('step2-hint1'),
        textOutput('step2-hint2'),
        downloadButton('step2-dl1',label = "Download study-specific NE signature"),
    )
  })
  delay(10,shinyjs::hide('step2-dl1'))
  if(small.NE.diff()) output$`step2-hint1`<-renderText('Low NE heterogeneity within sample set. Results may not be reliable!')
})
output$`step1-dl2` <- downloadHandler(
  filename=function(){'SCLC_signature_based_NE_scores.csv'},content = function(file){write.csv(isolate(NE.df()),file,row.names = F)}
)