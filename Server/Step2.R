shinyjs::hide('step2-dl1')
cor.df<-reactive({
  stat.df<-data.table(NE=SCLC.NE(),user.dat())[,lapply(.SD,function(x) unlist(cor.test(NE,x)[c('estimate','p.value')])),.SDcols=-1]
  stat.df<-data.table(colnames(user.dat()),t(stat.df));setnames(stat.df,c('gene','r','pv'));stat.df[,padj:=p.adjust(pv,'BH')]
  stat.df[,pv:=signif(pv,2)][,padj:=signif(padj,2)][,r:=round(r,2)];stat.df[order(r)]
})
low.cor<-reactive({min(cor.df()$r,na.rm = T) > -.5 | max(cor.df()$r,na.rm = T)<.5})
NE.sig<-reactive({
  data.frame(Symbol=c(cor.df()[order(r,decreasing = T)]$gene[1:25],cor.df()[order(r,decreasing = F)]$gene[25:1]),
             NE.weight=rep(c(1,0),each=25),nonNE.weight=rep(c(0,1),each=25))
})
observeEvent(input$`step2-btn1`,{
  output$`step3-tbl1`<-NULL;output$`step3-plot1`<-NULL
  show('step2-loader')
  toggle(id = 'step2-loader', condition = TRUE)
  Sys.sleep(1)
  output$`step2-tbl1`<-DT::renderDataTable(DT::datatable(data = isolate(cor.df()),options = list(lengthMenu = c(10,15))))
  if(low.cor()){
    output$`step2-hint2`<-renderText('Low correlation observed between NE score and gene expression data. NE signature may not be reliable!')
  }
  output$`step2-dl1` <- downloadHandler(
    filename=function(){'study-specific_NE_signature.csv'},content = function(file){write.csv(isolate(NE.sig()),file,row.names = F)}
  )
  shinyjs::show('step2-dl1')
  output$box3<-renderUI({
    box(width = 4,collapsed = T,
        title = "Calculate NE scores with study-specific NE signature", status = ifelse(small.NE.diff()|low.cor(),'warning',"primary"), solidHeader = TRUE,collapsible = F,
        actionButton('step3-btn1',label = "Calculate"),hr(),
        DT::dataTableOutput("step3-tbl1"),
        useShinyjs(),
        downloadButton('step3-dl1',label = "Download study-specific NE scores"),
        hr(),
        radioButtons('step3-rb1',label = "Compare SCLC and study-specific NE signatures",choices = c('scatter plot','heatmap'),inline = T),
        actionButton('step3-btn2',label = "Submit"),
        plotOutput("step3-plot1")
    )
  })
  delay(10,{shinyjs::hide('step3-dl1');shinyjs::hide('step3-rb1');shinyjs::hide('step3-btn2')})
})


