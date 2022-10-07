new.SCLC.NE<-reactive({
  comp.score(dat=t(user.dat()),sig = NE.sig(),logged = max(user.dat()[1,],na.rm = T)<100)
})
new.NE.df<-reactive({data.table(sample=rownames(user.dat()),`NE score`=round(new.SCLC.NE(),2))})
observeEvent(input$`step3-btn1`,{
  output$`step3-tbl1`<-DT::renderDataTable(DT::datatable(data = isolate(new.NE.df()),options = list(lengthMenu = c(5,10,15))))
  shinyjs::show('step3-dl1')
  shinyjs::show('step3-rb1')
  shinyjs::show('step3-btn2')
})
m.NE<-reactive({m.df<-merge(NE.df(),new.NE.df(),by='sample');colnames(m.df)[2:3]<-c('SCLC signature','study-specific signature');m.df})
hm<-reactive({
  dat1<-user.dat()[,intersect(switch(as.numeric(identical(colnames(user.dat())[1:4],toupper(colnames(user.dat())[1:4])))+1,sig.mouse,sig)$Symbol,colnames(user.dat()))]
  dat1[]<-apply(dat1,2,scale)
  r1<-as.numeric(data.table(NE=SCLC.NE(),dat1)[,lapply(.SD,function(x) cor.test(NE,x)$estimate),.SDcols=-1])
  dat2<-user.dat()[,match(NE.sig()$Symbol,colnames(user.dat()))]
  dat2[]<-apply(dat2,2,scale)
  r2<-as.numeric(data.table(NE=new.SCLC.NE(),dat2)[,lapply(.SD,function(x) cor.test(NE,x)$estimate),.SDcols=-1])
  draw(Heatmap(t(dat1[order(new.SCLC.NE()),]),left_annotation = HeatmapAnnotation(which = 'row',r=r1,col = list(r=circlize::colorRamp2(c(-1,-0.5,0,0.5,1),rev(brewer.pal(5,'RdBu'))))),
               top_annotation = HeatmapAnnotation(which = 'column',`NE score`=SCLC.NE()[order(new.SCLC.NE())],col = list(`NE score`=circlize::colorRamp2(c(-1,-0.5,0,0.5,1),brewer.pal(5,'PuOr')))),
               cluster_rows = F,cluster_columns = F,show_row_names = F,name = 'z-transformed expr',row_title = 'SCLC NE signature')%v%
         Heatmap(t(dat2[order(new.SCLC.NE()),]),left_annotation = HeatmapAnnotation(which = 'row',r=r2,col = list(r=circlize::colorRamp2(c(-1,-0.5,0,0.5,1),rev(brewer.pal(5,'RdBu'))))),
                 top_annotation = HeatmapAnnotation(which = 'column',`NE score`=sort(new.SCLC.NE()),col = list(`NE score`=circlize::colorRamp2(c(-1,-0.5,0,0.5,1),brewer.pal(5,'PuOr')))),
                 show_column_names = length(new.SCLC.NE())<100,column_names_gp = gpar(fontsize=max(round((132-length(new.SCLC.NE()))/9),4)),
                 cluster_rows = F,cluster_columns = F,show_row_names = F,name = 'z-transformed expr',row_title ='study-specific NE signature'),merge_legend=T,column_title='expression of NE signature genes')
})


observeEvent(input$`step3-btn2`,{
  if(input$`step3-rb1`=='scatter plot'){
    output$`step3-plot1`<-renderPlot(ggplot(isolate(m.NE()),aes(x=`SCLC signature`,y=`study-specific signature`))+geom_point(alpha=.3)+ggtitle('NE score')+theme_bw()+xlim(-1,1)+ylim(-1,1))
  }else{
    output$`step3-plot1`<-renderPlot(isolate(hm()))
  }
})
output$`step3-dl1` <- downloadHandler(
  filename=function(){'study-specific_signature_based_NE_scores.csv'},content = function(file){write.csv(isolate(new.NE.df()),file,row.names = F)}
)