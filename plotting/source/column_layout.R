args0<-commandArgs(FALSE)
self<-normalizePath(sub('^--file=','',grep('^--file=',args0,value=TRUE)[1]))
ss<-readLines(file.path(dirname(self),'draw_figures.R'),warn=FALSE,encoding='UTF-8')
eval(parse(text=ss[1:(match('figures<-list()',ss)-1)]))
out<-file.path(root,'review_layouts')
for (i in 4:6){
 spec<-list(c('8',1,2),c('9',3,4),c('10',2,4))[[i-3]]
 ps<-panels[[spec[1]]][as.integer(spec[-1])]
 p<-wrap_plots(lapply(1:2,function(j)tag(ps[[j]],j)),ncol=1)
 dest<-file.path(out,paste0('IPM_column_F',i))
 require_patchwork_panel_alignment(p,manifest_path=paste0(dest,'.layout.json'),report_path=paste0(dest,'.alignment.json'),width_in=78/25.4,height_in=106/25.4,audit_script=file.path(skill,'audit_panel_alignment.py'),python=python,tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=TRUE)
 grDevices::cairo_pdf(paste0(dest,'.pdf'),width=78/25.4,height=106/25.4,family='Arial');print(p);dev.off()
 pdftools::pdf_convert(paste0(dest,'.pdf'),dpi=130,filenames=paste0(dest,'.preview.png'))
}
