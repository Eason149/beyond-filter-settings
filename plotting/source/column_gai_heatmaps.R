self<-normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]))
ss<-readLines(file.path(dirname(self),'draw_figures.R'),warn=FALSE,encoding='UTF-8')
eval(parse(text=ss[1:(match('figures<-list()',ss)-1)]))
ps<-panels[['8']]
for(i in 1:4){
 ps[[i]]<-tag(ps[[i]],i)+labs(title=gsub('-','\u2013',windows[(i-1)%%2+1]))+
 theme(plot.title=element_text(size=7.2),plot.margin=margin(2,2,2,2))
 if(i>=3) ps[[i]]<-ps[[i]]+guides(fill=guide_colorbar(direction='horizontal',barwidth=unit(22,'mm'),barheight=unit(2,'mm'),display='rectangles'))
}
p<-wrap_plots(ps,ncol=2,heights=c(1,1.45),guides='collect') & theme(legend.position='bottom')
dest<-file.path(root,'review_layouts','IPM_column_F4_complete')
require_patchwork_panel_alignment(p,manifest_path=paste0(dest,'.layout.json'),report_path=paste0(dest,'.alignment.json'),width_in=80/25.4,height_in=135/25.4,audit_script=file.path(skill,'audit_panel_alignment.py'),python=python,tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=TRUE)
grDevices::cairo_pdf(paste0(dest,'.pdf'),width=80/25.4,height=135/25.4,family='Arial');print(p);dev.off()
pdftools::pdf_convert(paste0(dest,'.pdf'),dpi=155,filenames=paste0(dest,'.png'))
