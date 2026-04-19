library(terra)
library(stringr)

merge_uv<-function(u_dir,v_dir,outdir='merged_current'){
  u_files<-list.files(u_dir,pattern='re_current_u_\\d+_\\d+\\.tif$',full.names = TRUE)
  for(u_path in u_files){
    fname_u<-basename(u_path)
    depth_str<-str_extract(fname_u,"\\d+_\\d+")
    v_path<-file.path(v_dir,paste0("re_current_v_", depth_str, ".tif"))
    r_u<-rast(u_path)
    r_v<-rast(v_path)
    u_df<-as.data.frame(r_u,xy=TRUE)
    v_df<-as.data.frame(r_v,xy=TRUE)
    merged<-merge(u_df,v_df,by=c('x','y'),suffixes=c('_u','_v'))
    names(merged)<-c('lon','lat','u','v')
    merged$speed<-sqrt(merged$u^2+merged$v^2)
    out_path <- file.path(outdir, paste0("current_", depth_str, ".csv"))
    write.csv(merged,out_path)
  }
}

merge_uv('monthly_mean_by_depth','monthly_mean_by_depth',outdir="merged_current")


