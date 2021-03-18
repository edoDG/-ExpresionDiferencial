#------------------------------------------------------------------------------------
# 17 / 03 / 21                                        Expresion diferencial
#------------------------------------------------------------------------------------
#Edoardo Cruz

library("sleuth") #La libreria para realizar la actividad

#Una funcion con la cual podemos acceder a los daros necesarios para la actividad (de biomanRt)
tx2gene <- function(){
  
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") 
  #Aqu? lo que usamos es hsapiens_gene_ensembl" porque es donde vamos a usar los datos.
  
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = "ensembl_transcript_id",
                       ens_gene = "ensembl_gene_id", ext_gene = "external_gene_name")
  return(t2g)
}

t2g <- tx2gene() #Conectarse a la base de datos de ensambled
t2g


base_dir <-"~/Bioinfo/Archivo/" #Este es el directorio donde se encuentran tus archivos
#o sea los samples, esto se cambia 

#Selecci?n de los samples
#De que quiero los samples 1,2,3,10,11 y 12, entonces esos son los que ponemos en el concatenar
samples <- paste0("sample", c("1","2","3",
                                   "10","11","12")) 
samples #aqu? vemos que samples seleccionamos 

kal_dirs <- sapply(samples, function(id) file.path(base_dir, id)) 
#Entonces aqui lo que se hace es darle a cada ID la direcccion de donde esta pues el sample 

#Un data drame donde se colocan las direcciones, pues las samples (que ya fueron definidas) 
#y el nombre que se les van a colocar a las muestras
#en este caso la sample 1, 2 y 3 seran el control, mientras la 10, 11 y 12 las experimentales. 
s2c <- data.frame(path=kal_dirs, sample=samples, muestras <- c("control","control","control", "exp",
                                                              "exp", "exp"), stringsAsFactors=FALSE)

#
so <- sleuth_prep(s2c, full_model=~muestras, target_mapping = t2g, extra_bootstrap_summary = TRUE) 
#se almacenan datos para poder trabajar con ellos
#s2c es el dataframe
#formula que explica el dise?o del modelo
#La dataframe con los ID de los genes.

so <- sleuth_fit(so) #ajustar el error 
so <- sleuth_wt(so, which_beta = "muestrasexp")               
sleuth_live(so) #visualizacion con shiny

#Seleccionamos donde se encuentra la tabla con los resultados 
setwd("~/Bioinfo/Archivo/")
resultados<-read.table("test_table.csv",sep=",",
                       header=TRUE) #se cargan a R

#Busqueda y selección de los resultados que son significativos 
significativos<-which(resultados$qval<0.1)
significativos<-resultados[significativos,]

#Resultados que son aquellos que incrementaron 
upregulated<-which(significativos$b>0)
upregulated<-significativos[upregulated,]

#aquellos que decrecieron 
downregulated<-which(significativos$b<0)
downregulated<-significativos[downregulated,]

#vamos a guardarlos en tablas 
write.table(upregulated,file="~/Bioinfo/Archivo/UpRConvsexp",sep="\t")
write.table(downregulated,file="~/Bioinfo/Archivo/DownRConvsexp.txt",sep="\t")

#ver cuales son los que se regularon positivamente
upregulated$ext_gene