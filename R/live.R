#'Volcano plot in live Shiny session
#'
#'\code{plotVolcano_live} makes a volcano plot from limma results.
#'
#' @param pdat Protein \code{proteusData} object.
#' @param res Result table from  \code{\link{limmaDE}}.
#' @param max_points Maximum number of points that can be selected.
#'
#' @return A \code{shiny} session hosted locally.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' prodat.med <- normalizeData(prodat)
#' res <- limmaDE(prodat.med)
#' plotVolcano_live(prodat.med, res)
#'
#'@export
plotVolcano_live <- function(pdat, res, max_points=10){
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  manage.pkg<-function(package_name){
    if(!package_name%in%installed.packages()){
      install.packages(package_name)
    }
    library(package_name,character.only = TRUE)
  }
  
  # Load required packages
  lapply(c("shiny","ggplot2","dplyr","DT","gplots"), function(x) manage.pkg(x))
  
  
  # Marek plotVolcano function to be enhanced with Shiny.
  # pdat<-protdat.annot
  DEdat<-res
  DEdat$"-log10(P.Value)" <- -log10(DEdat$P.Value)
  uniprot_ids <- sapply(strsplit(as.character(DEdat$protein),"|",fixed = TRUE),function(x) x[2])
  urls <- paste0("https://www.uniprot.org/uniprot/",uniprot_ids)
  DEdat$urls <- urls
  
  
  #######################################################################
  
  ui <- shinyUI(fluidPage(
    
    # Application title
    titlePanel("PlotVolcano live version"),
    
    fluidRow(
      column(5, plotOutput("plotVolcano", height = "700px", width = "100%", brush = "plot_brush",hover="plot_hover")),
      column(7,
             fluidRow(htmlOutput("proteinInfo")),
             fluidRow(
               column(4,
                      radioButtons("intensityScale","Intesity Scale:",choices = c("Linear scale" = "","Log scale"="Log"),inline = TRUE)
               )
             ),
             fluidRow(
               column(4,
                      fluidRow(plotOutput("jitterPlot", height = "400px",width = "100%")),
                      fluidRow(htmlOutput("gap")),
                      fluidRow(tableOutput("significanceTable"))
               ),
               column(4,
                      fluidRow(tableOutput("replicateTable"))
               ),
               column(4,
                      fluidRow(plotOutput("heatMap",height = "700px" ,width = "100%"))
               )
             )
      ),
      fluidRow(
        column(6, htmlOutput("proteinTable"))
      )
    ),
    
    # Show main protein table
    fluidRow(
      column(width = 12,
             DT::dataTableOutput("allProteinTable"))
    )
  )
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output) {
    
    #function to fetch selected proteins from Volcano plot or table
    selectProtein <- function(data,max_hover=1){
      # print('selectProtein method')
      sel = -1
      tab_idx <- as.numeric(input$allProteinTable_rows_selected)
      # print(paste0('tab_idx= ',tab_idx))
      if(!is.null(input$plot_brush)){
        brushed <- na.omit(brushedPoints(DEdat, input$plot_brush))
        sel <-as.numeric(rownames(brushed))
      }else if(!is.null(input$plot_hover)){
        # print(input$plot_hover)
        near <- nearPoints(DEdat,input$plot_hover,threshold = 20,maxpoints = max_hover)
        sel <- as.numeric(rownames(near))
      }else if(length(tab_idx)>0){
        sel <- tab_idx
        # print('  tab')
      }
      # print('sel')
      # print(sel)
      return(sel)
    }
    
    #ProteinInfo
    output$proteinInfo <- renderUI({
      # print('proteinInfo method')
      sel <- selectProtein(pdat$tab)
      n <- length(sel)
      if (n == 1 && sel > 0){
        name <- paste0('<H3>', as.character(pdat$annotation$protein[sel]),'</H3>')
        descr <- as.character(pdat$annotation$name[sel])
        HTML(paste0(name,descr,'<hr/>'))
      }else if (n > 1 && n <= max_points){
        HTML(paste0('<H3>','selection of ', n, ' proteins', '</H3><hr/>'))                    
      }else if (n > max_points){
        HTML(paste0('<H3>','only ',max_points,' points can be selected', '</H3><hr/>'))
      }
    })
    
    output$gap <- renderUI({HTML('<br/>')})
    
    # replicateTable
    output$replicateTable <- renderTable({
      # print('replicateTable method')
      sel <- selectProtein(pdat$tab)
      if(length(sel) > 1 && sel > 0 && length(sel) <= max_points){
        data.frame(Sample=colnames(pdat$tab),Intensity=colMeans(pdat$tab[sel,],na.rm = TRUE))
      }
      if (length(sel) == 1 && sel > 0){
        data.frame(Sample=colnames(pdat$tab),Intensity=pdat$tab[sel,])
      }
      if (length(sel) > max_points && sel > 0){ return() }
    },digits = 0, width = "80px"
    )
    
    #significanceTable
    output$significanceTable <- renderTable({
      sel <- selectProtein(pdat$tab)
      if(length(sel) == 1 && sel > 0){
        data.frame(Contrast="1112-BMO", pvalue=as.numeric(DEdat$P.Value[sel]))
      }
    },digits = 4,width = "50px"
    )
    
    #heatMap
    output$heatMap <- renderPlot({
      # print('heatMap method')
      sel<-selectProtein(pdat$tab)
      if(length(sel) > 1 && sel > 0 && length(sel) <= max_points){
        d <- as.matrix(pdat$tab[sel,])
        mean <- rowMeans(d,na.rm = TRUE)
        mean[mean == 0] <- 1
        d <- d/mean
        d[is.nan(d)] <- NA
        row.names(d) <- pdat$annotation$protein[sel]
        # print(d)
        heatmap.2(d, na.rm=TRUE, dendrogram = "row",key=FALSE,keysize = 1,lhei = c(1,100),Colv = FALSE,srtRow = -35,cexRow = 1.0,na.color = "blue")
      }
    })
    
    #jitterPlot
    output$jitterPlot <- renderPlot({
      # print('jitterPlot method')
      sel <- selectProtein(pdat$tab)
      if(length(sel)>0 && sel > 0 && length(sel) <= max_points){
        dataIntensity<-pdat$tab[sel,]
        if (input$intensityScale == 'Log'){
          dataIntensity<-log10(dataIntensity)
        }
        if (length(sel) == 1){
          i <- dataIntensity
        }else{
          i <- colMeans(dataIntensity,na.rm = TRUE)
        }
        s <- sapply(as.data.frame(dataIntensity),function(x) sd(x,na.rm = TRUE)/sqrt(length(x)))
        n <- length(sel)
        p <- data.frame(
          intensity= i,
          lo=i-s,
          up=i+s,
          condition=factor(pdat$metadata$condition,levels=unique(pdat$condition)),
          replicates=factor(pdat$metadata$replicate)
        )
        #shape
        p$shape <- rep(21,length(p$intensity))
        p$shape[which(p$intensity==0)] <- 24
        pd <- position_dodge(width = 0.4)
        #colorblind friendly definition
        cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        ggplot(p, aes(x=condition, y=intensity, ymin=lo, ymax=up,colour=replicates,  shape= shape, fill=replicates)) +
          theme(text = element_text(size=20), legend.position = "bottom",legend.direction = "horizontal") +
          {if (input$intensityScale == '') ylim(0, NA)} +
          geom_point(position=pd, size=4) +
          {if(n > 1) geom_errorbar(position=pd, width = 0.1)} +
          scale_shape_identity() +  # necessary for shape mapping
          scale_fill_manual(values=cbPalette) +
          {if (input$intensityScale == 'Log') labs(x = 'Condition', y = 'Log Intensity') else labs(x = 'Condition', y = 'Intensity')}
      }
    })
    #Volcano plot
    output$plotVolcano <- renderPlot({
      # print('plotVolcano method')
      tab_idx <- as.numeric(input$allProteinTable_rows_selected)
      pVol <- plotVolcano(DEdat,binhex=FALSE)
      if (length(tab_idx) > 0){
        pVol <- pVol + geom_point(data=DEdat[tab_idx,],size=3,color='red')
      }
      pVol
    })
    
    #AllProteinTable
    output$allProteinTable <-DT::renderDataTable({
      # print('allProteinTable method')
      d <- data.frame(ProteinId=DEdat$protein,Protein_Name=pdat$annotation$`PROTEIN-NAMES`,mean_1112=formatC(DEdat$mean_1112),mean_BMO=formatC(DEdat$mean_BMO))
      datatable(
        d, class = 'cell-border strip hover'
      ) %>% formatStyle(0, cursor = 'pointer')
    })
    
    #Open browser for Uniprot
    observeEvent(input$allProteinTable_cell_clicked, {
      info = input$allProteinTable_cell_clicked
      # do nothing if not clicked yet, or the clicked cell is not in the 1st column
      if (is.null(info$value) || info$col != 1) return()
      browseURL(DEdat$urls[info$row])
    })
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
} 



#' Fold-change intensity diagram in live Shiny session
#'
#' \code{plotFID_live} makes a log10 fold change versus log10 sum intensity plot,
#' usually known as MA plot.
#'
#' @param pdat Protein \code{proteusData} object.
#' @param res Result table from  \code{\link{limmaDE}}.
#' @param max_points Maximum number of points that can be selected.
#'
#' @return A \code{shiny} session hosted locally.
#'
#' @examples
#' library(proteusUnlabelled)
#' data(proteusUnlabelled)
#' prodat.med <- normalizeData(prodat)
#' res <- limmaDE(prodat.med,sig.level = 0.05)
#' plotFID_live(prodat.med, res)
#'
#' @export
plotFID_live <- function(pdat, res, max_points=10){
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  
  manage.pkg <- function(package_name){
    if(!package_name%in%installed.packages()){
      install.packages(package_name)
    }
    library(package_name,character.only = TRUE)
  }
  
  # Load required packages
  lapply(c("shiny","ggplot2","dplyr","DT","gplots"), function(x) manage.pkg(x))
  
  # Marek plotFID function to be enhanced with Shiny.
  # pdat<-protdat.annot
  DEdat <- res
  DEdat$"-log10(P.Value)" <- -log10(DEdat$P.Value)
  uniprot_ids <- sapply(strsplit(as.character(DEdat$protein),"|",fixed = TRUE),function(x) x[2])
  urls <- paste0("https://www.uniprot.org/uniprot/",uniprot_ids)
  DEdat$urls <- urls
  
  #Generate the FID datasets for selection. Code from Marek plotFID function.
  condMeans <- function(cond) {
    m <- rowMeans(log10(pdat$tab)[,which(pdat$metadata$condition == cond), drop=FALSE], na.rm=TRUE)
    m[is.nan(m)] <- NA
    m
  }
  
  m1 <- condMeans(pdat$conditions[1])
  m2 <- condMeans(pdat$conditions[2])
  good <- !is.na(m1) & !is.na(m2)
  FDIdf <- data.frame(
    id = rownames(pdat$tab),
    x = (m1 + m2) / 2,
    y = m2 - m1,
    good = good
  )
  rownames(FDIdf) <- 1:nrow(FDIdf)
  
  mx <- 1.1 * max(abs(FDIdf$y), na.rm=TRUE)
  m <- rbind(m1[!good], m2[!good])
  FDIdf[!good, "x"] <- colSums(m, na.rm=TRUE)
  FDIdf[!good, "y"] <- ifelse(is.na(m[1,]), mx, -mx)
  
  #######################################################################
  
  ui <- shinyUI(fluidPage(
    
    # Application title
    titlePanel("PlotFID live version"),
    
    fluidRow(
      column(5, plotOutput("plotFID", height = "700px", width = "100%", brush = "plot_brush",hover="plot_hover")),
      column(7,
             fluidRow(htmlOutput("proteinInfo")),
             fluidRow(
               column(4,
                      radioButtons("intensityScale","Intesity Scale:",choices = c("Linear scale" = "","Log scale"="Log"),inline = TRUE)
               )
             ),
             fluidRow(
               column(4,
                      fluidRow(plotOutput("jitterPlot", height = "400px",width = "100%")),
                      fluidRow(htmlOutput("gap")),
                      fluidRow(tableOutput("significanceTable"))),
               column(4,
                      fluidRow(tableOutput("replicateTable"))),
               column(4,
                      fluidRow(plotOutput("heatMap",height = "700px", width = "100%")))
             )
      ),
      fluidRow(
        column(6, htmlOutput("proteinTable"))
      )
    ),
    
    # Show main protein table
    fluidRow(
      column(width = 12,
             DT::dataTableOutput("allProteinTable"))
    )
  )
  )
  
  # Define server logic required to draw a histogram
  server <- function(input, output) {
    
    #function to fetch selected proteins from Volcano plot or table
    selectProtein <- function(data,max_hover=1){
      # print('selectProtein method')
      sel = -1
      tab_idx <- as.numeric(input$allProteinTable_rows_selected)
      # print(paste0('tab_idx= ',tab_idx))
      if(!is.null(input$plot_brush)){
        brushed <- na.omit(brushedPoints(FDIdf, input$plot_brush))
        sel <- as.numeric(rownames(brushed))
      }else if(!is.null(input$plot_hover)){
        near <- nearPoints(FDIdf,input$plot_hover,threshold = 20,maxpoints = max_hover)
        sel <- as.numeric(rownames(near))
      }else if(length(tab_idx)>0){
        sel <- tab_idx
        # print('  tab')
      }
      # print('sel')
      # print(sel)
      return(sel)
    }
    
    #ProteinInfo
    output$proteinInfo <- renderUI({
      # print('proteinInfo method')
      sel <- selectProtein(pdat$tab)
      # print (paste0('protein length= ', length(pdat$annotation$protein)))
      n <- length(sel)
      if (n == 1 && sel > 0){
        name <- paste0('<H3>', as.character(pdat$annotation$protein[sel]),'</H3>')
        descr <- as.character(pdat$annotation$name[sel])
        HTML(paste0(name,descr,'<hr/>'))
      }else if (n > 1 && n <= max_points){
        HTML(paste0('<H3>','selection of ', n, ' proteins', '</H3><hr/>'))                    
      }else if (n > max_points){
        HTML(paste0('<H3>','only ',max_points,' points can be selected', '</H3><hr/>'))
      }
    })
    
    output$gap <- renderUI({HTML('<br/>')})
    
    # replicateTable
    output$replicateTable <- renderTable({
      # print('replicateTable method')
      sel <- selectProtein(pdat$tab)
      if(length(sel) > 1 && sel > 0 && length(sel) <= max_points){
        data.frame(Sample=colnames(pdat$tab),Intensity=colMeans(pdat$tab[sel,],na.rm = TRUE))
      }
      if (length(sel) == 1 && sel > 0){
        data.frame(Sample=colnames(pdat$tab),Intensity=pdat$tab[sel,])
      }
      if (length(sel) > max_points && sel > 0){ return() }
    },digits = 0, width = "80px"
    )
    
    #significanceTable
    output$significanceTable <- renderTable({
      sel <- selectProtein(pdat$tab)
      if(length(sel) == 1 && sel > 0){
        data.frame(Contrast="1112-BMO", pvalue=as.numeric(DEdat$P.Value[sel]))
      }
    },digits = 4,width = "50px"
    )
    
    #heatMap
    output$heatMap <- renderPlot({
      # print('heatMap method')
      sel<-selectProtein(pdat$tab)
      if(length(sel) > 1 && sel > 0 && length(sel) <= max_points){
        d <- as.matrix(pdat$tab[sel,])
        mean <- rowMeans(d,na.rm = TRUE)
        mean[mean == 0] <- 1
        d <- d/mean
        d[is.nan(d)] <- NA
        row.names(d) <- pdat$annotation$protein[sel]
        # print(d)
        heatmap.2(d, na.rm=TRUE, dendrogram = "row",key=FALSE,keysize = 1,lhei = c(1,100),Colv = FALSE,srtRow = -35,cexRow = 1.0,na.color = "blue")
      }
    })
    
    #jitterPlot
    output$jitterPlot <- renderPlot({
      # print('jitterPlot method')
      sel <- selectProtein(pdat$tab)
      if(length(sel)>0 && sel > 0 && length(sel) <= max_points){
        dataIntensity<-pdat$tab[sel,]
        if (input$intensityScale == 'Log'){
          dataIntensity<-log10(dataIntensity)
        }
        if (length(sel) == 1){
          i <- dataIntensity
        }else{
          i <- colMeans(dataIntensity,na.rm = TRUE)
        }
        s <- sapply(as.data.frame(dataIntensity),function(x) sd(x,na.rm = TRUE)/sqrt(length(x)))
        n <- length(sel)
        p <- data.frame(
          intensity= i,
          lo=i-s,
          up=i+s,
          condition=factor(pdat$metadata$condition,levels=unique(pdat$condition)),
          replicates=factor(pdat$metadata$replicate)
        )
        #shape
        p$shape <- rep(21,length(p$intensity))
        p$shape[which(p$intensity==0)] <- 24
        pd <- position_dodge(width = 0.4)
        #colorblind friendly definition
        cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        ggplot(p, aes(x=condition, y=intensity, ymin=lo, ymax=up, colour=replicates, shape= shape, fill=replicates)) +
          theme(text = element_text(size=20), legend.position = "bottom",legend.direction = "horizontal") +
          {if (input$intensityScale == '') ylim(0, NA)} +
          {if(n > 1) geom_errorbar(position=pd, width = 0.1)} +
          geom_point(position=pd, size=4) +
          scale_shape_identity() +  # necessary for shape mapping
          scale_fill_manual(values=cbPalette) +
          {if (input$intensityScale == 'Log') labs(x = 'Condition', y = 'Log Intensity') else labs(x = 'Condition', y = 'Intensity')}
      }
    })
    
    #FID plot
    output$plotFID <- renderPlot({
      # print('plotFID method')
      tab_idx <- as.numeric(input$allProteinTable_rows_selected)
      pFID <- plotFID(pdat,binhex=FALSE)
      if (length(tab_idx) > 0){
        pFID <- pFID + geom_point(data=FDIdf[tab_idx,],size=3,color='red')
      }
      pFID
    })
    
    #AllProteinTable
    output$allProteinTable <-DT::renderDataTable({
      # print('allProteinTable method')
      d <- data.frame(ProteinId=DEdat$protein,Protein_Name=pdat$annotation$`PROTEIN-NAMES`,mean_1112=formatC(DEdat$mean_1112),mean_BMO=formatC(DEdat$mean_BMO))
      datatable(
        d, class = 'cell-border strip hover'
      ) %>% formatStyle(0, cursor = 'pointer')
    })
    
    #Open browser for Uniprot
    observeEvent(input$allProteinTable_cell_clicked, {
      info = input$allProteinTable_cell_clicked
      # do nothing if not clicked yet, or the clicked cell is not in the 1st column
      if (is.null(info$value) || info$col != 1) return()
      browseURL(DEdat$urls[info$row])
    })
  }
  
  # Run the application 
  shinyApp(ui = ui, server = server)
}