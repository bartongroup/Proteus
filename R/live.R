# Fetch selected proteins from a plot or table
selectProtein <- function(tab, input, max_hover=1) {
  sel <- NULL
  tab_idx <- as.numeric(input$allProteinTable_rows_selected)
  if(!is.null(input$plot_brush)){
    brushed <- na.omit(brushedPoints(tab, input$plot_brush))
    sel <- as.numeric(rownames(brushed))
  } else if(!is.null(input$plot_hover)) {
    near <- nearPoints(tab, input$plot_hover, threshold = 20, maxpoints = max_hover)
    sel <- as.numeric(rownames(near))
  } else if(length(tab_idx) > 0) {
    sel <- tab_idx
  }
  return(sel)
}


#' proteinInfo
#'
#' Create a rendered table with protein annotations.
#'
#' @param tab Table used to create a plot
#' @param input Input variable from shiny server
#' @param pdat \code{proteusData} object with data and annotations
#'
#' @return Rendered table with protein(s) annotations
proteinInfo <- function(tab, input, pdat) {
  renderTable({
    sel <- selectProtein(tab, input)
    if(!is.null(sel)) {
      n <- length(sel)
      if (is.null(pdat$annotation)){
        data.frame(Error='No annotation found on the Proteus object. Consult vignette.')
      } else {
        if (n >= 1 && n <= max_points && sel > 0) {
          data.frame(pdat$annotation[sel, ])
        } else if (n > max_points) {
          data.frame(Error=paste('Only', max_points, 'points can be selected.'))
        }
      }
    }
  })
}


#' replicateTable
#'
#' Create a rendered table with replicate intensities. When multiple points are
#' selected, the returned table contains mean values. If \code{intensityScale}
#' in the input is 'log', data are log-10 transformed.
#'
#' @param tab Table used to create a plot
#' @param input Input variable from shiny server
#' @param pdat \code{proteusData} object with data and annotations
#'
#' @return Rendered table with replicate intensities.
replicateTable <- function(tab, input, pdat) {
  renderTable({
    sel <- selectProtein(tab, input)
    if(!is.null(sel)) {
      dat <- pdat$tab[sel,, drop=FALSE]
      if(input$intensityScale == 'Log'){
        dat <- log10(dat)
      }
      if(length(sel) <= max_points) {
        data.frame(Sample=colnames(dat), Intensity=sprintf("%.3g", colMeans(dat, na.rm = TRUE)))
      }
    }
  }, width = "80px")
}


#' significanceTable
#'
#' @param res Result from limma
#' @param input Input variable from shiny server
#'
#' @return A rendered table with pvalue and adjusted p-value
significanceTable <- function(res, input) {
  renderTable({
    sel <- selectProtein(res, input)
    if(!is.null(sel) && length(sel) == 1) {
      data.frame(`P-value`=sprintf("%.2g", res$P.Value[sel]), `adjusted P-value`=sprintf("%.2g", res$adj.P.Val[sel]), check.names = FALSE)
    }
  }, width = "100px")
}


jitterPlot <- function(tab, input, pdat, max_points) {
  renderPlot({
    sel <- selectProtein(tab, input)
    if(!is.null(sel) && length(sel) <= max_points) {
      dat <- pdat$tab[sel,, drop=FALSE]
      if(input$intensityScale == 'Log'){
        dat <- log10(dat)
      }
      m <- colMeans(dat, na.rm = TRUE)
      s <- apply(dat, 2, function(x) sd(x, na.rm = TRUE) / sqrt(na.omit(length(x))))
      n <- length(sel)
      p <- data.frame(
        intensity = m,
        lo = m - s,
        up = m + s,
        condition = factor(pdat$metadata$condition, levels=pdat$conditions),
        replicate = as.factor(pdat$metadata$replicate)
      )
      p$shape <- rep(21, length(p$intensity))
      p$shape[which(p$intensity==0)] <- 24
      pd <- position_dodge(width = 0.4)
      # colorblind friendly definition (alas, only 7 colours, so doesn't work for more replicates)
      # cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
      ggplot(p, aes(x=condition, y=intensity, ymin=lo, ymax=up, colour=replicate, shape=shape, fill=replicate)) +
        theme(text = element_text(size=20), legend.position = "bottom", legend.direction = "horizontal") +
        {if (input$intensityScale == '') ylim(0, NA)} +
        geom_point(position=pd, size=4) +
        {if(n > 1) geom_errorbar(position=pd, width = 0.1)} +
        scale_shape_identity() +  # necessary for shape mapping
        #scale_fill_manual(values=cbPalette) +
        {if (input$intensityScale == 'Log') labs(x = 'Condition', y = 'Log Intensity') else labs(x = 'Condition', y = 'Intensity')}
    }
  })

}



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
plotVolcano_live <- function(pdat, res, max_points=100){
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")

  res$"-log10(P.Value)" <- -log10(res$P.Value)

  #######################################################################

  ui <- shinyUI(fluidPage(

    # Application title
    titlePanel("PlotVolcano live version"),

    fluidRow(
      column(5, plotOutput("plotVolcano", height = "700px", width = "100%", brush = "plot_brush",hover="plot_hover")),
      column(7,
             fluidRow(tableOutput("proteinInfo")),
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
               column(3,
                      fluidRow(tableOutput("replicateTable"))
               )
#               column(5,
#                      fluidRow(plotOutput("heatMap",height = "700px" ,width = "100%"))
#               )
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
    output$proteinInfo <- proteinInfo(res, input, pdat)
    output$gap <- renderUI({HTML('<br/>')})
    output$replicateTable <- replicateTable(res, input, pdat)
    output$significanceTable <- significanceTable(res, input)
    output$jitterPlot <- jitterPlot(res, input, pdat, max_points)

    # Volcano plot
    output$plotVolcano <- renderPlot({
      tab_idx <- as.numeric(input$allProteinTable_rows_selected)
      pVol <- plotVolcano(res, binhex=FALSE)
      if(length(tab_idx) > 0) {
        pVol <- pVol + geom_point(data=res[tab_idx,], size=3, color='red')
      }
      pVol
    })

    # AllProteinTable
    output$allProteinTable <- DT::renderDataTable({
      cols <- c("protein", "logFC", "adj.P.Val", grep("mean_", colnames(res), value=TRUE))
      d <- res[, cols]
      d[, 2:ncol(d)] <- sapply(d[, 2:ncol(d)], function(x) sprintf("%.3g", x))
      datatable(
        d,
        class = 'cell-border strip hover'
      ) %>% formatStyle(0, cursor = 'pointer')
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
plotFID_live <- function(pdat, res, max_points=100){
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  res$"-log10(P.Value)" <- -log10(res$P.Value)

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
             fluidRow(tableOutput("proteinInfo")),
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
               column(3,
                      fluidRow(tableOutput("replicateTable"))),
               column(5,
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
      sel = -1
      tab_idx <- as.numeric(input$allProteinTable_rows_selected)
      if(!is.null(input$plot_brush)){
        brushed <- na.omit(brushedPoints(FDIdf, input$plot_brush))
        sel <- as.numeric(rownames(brushed))
      }else if(!is.null(input$plot_hover)){
        near <- nearPoints(FDIdf,input$plot_hover,threshold = 20,maxpoints = max_hover)
        sel <- as.numeric(rownames(near))
      }else if(length(tab_idx)>0){
        sel <- tab_idx
      }
      return(sel)
    }

    #ProteinInfo
    output$proteinInfo <- renderTable({
      sel <- selectProtein(pdat$tab)
      n <- length(sel)
      if (!'annotation' %in% names(pdat)){
        data.frame(Error='no annotation found on the Proteus object. Consult vignette.')
      }else{
        if (n >= 1 && n <= max_points && sel > 0){
          data.frame(pdat$annotation[sel,])
        }else if (n > max_points){
          data.frame(Error=paste0('only ',max_points,' points can be selected.'))
        }
      }
    })

    output$gap <- renderUI({HTML('<br/>')})

    # replicateTable
    output$replicateTable <- renderTable({
      sel <- selectProtein(pdat$tab)
      if(length(sel) > 1 && sel > 0 && length(sel) <= max_points){
        data.frame(Sample=colnames(pdat$tab),Intensity=colMeans(pdat$tab[sel,],na.rm = TRUE))
      }
      if (length(sel) == 1 && sel > 0){
        data.frame(Sample=colnames(pdat$tab),Intensity=pdat$tab[sel,])
      }
    },digits = 0, width = "80px"
    )

    #significanceTable
    output$significanceTable <- renderTable({
      sel <- selectProtein(pdat$tab)
      if(length(sel) == 1 && sel > 0){
        data.frame(Contrast="1112-BMO", pvalue=as.numeric(res$P.Value[sel]))
      }
    },digits = 4,width = "50px"
    )

    #heatMap
    output$heatMap <- renderPlot({
      sel<-selectProtein(pdat$tab)
      if(length(sel) > 1 && sel > 0 && length(sel) <= max_points){
        d <- as.matrix(pdat$tab[sel,])
        mean <- rowMeans(d,na.rm = TRUE)
        mean[mean == 0] <- 1
        d <- d/mean
        d[is.nan(d)] <- NA
        row.names(d) <- rownames(pdat$tab)[sel]
        heatmap.2(d, na.rm=TRUE, dendrogram = "row",key=FALSE,keysize = 1,lhei = c(1,100),Colv = FALSE,srtRow = -35,cexRow = 1.0,na.color = "blue")
      }
    })

    #jitterPlot
    output$jitterPlot <- renderPlot({
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
      tab_idx <- as.numeric(input$allProteinTable_rows_selected)
      pFID <- plotFID(pdat,binhex=FALSE)
      if (length(tab_idx) > 0){
        pFID <- pFID + geom_point(data=FDIdf[tab_idx,],size=3,color='red')
      }
      pFID
    })

    #AllProteinTable
    output$allProteinTable <-DT::renderDataTable({
      d <- data.frame(ProteinId=res$protein,mean_1112=formatC(res$mean_1112),mean_BMO=formatC(res$mean_BMO))
      datatable(
        d, class = 'cell-border strip hover'
      ) %>% formatStyle(0, cursor = 'pointer')
    })
  }

  # Run the application
  shinyApp(ui = ui, server = server)
}
