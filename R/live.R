#' @import shiny
#' @import ggplot2

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


#' Protein information
#'
#' Internal function to create a rendered table with protein annotations.
#'
#' @param tab Table used to create a plot
#' @param input Input variable from shiny server
#' @param pdat \code{proteusData} object with data and annotations
#' @param max_points Maximum number of points to select
#'
#' @return Rendered table with protein(s) annotations
proteinInfo <- function(tab, input, pdat, max_points) {
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


#' Replicate table
#'
#' Internal function to create a rendered table with replicate intensities. When
#' multiple points are selected, the returned table contains mean values. If
#' \code{intensityScale} in the input is 'log', data are log-10 transformed.
#'
#' @param tab Table used to create a plot
#' @param input Input variable from shiny server
#' @param pdat \code{proteusData} object with data and annotations
#' @param max_points Maximum number of points to select
#'
#' @return Rendered table with replicate intensities.
replicateTable <- function(tab, input, pdat, max_points) {
  renderTable({
    sel <- selectProtein(tab, input)
    if(!is.null(sel)) {
      dat <- pdat$tab[sel,, drop=FALSE]
      if(input$intensityScale == 'Log'){
        dat <- log10(dat)
      }
      if(length(sel) <= max_points) {
        df <- data.frame(Sample=colnames(dat), Intensity=signif(colMeans(dat, na.rm = TRUE), 3))
        df$Intensity[is.nan(df$Intensity)] <- NA
        df
      }
    }
  }, width = "80px")
}


#' Significance table
#'
#' @param tab Table used to create a plot
#' @param res Result from limma
#' @param input Input variable from shiny server
#'
#' @return A rendered table with p-value and adjusted p-value
significanceTable <- function(tab, res, input) {
  renderTable({
    sel <- selectProtein(tab, input)
    if(!is.null(sel) && length(sel) == 1) {
      data.frame(`P-value`=sprintf("%.2g", res$P.Value[sel]), `adjusted P-value`=sprintf("%.2g", res$adj.P.Val[sel]), check.names = FALSE)
    }
  }, width = "100px")
}


#' Jitter plot
#'
#' Internal function to create a jitter plot with replicate intensities versus
#' condition.
#'
#' @param tab Table used to create a plot
#' @param input Input variable from shiny server
#' @param pdat \code{proteusData} object with data and annotations
#' @param max_points Maximum number of points to select
#'
#' @return A ggplot2 object.
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
      pd <- ggplot2::position_dodge(width = 0.4)
      ggplot(p, aes_(x=~condition, y=~intensity, ymin=~lo, ymax=~up, colour=~replicate, shape=~shape, fill=~replicate)) +
        theme(text = element_text(size=20), legend.position = "none", legend.direction = "horizontal") +
        {if (input$intensityScale == '') ylim(0, NA)} +
        geom_point(position=pd, size=4, na.rm=TRUE) +
        {if(n > 1) geom_errorbar(position=pd, width = 0.1)} +
        scale_shape_identity() +  # necessary for shape mapping
        viridis::scale_fill_viridis(discrete=TRUE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
        {if (input$intensityScale == 'Log') labs(x = 'Condition', y = 'Log Intensity') else labs(x = 'Condition', y = 'Intensity')}
    }
  })

}


#' All protein table
#'
#' Internal function to create a table with all proteins (bottom of the web page).
#'
#' @param res A data frame with results from limmaDE.
#'
#' @return A rendered table.
allProteinTable <- function(res) {
  DT::renderDataTable({
    # assume first column is id ("protein" or "peptide")
    idcol <- names(res)[1]
    cols <- c(idcol, "logFC", "adj.P.Val", grep("mean_", colnames(res), value=TRUE))
    d <- res[, cols]
    d[, 2:ncol(d)] <- sapply(d[, 2:ncol(d)], function(x) signif(x, 3))
    d <- DT::datatable(d, class = 'cell-border strip hover')
    DT::formatStyle(d, 0, cursor = 'pointer')
  })
}



#' Volcano plot in live Shiny session
#'
#' \code{plotVolcano_live} makes a volcano plot from limma results.
#'
#' @param pdat Protein \code{proteusData} object.
#' @param res Result table from  \code{\link{limmaDE}}.
#' @param max_points Maximum number of points that can be selected.
#'
#' @return A \code{shiny} session hosted locally.
#'
#' @examples
#' library(shiny)
#' library(proteusLabelFree)
#' data(proteusLabelFree)
#' prodat.med <- normalizeData(prodat)
#' res <- limmaDE(prodat.med)
#' \dontrun{
#' plotVolcano_live(prodat.med, res)
#' }
#'
#' @export
plotVolcano_live <- function(pdat, res, max_points=100){
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")

  res$"-log10(P.Value)" <- -log10(res$P.Value)

  #######################################################################

  ui <- shinyUI(fluidPage(

    # Application title
    titlePanel("PlotVolcano live"),

    fluidRow(
      column(5, plotOutput("plotVolcano", height = "600px", width = "100%", brush = "plot_brush",hover="plot_hover")),
      column(7,
        fluidRow(tableOutput("proteinInfo")),
          fluidRow(
            column(4,
              radioButtons("intensityScale","Intesity Scale:",choices = c("Linear scale" = "","Log scale"="Log"),inline = TRUE)
              )
             ),
             fluidRow(
               column(4,
                fluidRow(plotOutput("jitterPlot", height = "300px",width = "100%")),
                fluidRow(htmlOutput("gap")),
                fluidRow(tableOutput("significanceTable"))
               ),
               column(3,
                    fluidRow(tableOutput("replicateTable"))
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
    output$proteinInfo <- proteinInfo(res, input, pdat, max_points)
    output$gap <- renderUI({HTML('<br/>')})
    output$replicateTable <- replicateTable(res, input, pdat, max_points)
    output$significanceTable <- significanceTable(res, res, input)
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

    output$allProteinTable <- allProteinTable(res)
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
#' library(shiny)
#' library(proteusLabelFree)
#' data(proteusLabelFree)
#' prodat.med <- normalizeData(prodat)
#' res <- limmaDE(prodat.med,sig.level = 0.05)
#' \dontrun{
#' plotFID_live(prodat.med, res)
#' }
#'
#' @export
plotFID_live <- function(pdat, res, max_points=100){
  if(!is(pdat, "proteusData")) stop ("Input data must be of class proteusData.")
  res$"-log10(P.Value)" <- -log10(res$P.Value)

  # Generate the fold-change/intensity dataset. The same as in the FID plot.
  condMeans <- function(cond) {
    m <- rowMeans(log10(pdat$tab)[,which(pdat$metadata$condition == cond), drop=FALSE], na.rm=TRUE)
    m[is.nan(m)] <- NA
    m
  }
  m1 <- condMeans(pdat$conditions[1])
  m2 <- condMeans(pdat$conditions[2])
  good <- !is.na(m1) & !is.na(m2)
  fi <- data.frame(
    id = rownames(pdat$tab),
    x = (m1 + m2) / 2,
    y = m2 - m1,
    good = good
  )
  rownames(fi) <- 1:nrow(fi)

  mx <- 1.1 * max(abs(fi$y), na.rm=TRUE)
  m <- rbind(m1[!good], m2[!good])
  fi[!good, "x"] <- colSums(m, na.rm=TRUE)
  fi[!good, "y"] <- ifelse(is.na(m[1,]), mx, -mx)

  #######################################################################

  ui <- shinyUI(fluidPage(

    # Application title
    titlePanel("PlotFID live"),

    fluidRow(
      column(5, plotOutput("plotFID", height = "600px", width = "100%", brush = "plot_brush",hover="plot_hover")),
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
                      fluidRow(tableOutput("replicateTable")))
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
    output$proteinInfo <- proteinInfo(fi, input, pdat, max_points)
    output$gap <- renderUI({HTML('<br/>')})
    output$replicateTable <- replicateTable(fi, input, pdat, max_points)
    output$significanceTable <- significanceTable(fi, res, input)
    output$jitterPlot <- jitterPlot(fi, input, pdat, max_points)

    #FID plot
    output$plotFID <- renderPlot({
      tab_idx <- as.numeric(input$allProteinTable_rows_selected)
      pFID <- plotFID(pdat,binhex=FALSE)
      if (length(tab_idx) > 0){
        pFID <- pFID + geom_point(data=fi[tab_idx,], size=3, color='red')
      }
      pFID
    })

    output$allProteinTable <- allProteinTable(res)
  }

  # Run the application
  shinyApp(ui = ui, server = server)
}
