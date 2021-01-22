####################################################################################################################################
### Filename:    GUI.R
### Description: Functions to provide a graphical user interface;
###              'hrm.GUI' is only dealing with the input provided by the user and then calls the function 'hrm_test'
###
###
####################################################################################################################################


#' Function for displaying and saving the results of the function 'hrm_test' within the function 'hrm.GUI'
#'
#' @param result data.frame from the ouput of the function hrm_test
#' @param factors list containing the column names of the factors, first elements are the wholeplot factors, second the subplot factors
#' @param dec decimal mark
#' @param sep seperation mark for saving
#' @keywords internal
gui.results <- function(result, factors, dec, sep) {

  result[,2:6] <- format(result[,2:6], digits = 2, decimal.mark = dec)
  result[ ,6] <- ifelse (result[, 6]<0.001 , gsub(" ", "", paste("<0", dec, "001"), fixed = TRUE),result[, 6])


  # window for results
  windowR <- RGtk2::gtkWindow()
  windowR["title"] <- "HRM"

  quit_cb <- function(widget, window){
    window$destroy()
  }
  save_LaTeX_cb <- function(widget, window){
    directory <- NULL
    tryCatch(directory <- tclvalue(tkchooseDirectory(initialdir=getwd())), error = function(e){  "" }, warning = function(w) "")
    if(!is.null(directory) & !is.na(directory) & directory != ""){
      tryCatch({
        fileConn <- file(paste(directory,"\\result_HRM_analysis.tex",sep=""))
        writeLines(print(xtable(result), include.rownames = FALSE), fileConn)
        close(fileConn)
      }, error = function(e) { GUI_error(e, "An error occured while saving the results:")})
    }
  }
  save_cb <- function(widget, window){
    directory <- NULL
    tryCatch(directory <- tclvalue(tkchooseDirectory(initialdir=getwd())), error = function(e){  "" }, warning = function(w) "")
    if(!is.null(directory) & !is.na(directory) & directory != ""){
      tryCatch({
        write.table(print(result, row.names = FALSE), file = paste(directory,"\\result_HRM_analysis.txt",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
      }, error = function(e) { GUI_error(e, "An error occured while saving the results:")})
    }
  }
  actions <- list(
    list("FileMenu", NULL, "_File"),
    #  list("Open", "gtk-open", "_Open File", "<control>O", "Open CSV", quit_cb),
    list("Save", "gtk-save", "_Save as LaTeX Table", "<control>S", "Save Table", save_LaTeX_cb),
    list("Save2", "gtk-save", "_Save File", "<control>S", "Save Results", save_cb),
    list("Exit", "gtk-quit", "E_xit", "<control>X", "Exit", quit_cb)
  )
  action_group <- RGtk2::gtkActionGroup("spreadsheetActions")
  action_group$addActions(actions, windowR)

  uiManager <- RGtk2::gtkUIManager()
  uiManager$insertActionGroup(action_group, 0)
  merge <- uiManager$newMergeId()
  # File Menu
  uiManager$addUi(merge.id = merge, path = "/", name = "menubar",
                  action = NULL, type = "menubar", top = FALSE)
  uiManager$addUi(merge, "/menubar", "file", "FileMenu", "menu", FALSE)
  #uiManager$addUi(merge, "/menubar/file", "open", "Open", "menuitem", FALSE)
  uiManager$addUi(merge, "/menubar/file", "save", "Save", "menuitem", FALSE)
  uiManager$addUi(merge, "/menubar/file", "save2", "Save2", "menuitem", FALSE)
  uiManager$addUi(merge, "/menubar/file", NULL, NULL, "separator", FALSE)
  uiManager$addUi(merge, "/menubar/file", "exit", "Exit", "menuitem", FALSE)

  # TooLbar
  uiManager$addUi(merge, "/", "toolbar", NULL, "toolbar", FALSE)
  #uiManager$addUi(merge, "/toolbar", "open", "Open", "toolitem", FALSE)
  uiManager$addUi(merge, "/toolbar", "save", "Save", "toolitem", FALSE)
  uiManager$addUi(merge, "/toolbar", "save2", "Save2", "toolitem", FALSE)
  uiManager$addUi(merge, "/toolbar", "exit", "Exit", "toolitem", FALSE)

  menubar <- uiManager$getWidget("/menubar")
  toolbar <- uiManager$getWidget("/toolbar")
  windowR$addAccelGroup(uiManager$getAccelGroup())

  vBox <- RGtk2::gtkVBox()
  windowR$add(vBox)
  vBox$packStart(menubar, expand = FALSE, fill = FALSE, padding = 0)
  vBox$packStart(toolbar, FALSE, FALSE, 0)

  hbox <- RGtk2::gtkHBoxNew(homogeneous = FALSE, spacing = 0)
  #windowR$add(hbox)
  vBox$packStart(hbox)
  vboxLoad <- RGtk2::gtkVBoxNew(homogeneous = FALSE, spacing = 0)
  vboxLoad$setSizeRequest(800,500)
  hbox$add(vboxLoad)

  scroll <- RGtk2::gtkScrolledWindow()
  vbox2 <- RGtk2::gtkVBoxNew(homogeneous = FALSE, spacing = 0)

  scroll$addWithViewport(vbox2)
  vboxLoad$add(scroll)

  frameR <- RGtk2::gtkFrameNew("Results")
  vbox2$add(frameR)

  vBoxR <- RGtk2::gtkVBoxNew()
  vBoxR$setBorderWidth(20)
  frameR$add(vBoxR)   #add vBox to the frame

  hBoxR0 <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBoxR0$setBorderWidth(20)
  vBoxR$packStart(hBoxR0, F, F, 0)

  l <- RGtk2::gtkLabelNew(paste("wholeplot-factors: ", paste0(factors[[1]], sep=", ", collapse = "") ))
  hBoxR0$packStart(l)

  l <- RGtk2::gtkLabelNew(paste("subplot-factors: ", paste0(factors[[2]], sep=", ", collapse = "") ))
  hBoxR0$packStart(l)

  hBoxR <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBoxR$setBorderWidth(20)
  vBoxR$packStart(hBoxR, F, F, 0)
  for(i in 1:7){
    l <- RGtk2::gtkLabelNew(colnames(result)[i])
    RGtk2::gtkLabelSetWidthChars(l, 11)
    if(i == 1){
      RGtk2::gtkLabelSetWidthChars(l, 20)
    }
    hBoxR$packStart(l, F, F, 0)

  }

  # first Row
  for(i in 1:dim(result)[1]){
    hBoxR <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
    hBoxR$setBorderWidth(2)
    vBoxR$packStart(hBoxR, F, F, 0)
    l <- RGtk2::gtkLabelNew(result[i,1])
    RGtk2::gtkLabelSetWidthChars(l, 20)
    hBoxR$packStart(l, F, F, 0)
    RGtk2::gtkLabelSetLineWrap(l, TRUE) # automatic line breaks, if text ist too long

    t <- RGtk2::gtkEntryNew()
    t$setWidthChars(10)
    RGtk2::gtkEntrySetText(t, result[i,2])
    hBoxR$packStart(t, F, F, 0)

    t <- RGtk2::gtkEntryNew()
    t$setWidthChars(10)
    RGtk2::gtkEntrySetText(t, result[i,3])
    hBoxR$packStart(t, F, F, 0)

    t <- RGtk2::gtkEntryNew()
    t$setWidthChars(10)
    RGtk2::gtkEntrySetText(t, result[i,4])
    hBoxR$packStart(t, F, F, 0)

    t <- RGtk2::gtkEntryNew()
    t$setWidthChars(10)
    RGtk2::gtkEntrySetText(t, result[i,5])
    hBoxR$packStart(t, F, F, 0)

    t <- RGtk2::gtkEntryNew()
    t$setWidthChars(10)
    RGtk2::gtkEntrySetText(t, result[i,6])
    hBoxR$packStart(t, F, F, 0)

    t <- RGtk2::gtkLabelNew()
    t$setWidthChars(10)
    RGtk2::gtkLabelSetText(t, result[i,7])
    hBoxR$packStart(t, F, F, 0)
  }

  hBoxR <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBoxR$setBorderWidth(20)
  vBoxR$packStart(hBoxR, F, F, 0)
  t <- RGtk2::gtkLabelNew()
  RGtk2::gtkLabelSetText(t, paste("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"))
  hBoxR$packStart(t, F, F, 0)


  # Functions for Handling the events to save the results
  # save results as LaTeX Code
  saveResults <- function(object, user.data){
    result <- user.data
    directory <- NULL
    tryCatch(directory <- tclvalue(tkchooseDirectory(initialdir=getwd())), error = function(e){  "" }, warning = function(w) "")
    if(!is.null(directory) & !is.na(directory) & directory != ""){
      tryCatch({
        fileConn <- file(paste(directory,"\\result_HRM_analysis.tex",sep=""))
        writeLines(print(xtable(result), include.rownames = FALSE), fileConn)
        close(fileConn)
      }, error = function(e) { GUI_error(e, "An error occured while saving the results:")})
    }
  }
  # save results as plain text
  saveResults2 <- function(object, user.data){
    result <- user.data
    directory <- NULL
    tryCatch(directory <- tclvalue(tkchooseDirectory(initialdir=getwd())), error = function(e){  "" }, warning = function(w) "")
    if(!is.null(directory) & !is.na(directory) & directory != ""){
      tryCatch({
        write.table(print(result, row.names = FALSE), file = paste(directory,"\\result_HRM_analysis.txt",sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
      }, error = function(e) { GUI_error(e, "An error occured while saving the results:")})
    }
  }


  hBoxR <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBoxR$setBorderWidth(20)
  vBoxR$packStart(hBoxR, F, F, 0)


  saveButton <- RGtk2::gtkButton("Save Results as LaTeX Code")
  hBoxR$packStart(saveButton, T, T, 0)

  saveButton2 <- RGtk2::gtkButton("Save Results as Plain Text")
  hBoxR$packStart(saveButton2, T, T, 0)

  closeButton <- RGtk2::gtkButton("Close")
  hBoxR$packStart(closeButton, T, T, 0)


  # Event Handling
  RGtk2:: gSignalConnect(closeButton, "clicked", windowR$destroy)
  RGtk2:: gSignalConnect(saveButton, "clicked", saveResults, data = result)
  RGtk2:: gSignalConnect(saveButton2, "clicked", saveResults2, data = result)


}

#' Graphical User Interface for Testing Multi-Factor High-Dimensional Repeated Measures
#'
#' @description Graphical User Interface (R Package RGtk2 needed) for the Function 'hrm_test': Test for main effects and interaction effects of one or two between-subject factors and one, two or three within-subject factors (at most four factors can be used).
#' @return The results can be saved as LaTeX Code or as plain text. Additionally a plot of the group profiles an be saved when using one whole- and one subplot factor.
#' @keywords export
hrm_GUI <- function(){

  # variable to temporarily store the loaded data from the user
  tmp <- NULL

  # loading the RGtk2 package; stops if it cannot be loaded.
  requireNamespace("RGtk2", quietly = TRUE)
  if(!("package:RGtk2" %in% search())){attachNamespace("RGtk2")}
  if(!isNamespaceLoaded("RGtk2")){
    stop("The package 'RGkt2' is needed for using the graphical user interface.\nPlease install this package or use the function 'hrm_test' from package 'HRM' to perform the analysis without GUI.")
  }
  # loading the cairoDevice package; stops if it cannot be loaded.
  requireNamespace("cairoDevice", quietly = TRUE)
  if(!("package:cairoDevice" %in% search())){attachNamespace("cairoDevice")}
  if(!isNamespaceLoaded("cairoDevice")){
    stop("The package 'cairoDevice' is needed for using the graphical user interface.\nPlease install this package or use the function 'hrm_test' from package 'HRM' to perform the analysis without GUI.")
  }
  # loading the RGtk2Extras package; stops if it cannot be loaded.
  # requireNamespace("RGtk2Extras", quietly = TRUE)
  # if(!("package:RGtk2Extras" %in% search())){attachNamespace("RGtk2Extras")}
  # if(!isNamespaceLoaded("RGtk2Extras")){
  #   stop("The package 'RGtk2Extras' is needed for using the graphical user interface.\nPlease install this package or use the function 'hrm_test' from package 'HRM' to perform the analysis without GUI.")
  # }

  # Functions for Menubar
  quit_cb <- function(widget, window){
    window$destroy()
  }
  open_cb <- function(widget, window){
    getDirectory(NULL, NULL)
  }
  results_cb <- function(widget, window){
    calculate(NULL, NULL)
  }
  # view_cb <- function(widget, window){
  #   RGtk2Extras::dfview(tmp)
  # }
  #list("View Data", "gtk-open", "_View Data", "<control>V", "View CSV", view_cb),

  window <- RGtk2::gtkWindow()
  window["title"] <- "HRM"

  actions <- list(
    list("FileMenu", NULL, "_File"),
    list("Open", "gtk-open", "_Open File", "<control>O", "Open CSV", open_cb),
    list("Get Results", "gtk-open", "Get _Results", "<control>R", "Get Results", results_cb),
    list("Exit", "gtk-quit", "E_xit", "<control>X", "Exit", quit_cb)
  )
  action_group <- RGtk2::gtkActionGroup("spreadsheetActions")
  action_group$addActions(actions, window)

  uiManager <- RGtk2::gtkUIManager()
  uiManager$insertActionGroup(action_group, 0)
  merge <- uiManager$newMergeId()
  # File Menu
  uiManager$addUi(merge.id = merge, path = "/", name = "menubar",
                  action = NULL, type = "menubar", top = FALSE)
  uiManager$addUi(merge, "/menubar", "file", "FileMenu", "menu", FALSE)
  uiManager$addUi(merge, "/menubar/file", "open", "Open", "menuitem", FALSE)
  #uiManager$addUi(merge, "/menubar/file", "view", "View Data", "menuitem", FALSE)
  uiManager$addUi(merge, "/menubar/file", "Calculate", "Get Results", "menuitem", FALSE)
  uiManager$addUi(merge, "/menubar/file", NULL, NULL, "separator", FALSE)
  uiManager$addUi(merge, "/menubar/file", "exit", "Exit", "menuitem", FALSE)

  # TooLbar
  uiManager$addUi(merge, "/", "toolbar", NULL, "toolbar", FALSE)
  uiManager$addUi(merge, "/toolbar", "open", "Open", "toolitem", FALSE)
  #uiManager$addUi(merge, "/toolbar", "view", "View Data", "toolitem", FALSE)
  uiManager$addUi(merge, "/toolbar", "Calculate", "Get Results", "toolitem", FALSE)
  uiManager$addUi(merge, "/toolbar", "exit", "Exit", "toolitem", FALSE)

  menubar <- uiManager$getWidget("/menubar")
  toolbar <- uiManager$getWidget("/toolbar")
  window$addAccelGroup(uiManager$getAccelGroup())

  vBoxM <- RGtk2::gtkVBox()
  window$add(vBoxM)
  vBoxM$packStart(menubar, expand = FALSE, fill = FALSE, padding = 0)
  vBoxM$packStart(toolbar, FALSE, FALSE, 0)

  frame <- RGtk2::gtkFrameNew("Loading Data")
  #window$add(frame)
  vBoxM$packStart(frame)
  # Creating one vertical Box, which consists of multiple horizontal Boxes

  vBox <- RGtk2::gtkVBoxNew()
  vBox$setBorderWidth(20)
  frame$add(vBox)   #add vBox to the frame

  # first Row

  hBox <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBox$setBorderWidth(20)
  vBox$packStart(hBox, F, F, 0)


  labelPath <- RGtk2::gtkLabelNewWithMnemonic("File-Path:")
  hBox$packStart(labelPath, F, F, 0)

  pathEntry <- RGtk2::gtkEntryNew()
  pathEntry$setWidthChars(60)
  hBox$packStart(pathEntry, F, F, 0)

  loadButton <- RGtk2::gtkButton("Load Data")
  hBox$packStart(loadButton, F, F, 10)

  labelPath2 <- RGtk2::gtkLabelNewWithMnemonic("Data needs to be in a long table format,\ni.e. all measurements have to be in one column.")
  hBox$packStart(labelPath2, F, F, 0)

  # additional for first row

  hBox1 <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBox1$setBorderWidth(20)
  vBox$packStart(hBox1, F, F, 0)


  labelSep <- RGtk2::gtkLabelNewWithMnemonic("Seperator:")
  hBox1$packStart(labelSep, F, F, 0)

  sepEntry <- RGtk2::gtkEntryNew()
  sepEntry$setWidthChars(1)
  hBox1$packStart(sepEntry, F, F, 0)
  RGtk2::gtkEntrySetText(sepEntry,",")

  labelDec <- RGtk2::gtkLabelNewWithMnemonic("Decimal:")
  hBox1$packStart(labelDec, F, F, 0)

  decEntry <- RGtk2::gtkEntryNew()
  decEntry$setWidthChars(1)
  hBox1$packStart(decEntry, F, F, 0)
  RGtk2::gtkEntrySetText(decEntry,".")

  labelHeader <- RGtk2::gtkLabelNewWithMnemonic("Header:")
  hBox1$packStart(labelHeader, F, F, 0)

  headerCheck <- RGtk2::gtkCheckButtonNew()
  hBox1$packStart(headerCheck, F, F, 0)
  headerCheck$active <- TRUE

  # second Row

  hBox2 <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBox2$setBorderWidth(20)
  vBox$packStart(hBox2, F, F, 0)

  labelFormula <- RGtk2::gtkLabelNewWithMnemonic("Formula:")
  hBox2$packStart(labelFormula, F, F, 0)

  formulaEntry <- RGtk2::gtkEntryNew()
  formulaEntry$setWidthChars(60)
  hBox2$packStart(formulaEntry, F, F, 0)

  labelFormula2 <- RGtk2::gtkLabelNewWithMnemonic("e.g. measurement ~ groupfactor * timefactor")
  hBox2$packStart(labelFormula2, F, F, 0)

  # third Row

  hBox3 <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBox3$setBorderWidth(20)
  vBox$packStart(hBox3, F, F, 0)

  labelColumns <- RGtk2::gtkLabelNewWithMnemonic("Columns:")
  hBox3$packStart(labelColumns, F, F, 0)

  model<-RGtk2::rGtkDataFrame(c("No Data Loaded"))
  columnsCombo <- RGtk2::gtkComboBox() #text label
  crt <- RGtk2::gtkCellRendererText()
  columnsCombo$packStart(crt)
  columnsCombo$addAttribute(crt, "text", 0)

  RGtk2::gtkComboBoxSetActive(columnsCombo,0)
  hBox3$packStart(columnsCombo, F, F, 0)

  explButton <- RGtk2::gtkButton("Explained by (~)")
  hBox3$packStart(explButton, T, T, 0)

  plusButton <- RGtk2::gtkButton("Additive (+)")
  hBox3$packStart(plusButton, T, T, 0)

  interactionButton <- RGtk2::gtkButton("Interaction Only (:)")
  hBox3$packStart(interactionButton, T, T, 0)

  fullButton <- RGtk2::gtkButton("Additive and Interaction (*)")
  hBox3$packStart(fullButton, T, T, 0)

  # fourth Row

  hBox4 <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBox4$setBorderWidth(20)
  vBox$packStart(hBox4, F, F, 0)

  labelSubject <- RGtk2::gtkLabelNewWithMnemonic("Subject Column:")
  hBox4$packStart(labelSubject, F, F, 0)

  subjectEntry <- RGtk2::gtkEntryNew()
  subjectEntry$setWidthChars(52)

  modelS<-RGtk2::rGtkDataFrame(c("No Data Loaded"))
  columnsSCombo <- RGtk2::gtkComboBox()
  crt2 <- RGtk2::gtkCellRendererText()
  columnsSCombo$packStart(crt2)
  columnsSCombo$addAttribute(crt2, "text", 0)

  # labelSubject2 = RGtk2::gtkLabelNewWithMnemonic("Choose a Column:")
  # hBox4$packStart(labelSubject2, F, F, 0)

  RGtk2::gtkComboBoxSetActive(columnsSCombo,0)
  hBox4$packStart(columnsSCombo, F, F, 0)

  # fifth Row

  hBox5 <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBox5$setBorderWidth(20)
  vBox$packStart(hBox5, F, F, 0)

  labelType1 <- RGtk2::gtkLabelNewWithMnemonic("Type-I Error Rate:")
  hBox5$packStart(labelType1, F, F, 0)

  type1Entry <- RGtk2::gtkEntryNew()
  type1Entry$setWidthChars(52)
  hBox5$packStart(type1Entry, F, F, 0)
  RGtk2::gtkEntrySetText(type1Entry, "0.05")

  # last Row

  hBox6 <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBox6$setBorderWidth(20)
  vBox$packStart(hBox6, F, F, 0)

  okButton <- RGtk2::gtkButton("OK")
  hBox6$packStart(okButton, T, T, 0)

  closeButton <- RGtk2::gtkButton("Close")
  hBox6$packStart(closeButton, T, T, 0)



  # function for loading the data file
  getDirectory <- function(object, user.data){
    directory <- NULL
    # get the path to the file
    tryCatch({directory <- tclvalue(tkgetOpenFile(initialdir=getwd()))
    stopifnot(directory!="")

    RGtk2::gtkEntrySetText(pathEntry, directory)

    sep <- RGtk2::gtkEntryGetText(sepEntry)
    dec <- RGtk2::gtkEntryGetText(decEntry)
    header <- headerCheck$active}, error = function(e) "", warning = function(w) "")

    # read the file given there is a valid path
    if(!is.null(directory) & directory != ""){
      tryCatch({tmp <<- read.table(directory, sep=sep, dec=dec, header=header)
      model<-RGtk2::rGtkDataFrame(c("choose", colnames(tmp)))
      modelS<-RGtk2::rGtkDataFrame(c("choose", colnames(tmp)))
      RGtk2::gtkComboBoxSetModel(columnsCombo, model)
      RGtk2::gtkComboBoxSetModel(columnsSCombo, modelS)}, error = function(e) { GUI_error(e, "There was a problem to load the data:")} )

      RGtk2::gtkEntrySetText(formulaEntry, "")
    }
  }

  setFormulaCombo <- function(object, user.data){
    column <- colnames(tmp)[RGtk2::gtkComboBoxGetActive(object)]
    RGtk2::gtkEditableInsertText(formulaEntry, paste(column,""), position = RGtk2::gtkEntryGetTextLength(formulaEntry))
  }

  setFormulaComboS <- function(object, user.data){
    column <- colnames(tmp)[RGtk2::gtkComboBoxGetActive(object)]
    RGtk2::gtkEntrySetText(subjectEntry, gsub(" ", "", paste(column, ""), fixed = TRUE))
  }

  setFormulaPlus <- function(object, user.data){
    RGtk2::gtkEditableInsertText(formulaEntry, paste("+",""), position = RGtk2::gtkEntryGetTextLength(formulaEntry))
  }

  setFormulaInteraction <- function(object, user.data){
    RGtk2::gtkEditableInsertText(formulaEntry, paste(":",""), position = RGtk2::gtkEntryGetTextLength(formulaEntry))
  }

  setFormulaFull <- function(object, user.data){
    RGtk2::gtkEditableInsertText(formulaEntry, paste("*",""), position = RGtk2::gtkEntryGetTextLength(formulaEntry))
  }

  setFormulaExpl <- function(object, user.data){
    RGtk2::gtkEditableInsertText(formulaEntry, paste("~",""), position = RGtk2::gtkEntryGetTextLength(formulaEntry))
  }

  calculate <- function(object, user.data){

    errorOccured <- 0
    alpha <- NA

    # verifying the input for the formula, subject and alpha error and calculation of the test statistics
    if(is.character(RGtk2::gtkEntryGetText(formulaEntry)) & is.character(RGtk2::gtkEntryGetText(subjectEntry))) {
      if(nchar(RGtk2::gtkEntryGetText(formulaEntry))>0 & nchar(RGtk2::gtkEntryGetText(subjectEntry))>0){

        tryCatch(formula <- as.formula(RGtk2::gtkEntryGetText(formulaEntry)), error = function(e) {GUI_error(e, "There is a problem with your formula:")
          errorOccured <<- 1} )

        if(errorOccured == 0){
          factors <- attributes(terms.formula(formula))$term.labels

          tryCatch({
            subject <- RGtk2::gtkEntryGetText(subjectEntry)
            tmp[,subject] <- as.factor(tmp[,subject])
          }, error = function(e) {GUI_error(e, "Please check the column name for the subject:")
            errorOccured <<- 1} )

          tmpSubset <- subset(tmp, tmp$subject == tmp$subject[1])
          groupFactor <- NULL
          timeFactor <- NULL
          nfactors <- 0 # to count, how many factors there are used
          for(i in 1:length(factors)){
            factorsSplit <- unlist(strsplit(factors[i], ":"))
            if(length(factorsSplit)==1 & errorOccured == 0){
              nfactors <- nfactors + 1
              tryCatch({
                tmp[,factorsSplit] <- as.factor(tmp[,factorsSplit])
                if(nlevels(tmp[,factorsSplit]) == dim(tmpSubset)[1] ){
                  timeFactor <- factorsSplit
                } else if(nlevels(tmp[,factorsSplit]) < dim(tmpSubset)[1] ){
                  groupFactor <- factorsSplit
                }
              }, error = function(e) {GUI_error(e, "There is a problem with your formula for the explaining variables.")
                errorOccured <<- 1} )
            }
          }
          tryCatch(alpha <- as.double(RGtk2::gtkEntryGetText(type1Entry)), warning = function(w) { GUI_error(w, "The Type-I error rate needs to be numeric.")}, error = function(e) { GUI_error(e, "The Type-I error rate needs to be numeric.")})
          tryCatch(responseVariable <- as.character(terms.formula(formula)[[2]]), warning = function(w) "", error = function(e) "" )

          tryCatch({
            if(!is.numeric(tmp[,responseVariable])){
              GUI_error(NULL, "The response variable needs to be numeric.")
              errorOccured <- 1
            }
          }, warning = function(w) "", error = function(e) {GUI_error(e, "There is a problem with your formula for the response variable.")
            errorOccured <<- 1})
        }


        if(!is.na(alpha) & errorOccured == 0){
          if(alpha > 0 & alpha < 1) {
            # if the input by the user is fine, then do the caluclation
            tryCatch({

                object.hrm <-  hrm_test(formula = formula, data = tmp, alpha = alpha, subject = subject )
                result <- object.hrm$result

                # determin which columns are whole- and subplot factors
                dat <- model.frame(formula, tmp)
                dat2 <- data.frame(dat,subject=tmp[,subject])
                m <- ncol(dat)
                # find out, in which columns are the wholeplot or subplot factors
                s1<-subset(dat2, dat2$subject==dat2$subject[1])
                measurements <- dim(s1)[1]
                countSubplotFactor <- 1
                wholeplot<-rep(-1, m)
                subplot<-rep(-1, m)
                for(i in 2:m){
                  if(length(unique(s1[,i]))==nlevels(dat2[,i])){
                    subplot[i]<-1
                    countSubplotFactor <- countSubplotFactor*nlevels(s1[,i])
                  }
                  else{
                    wholeplot[i]<-1
                  }
                }
                wholeplot <- which(wholeplot==1)
                subplot <- which( subplot==1)
                factors <- list(colnames(dat)[wholeplot], colnames(dat)[subplot])

                # showing results
                gui.results(result, factors, RGtk2::gtkEntryGetText(decEntry), RGtk2::gtkEntryGetText(sepEntry))
              }, error = function(e) {GUI_error(e, NULL)
                errorOccured <<- 1})

            # if there are only two factors; 1 whole- and 1 subplot-facor, then plot the profiles
            if(nfactors == 2 & !is.null(groupFactor) & !is.null(timeFactor) & errorOccured == 0){
              tryCatch(responseVariable <- as.character(terms.formula(formula)[[2]]), warning = function(w) "", error = function(e) "" )
              if(is.character(responseVariable)){
                print("Profiles are being plotted ...")
                GUI_plot()
                tryCatch({
                  print(plot.HRM(object.hrm))
                }, error = function(e) "", warning = function(w)  "")
              }
            } else if(nfactors == 1 & !is.null(timeFactor) & errorOccured == 0) {
                GUI_plot()
                tryCatch({
                  print(plot.HRM(object.hrm))
                }, error = function(e) "", warning = function(w)  "")
            }
          } else {
            GUI_error(NULL, "The type-I error rate needs to be within the interval (0, 1).")
          }
        }
        if(RGtk2::gtkEntryGetText(type1Entry) == ""){
          GUI_error(NULL, "The type-I error rate is missing.")
        }

      } else {
        GUI_error(NULL, "Formula or column for the subject is missing.")
      }
    } else {
      GUI_error(NULL, "Please check your input for the formula or the subject column.")
    }


  }


  RGtk2:: gSignalConnect(loadButton, "clicked", getDirectory)
  RGtk2:: gSignalConnect(columnsCombo, "changed", setFormulaCombo)
  RGtk2:: gSignalConnect(columnsSCombo, "changed", setFormulaComboS)
  RGtk2:: gSignalConnect(plusButton, "clicked", setFormulaPlus)
  RGtk2:: gSignalConnect(interactionButton, "clicked", setFormulaInteraction)
  RGtk2:: gSignalConnect(fullButton, "clicked", setFormulaFull)
  RGtk2:: gSignalConnect(explButton, "clicked", setFormulaExpl)
  RGtk2:: gSignalConnect(okButton, "clicked", calculate)
  RGtk2:: gSignalConnect(closeButton, "clicked", window$destroy)

}

#' Function for presenting error messages to the user
#'
#' @param e exception
#' @param msg message to be presented to the end user
#' @keywords internal
GUI_error <- function(e, msg){
  windowE <- RGtk2::gtkWindow()
  windowE["title"] <- "HRM"

  frameE <- RGtk2::gtkFrameNew("Error")
  windowE$add(frameE)

  vBoxE <- RGtk2::gtkVBoxNew()
  vBoxE$setBorderWidth(20)
  frameE$add(vBoxE)   #add vBox to the frame

  hBoxE <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBoxE$setBorderWidth(20)
  vBoxE$packStart(hBoxE, F, F, 0)

  hBoxE$packStart(RGtk2::gtkLabelNewWithMnemonic(msg), F, F, 0)

  hBoxE2 <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBoxE2$setBorderWidth(20)
  vBoxE$packStart(hBoxE2, F, F, 0)
  hBoxE2$packStart(RGtk2::gtkLabelNewWithMnemonic(e$message), F, F, 0)

  hBoxE3 <- RGtk2::gtkHBoxNew(spacing= 10) #distance between elements
  hBoxE3$setBorderWidth(20)
  vBoxE$packStart(hBoxE3, F, F, 0)

  closeButton <- RGtk2::gtkButton("OK")
  hBoxE3$packStart(closeButton, T, T, 0)
  RGtk2:: gSignalConnect(closeButton, "clicked", windowE$destroy)

}

GUI_plot <- function(){
  win <- RGtk2::gtkWindow(show = FALSE)
  win["title"] <- "HRM"
  win$setDefaultSize(400, 400)

  quit_cb <- function(widget, window){
    window$destroy()
  }
  save_cb <- function(widget, window){
    tryCatch(directory <- tclvalue(tkchooseDirectory(initialdir=getwd())), error = function(e){  "" }, warning = function(w) "")
    if(!is.null(directory) & !is.na(directory) & directory != ""){
      tryCatch({ggsave(path = "C:\\Users\\b1011921\\Desktop\\", filename = "Plot_HRM.pdf")
      }, error = function(e) { GUI_error(e, "An error occured while saving the results:")})
    }
  }
  actions <- list(
    list("FileMenu", NULL, "_File"),
    #  list("Open", "gtk-open", "_Open File", "<control>O", "Open CSV", quit_cb),
    list("Save", "gtk-save", "_Save File", "<control>S", "Save CSV", save_cb),
    list("Exit", "gtk-quit", "E_xit", "<control>X", "Exit", quit_cb)
  )
  action_group <- RGtk2::gtkActionGroup("spreadsheetActions")
  action_group$addActions(actions, win)

  uiManager <- RGtk2::gtkUIManager()
  uiManager$insertActionGroup(action_group, 0)
  merge <- uiManager$newMergeId()
  # File Menu
  uiManager$addUi(merge.id = merge, path = "/", name = "menubar",
                  action = NULL, type = "menubar", top = FALSE)
  uiManager$addUi(merge, "/menubar", "file", "FileMenu", "menu", FALSE)
  #uiManager$addUi(merge, "/menubar/file", "open", "Open", "menuitem", FALSE)
  uiManager$addUi(merge, "/menubar/file", "save", "Save", "menuitem", FALSE)
  uiManager$addUi(merge, "/menubar/file", NULL, NULL, "separator", FALSE)
  uiManager$addUi(merge, "/menubar/file", "exit", "Exit", "menuitem", FALSE)

  # Toobar
  uiManager$addUi(merge, "/", "toolbar", NULL, "toolbar", FALSE)
  #uiManager$addUi(merge, "/toolbar", "open", "Open", "toolitem", FALSE)
  uiManager$addUi(merge, "/toolbar", "save", "Save", "toolitem", FALSE)
  uiManager$addUi(merge, "/toolbar", "exit", "Exit", "toolitem", FALSE)

  menubar <- uiManager$getWidget("/menubar")
  toolbar <- uiManager$getWidget("/toolbar")
  win$addAccelGroup(uiManager$getAccelGroup())

  graphics <- RGtk2::gtkDrawingArea()
  vBox <- RGtk2::gtkVBox()
  vBox$packStart(menubar, expand = FALSE, fill = FALSE, padding = 0)
  vBox$packStart(toolbar, FALSE, FALSE, 0)
  vBox$packStart(graphics, expand = TRUE, fill = TRUE, padding = 0)
  win$add(vBox)
  cairoDevice::asCairoDevice(graphics)
  win$show()
}

