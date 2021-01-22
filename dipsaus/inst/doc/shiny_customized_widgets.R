## ---- include = FALSE, echo=FALSE---------------------------------------------

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup, echo=FALSE--------------------------------------------------------
library(shiny)
library(dipsaus)

## ---- results='asis', echo=FALSE----------------------------------------------
shiny::includeCSS(system.file('www/shared/bootstrap/css/bootstrap.min.css', 
                              package = 'shiny'))
shiny::includeCSS(system.file('rmarkdown/templates/html_vignette/resources/vignette.css', 
                              package = 'rmarkdown'))
shiny::includeCSS(system.file('shiny-addons/dipsaus/dipsaus-dipterix-lib.js', 
                              package = 'dipsaus'))

## ---- eval=FALSE--------------------------------------------------------------
#  
#  # UI function
#  actionButtonStyled(inputId, label, icon = NULL, width = NULL,
#                     btn_type = "button", type = "primary", class = "", ...)
#  
#  # Update function
#  updateActionButtonStyled(session, inputId, label = NULL, icon = NULL,
#                           type = NULL, disabled = NULL, ...)

## ---- results='asis', echo=FALSE----------------------------------------------

btypes <- c('default', 'primary', 'info', 'success', 'warning', 'danger')
btypes <- rbind(btypes, stringr::str_to_title(btypes))

tags$table(
  class = "table table-bordered text-center",
  tags$tbody(
    tags$tr(tagList(
      apply(btypes, 2, function(x){
        tags$th(x[2], br(), tags$small(paste0('type="', x[1], '"')))
      }),
      tags$th('')
    )),
    tags$tr(
      apply(btypes, 2, function(x){
        tags$td(actionButtonStyled('btn', x[2], type = x[1], btn_type = 'a', width = '100%'))
      }),
      tags$td('[default]')
    ),
    tags$tr(
      apply(btypes, 2, function(x){
        tags$td(actionButtonStyled('btn', x[2], type = x[1], btn_type = 'a', 
                                   width = '100%', class = 'btn-lg'))
      }),
      tags$td('class="btn-lg"')
    ),
    tags$tr(
      apply(btypes, 2, function(x){
        tags$td(actionButtonStyled('btn', x[2], type = x[1], btn_type = 'a', 
                                   width = '100%', class = 'btn-sm'))
      }),
      tags$td('class="btn-sm"')
    ),
    tags$tr(
      apply(btypes, 2, function(x){
        tags$td(actionButtonStyled('btn', x[2], type = x[1], btn_type = 'a', 
                                   width = '100%', class = 'btn-xs'))
      }),
      tags$td('class="btn-xs"')
    ),
    tags$tr(
      apply(btypes, 2, function(x){
        tags$td(actionButtonStyled('btn', x[2], type = x[1], btn_type = 'a', 
                                   width = '100%', disabled=TRUE))
      }),
      tags$td('disabled=TRUE')
    )
  )
)



## ---- eval=FALSE--------------------------------------------------------------
#  compoundInput2(
#    'compound', 'Group Label', label_color = 1:10,
#    components = div(
#      textInput('txt', 'Text'),
#      selectInput('sel', 'Select', choices = 1:10, multiple = TRUE),
#      sliderInput('sli', 'Slider', max=1, min=0, val=0.5)
#    ), max_ncomp = 10, min_ncomp = 0, initial_ncomp = 1
#  )

## ---- eval = FALSE------------------------------------------------------------
#  # Bad example
#  observeEvent(input$A, {
#    updateSliderInput(session, 'B', value = input$A)
#  })
#  observeEvent(input$B, {
#    updateTextInput(session, 'A', value = input$B)
#  })

## ---- eval = FALSE------------------------------------------------------------
#  sync_shiny_inputs(input, session, inputIds = c('A', 'B'), uniform = list(
#    function(a){as.numeric(a)},
#    function(b){ b }
#  ), updates = list(
#    function(a){updateTextInput(session, 'A', value = a)},
#    function(b){updateSliderInput(session, 'B', value = b)}
#  ))

