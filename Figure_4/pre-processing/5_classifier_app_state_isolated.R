library(tidyverse)
library(dplyr)
library(corrr)
library(caret)
library(cluster)
library(factoextra)
library(EBImage)
library(umap)
library(randomForest)
library(e1071)
library(RColorBrewer)
library(plotly)
library(data.table)
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(pryr)
library(DT)


# In order to increase the speed at which the classifier learns(this is that it increase its accuracy with less training images), after around 150-200 cells, the images
# shown to the user will filtered, so that cells which are difficult for the model to predict will be enriched. By doing this we hope to correct for class imbalance, 
# and also make the learning process of the model much faster. This means that the first 150-200 cells will be drawn randomly, and afterwards a prior filtering will
#be applied. 

images_folder <- "single_cell_isolated_images/"

list_images_files <- list.files(images_folder)

load("objects/numeric_features.RData")

#Gene annotation

query_genes <- readRDS("objects/query_genes.rds") %>% separate(plate_barcode, into = c("barcode", "screen"), sep = "_")

target_genes <- readRDS("objects/target_genes.rds") %>% mutate(gene_symbol = as.character(gene_symbol), plate = as.character(plate), well = as.character(well))

#Colours UMAP

google_green <- "#3cba54"

google_blue <- "#4885ed"

google_yellow <- "#f4c20d"

google_red <- "#db3236"

lacoste_green <- "#004526"

colours <- c(google_green, lacoste_green, google_yellow, google_red, google_blue)


ui <- fluidPage(theme = "journal",
                fluidRow(
                  
                  #I will initially define 6 groups into which the cells will be classified:
                  #Normal cells: they show the characteristic round shape and regular nucleus and body size. 
                  #Irregular_nucleus: they can be detected either by looking at the DNA or at the tubulin. 
                  #Elongated: so called "stretchy" cells, they are much longer than normal cells. Generally the tubulin signal is very pronounced. 
                  #Big: Cells which show a body or nucleus which is much bigger than the average. 
                  #Condensed: cells which have a very tiny nucleus and small body will be included in this group. 
                  #Skip: cells which are wrongly segmented or are very difficult to label could be skipped. 
                  
                  column(12, offset = 2, actionGroupButtons(inputIds = c("Normal", "Irregular_nucleus", "Elongated", "Big", "Condensed", "Skip"),
                                                            label =c("Normal", "Irregular nucleus", "Elongated","Big", "Condensed", "Skip")))),
                  
                  fluidRow(
                  
                  #radio buttons with the group names will also be added. When a user selects one of the groups the images of the first 40 classified cells from the 
                  #group will be displayed in the group images tab.
                  #Another set of radio buttons will be shown which allow the user to choose whether the cells shown are completely random or if only cells which the 
                  #classifier does not label with high confidence are shown. The latter option can only be done after the firs 100 cells are classified. 
                  
                  
                  column(8, offset = 1, radioButtons(inputId = "group",
                                                     label = "Cell state",
                                                     choices = list("Big", "Condensed", "Elongated", "Irregular_nucleus", "Normal"),
                                                     inline = TRUE)),
                  column(4, offset = 2, radioButtons(inputId = "cell_selection",
                                         label = "Cell type",
                                         choices = list("Random", "Difficult"),
                                         inline = T))),
                  
                
                #In the main panel we would like to have in the main tab the three images with the buttons, and then perhaps another tab with a table 
                #containing the information of the classified cells (e.g. features, cluster...). We will show the non-normalized images on the top of the tab
                #and the normalized ones at the bottom. the raw image allows comparisons between intensities among different cells. This could make it easier for instance
                #to distinguish an interphase cell from an S-phase cell based on the DNA intensity or distinguish a mitotic cell from one in  G2 by looking at the tubulin
                #intensity. 
                #After the first 100 cells only cells which are difficult for the model to classify are shown. In order to make the algorithm more accurate we could show
                #the probabilities that the model assigns to each of the groups to facilitate the classification. 
                #A third tab will show for the selected group 50 cells in the dapi, actin and tubulin channel. This could be useful to determine the quality of the 
                #classification and how biologically meaningful is it. 
                #A foruth tab will show how many cells have been classified for each group.
                #A fifth tab will show the the efficiency of the Random forest by plotting the classification error as the number of classified cells increases. 
                #It is intersting to see which classes are classified more accurately by the model. This can be accomplished by showing the confusion matrix.
                #It will also be good to show which features are important to build the model, therefore in one of the tabs the plot of feature importance will be shown. 
                
                mainPanel(tabsetPanel(
                  tabPanel("Images", fluidRow(column(10, displayOutput("images")),
                                              column(2,
                                                sliderInput("dapi_adjust", "Dapi:", min = 0.1, max = 10, value = 2, step = 0.01),
                                                sliderInput("actin_adjust", "Actin:", min = 0.1, max = 10, value = 1.5, step = 0.01),
                                                sliderInput("tubulin_adjust", "Tubulin:", min = 0.1, max = 10, value = 1.5, step = 0.01),
                                                tableOutput("features"))),
                           
                           
                           
                #The screen, plate, number, well, dist.10.nn, query name and target name will be displayed for the shown image. 
                #After the first 100 cells only cells which are difficult for the model to classify are shown. In order to make the algorithm more accurate 
                #we could show the probabilities that the model assigns to each of the groups to facilitate the classification. 
                           
                           fluidRow(
                                    #column(2, textOutput("image_number")),
                                    #column(2, textOutput("used_cells")),
                                    column(6, tableOutput("class_probabilities")),
                                    column(6, tableOutput("cell_info"))
                                    )),
                
                #The second tab contains the dataframe with the labelled cells so far. 
                  
                  tabPanel("Classified cells", DT::dataTableOutput("table")),
                
                #A third tab will show for the selected group 40 cells in the dapi, actin and tubulin channel. This could be useful to determine the quality of the 
                #classification and how biologically meaningful is it. 

                  tabPanel("Group images", fluidRow(column(12,displayOutput("dapi_group"))),
                                           fluidRow(column(12,displayOutput("actin_group"))),
                                           fluidRow(column(12,displayOutput("tubulin_group")))),
                
                #A foruth tab will show how many cells have been classified for each group.
                
                  tabPanel("Group numbers",fluidRow(column(12, DT::dataTableOutput("number_cells")))),
                
                #A fifth tab will show the the efficiency of the Random forest by plotting the classification error as the number of classified cells increases.
                
                  tabPanel("Classification error", plotOutput("model_error")),
                
                #It is intersting to see which classes are classified more accurately by the model. This can be accomplished by showing the confusion matrix.
              
                  tabPanel("Confusion matrix", plotlyOutput("heat_map")),
                
                #I will show the distribution of the cells in a UMAP and we will fill by the group.
                
                  tabPanel("UMAP", plotlyOutput("umap")))
                              
              ))

server <- function(input, output, session) {
  
  #A data table will be shown with the user classified cells. It has to be a reactive object in order to bind the rows corresponding to the classified cells. The table
  # will be referred throughout the app as classified_cells$df. 
  #Every time a cell is classified, the file index is save in the used_numbers vector. Before the new cell is loaded it checks whether it has been shown before 
  #(whether its index is contained in the classifier_cells$used_numbers object). This makes sure that every cell is shown only once to the user. 
  
  classified_cells <- reactiveValues(df = as.data.frame(matrix(ncol = 0, nrow = 0)), used_numbers = vector())
  
  #Every time we click a button to label a cell, a new random cell will be selected. We will sample a number and load the corresponding image. This number should
  #be a ReactiveValues object.
  
  image_number <- reactiveValues(vector = sample(1:length(list_images_files), 1))
  
  #We would like to show a sample size of the images of each group. That is why a reactive list will be created. Every time an image is shown, it will be saved 
  #in the corresponding group entry. The images will be shown in the tab Group images as a tiled image. 
  
  group_images_dapi <- reactiveValues(list = list())
  
  group_images_actin <- reactiveValues(list = list())
  
  group_images_tubulin <- reactiveValues(list = list())
  
  # A reactive table will be created which will be updated with the number of cells of each group.
  
  number_cells <- reactiveValues(table = as.data.frame(matrix(ncol = 0, nrow = 0)))
  
  #The classifier has to be a reactive variable because two different chunks of code will use it and update it. After 100 cells the model will predict the class of the 
  #chosen cells. The probabilities of each class will be saved and shown to the user. 
  
  classifier <- reactiveValues(model=NULL, probabilities = NULL)
  
  #We would like to plot the classification error of the model as we increase the number of cells, therefore we will create a reactive dataframe in which we will
  #store the error and the number of cells every iteration
  
  error <- reactiveValues(table = as.data.frame(matrix(ncol = 0, nrow = 0)))
  
 
  
  #Action buttons do not store any information when they are clicked, they can only trigger functions. We would like to know which button has been clicked from the list
  #that is why a diffreren observe event will be created for every button. when the button is clicked in the cell row a group column is added with the name of the group.
  #The row will then be added to the training set dataframe and displayed as data.table in the classified cell tab. 
  
  #################################################################################################################################################################
  
  #NORMAL GROUP
  
  observeEvent(input$Normal,{

    cell_row <- readRDS(paste0(images_folder, list_images_files[image_number$vector]))$features %>% as.data.frame() %>% mutate(field = as.character(field))
    
    if(colnames(cell_row)[1] == "barcode"){
      
      cell_row <- cell_row %>% select(-barcode)
      
    }
    
    cell_row["good_segm"] <- "Yes"
    
    cell_row["context"] <- "Isolated"

    cell_row["group"] <- "Normal"

    cell_row <- cell_row %>% dplyr::select(screen, well, plate, field, number, group, everything()) 
    
    classified_cells$df <- classified_cells$df %>% rbind(cell_row)
    
    classified_cells$used_numbers <- c(classified_cells$used_numbers, image_number$vector)
    
    # if(nrow(classified_cells$df %>% filter(group == "Interphase"))<250){
    # 
    #   classified_cells$df <- classified_cells$df %>% rbind(cell_row)
    # }

    #The array corresponding to the images will be append to the group entry of the lists. Each image has to have a number and this is given
    #by the term [[length(group_images_dapi$list[[input$group]])+1]]. So that it will start with 1 and go all the way to 40. It plays a similar role as if we were doing a for
    #loop and include a count variable.
    #The next step will be to save the images. One way to do it is once a group has reached 40 members, the images are saved. They are saved as RGB images. Each channel image
    #is saved in a different tif file called dapi, actin, tubulin. The image is only saved one time when a group reaches 40 cells.

    if(length(group_images_dapi$list[["Normal"]])<40){

      group_images_dapi$list[["Normal"]][[length(group_images_dapi$list[["Normal"]])+1]]<-  rgbImage(blue=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi[,,3]))

      group_images_actin$list[["Normal"]][[length(group_images_dapi$list[["Normal"]])+1]] <-  rgbImage(red=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin[,,1]))

      group_images_tubulin$list[["Normal"]][[length(group_images_dapi$list[["Normal"]])+1]] <- rgbImage(green=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin[,,2]))

    }

    if(nrow(classified_cells$df %>% filter(group == "Normal"))==40){

      writeImage(rgbImage(blue = tile(EBImage::combine(group_images_dapi$list[["Normal"]]), fg.col = "white"),
                          red = tile(EBImage::combine(group_images_actin$list[["Normal"]]), fg.col = "white"),
                          green = tile(EBImage::combine(group_images_tubulin$list[["Normal"]]), fg.col = "white")),
                 files = c("images_groups/Normal/actin.tif",
                           "images_groups/Normal/tubulin.tif",
                           "images_groups/Normal/dapi.tif"))

    }

   })
  

  #################################################################################################################################################################
  
  # IRREGULAR NUCLEUS GROUP
  
  observeEvent(input$Irregular_nucleus,{
    
    cell_row <- readRDS(paste0(images_folder, list_images_files[image_number$vector]))$features %>% as.data.frame() %>% mutate(field = as.character(field))
    
    if(colnames(cell_row)[1] == "barcode"){
      
      cell_row <- cell_row %>% select(-barcode)
      
    }
    
    cell_row["good_segm"] <- "Yes"
    
    cell_row["context"] <- "Isolated"
    
    cell_row["group"] <- "Irregular_nucleus"
    
    cell_row <- cell_row %>% dplyr::select(screen, well, plate, field, number, group, everything()) 
    
    classified_cells$df <- classified_cells$df %>% rbind(cell_row)
    
    classified_cells$used_numbers <- c(classified_cells$used_numbers, image_number$vector)
    
    
    #The array corresponding to the images will be append to the group entry of the lists. Each image has to have a number and this is given
    #by the term [[length(group_images_dapi$list[[input$group]])+1]]. So that it will start with 1 and go all the way to 40. It plays a similar role as if we were doing a for
    #loop and include a count variable.
    #The next step will be to save the images. One way to do it is once a group has reached 40 members, the images are saved. They are saved as RGB images. Each channel image
    #is saved in a different tif file called dapi, actin, tubulin. The image is only saved one time when a group reaches 40 cells.
    
    if(length(group_images_dapi$list[["Irregular_nucleus"]])<40){
      
      group_images_dapi$list[["Irregular_nucleus"]][[length(group_images_dapi$list[["Irregular_nucleus"]])+1]]<-  rgbImage(blue=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi[,,3]))
      
      group_images_actin$list[["Irregular_nucleus"]][[length(group_images_dapi$list[["Irregular_nucleus"]])+1]] <-  rgbImage(red=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin[,,1]))
      
      group_images_tubulin$list[["Irregular_nucleus"]][[length(group_images_dapi$list[["Irregular_nucleus"]])+1]] <- rgbImage(green=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin[,,2]))
      
    }
    
    if(nrow(classified_cells$df %>% filter(group == "Irregular_nucleus"))==40){
      
      writeImage(rgbImage(blue = tile(EBImage::combine(group_images_dapi$list[["Irregular_nucleus"]]), fg.col = "white"),
                          red = tile(EBImage::combine(group_images_actin$list[["Irregular_nucleus"]]), fg.col = "white"),
                          green = tile(EBImage::combine(group_images_tubulin$list[["Irregular_nucleus"]]), fg.col = "white")),
                 files = c("images_groups/Irregular_nucleus/actin.tif",
                           "images_groups/Irregular_nucleus/tubulin.tif",
                           "images_groups/Irregular_nucleus/dapi.tif"))
      
    }
    
  })
  
  #################################################################################################################################################################
  
  #ELONGATED
  
  observeEvent(input$Elongated,{
    
    cell_row <- readRDS(paste0(images_folder, list_images_files[image_number$vector]))$features %>% as.data.frame() %>% mutate(field = as.character(field))
    
    if(colnames(cell_row)[1] == "barcode"){
      
      cell_row <- cell_row %>% select(-barcode)
      
    }
    
    cell_row["good_segm"] <- "Yes"
    
    cell_row["context"] <- "Isolated"
    
    cell_row["group"] <- "Elongated"
    
    cell_row <- cell_row %>% dplyr::select(screen, well, plate, field, number, group, everything()) 
    
    classified_cells$df <- classified_cells$df %>% rbind(cell_row)
    
    classified_cells$used_numbers <- c(classified_cells$used_numbers, image_number$vector)
    
    
    if(length(group_images_dapi$list[["Elongated"]])<40){
      
      group_images_dapi$list[["Elongated"]][[length(group_images_dapi$list[["Elongated"]])+1]]<-  rgbImage(blue=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi[,,3]))
      
      group_images_actin$list[["Elongated"]][[length(group_images_dapi$list[["Elongated"]])+1]] <-  rgbImage(red=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin[,,1]))
      
      group_images_tubulin$list[["Elongated"]][[length(group_images_dapi$list[["Elongated"]])+1]] <- rgbImage(green=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin[,,2]))
      
    }
    
    if(nrow(classified_cells$df %>% filter(group == "Elongated"))==40){
      
      writeImage(rgbImage(blue = tile(EBImage::combine(group_images_dapi$list[["Elongated"]]), fg.col = "white"),
                          red = tile(EBImage::combine(group_images_actin$list[["Elongated"]]), fg.col = "white"),
                          green = tile(EBImage::combine(group_images_tubulin$list[["Elongated"]]), fg.col = "white")),
                 files = c("images_groups/Elongated/actin.tif",
                           "images_groups/Elongated/tubulin.tif",
                           "images_groups/Elongated/dapi.tif"))
      
    }
    
  })
  
  #################################################################################################################################################################
  
  # BIG GROUP
  
  observeEvent(input$Big,{

    cell_row <- readRDS(paste0(images_folder, list_images_files[image_number$vector]))$features %>% as.data.frame() %>% mutate(field = as.character(field))
    
    if(colnames(cell_row)[1] == "barcode"){
      
      cell_row <- cell_row %>% select(-barcode)
      
    }

    cell_row["good_segm"] <- "Yes"
    
    cell_row["context"] <- "Isolated"
    
    cell_row["group"] <- "Big"

    cell_row <- cell_row %>% dplyr::select(screen, well, plate, field, number, group, everything()) 

    classified_cells$df <- classified_cells$df %>% rbind(cell_row)
    
    classified_cells$used_numbers <- c(classified_cells$used_numbers, image_number$vector)

    if(length(group_images_dapi$list[["Big"]])<40){

      group_images_dapi$list[["Big"]][[length(group_images_dapi$list[["Big"]])+1]]<-  rgbImage(blue=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi[,,3]))

      group_images_actin$list[["Big"]][[length(group_images_dapi$list[["Big"]])+1]] <-  rgbImage(red=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin[,,1]))

      group_images_tubulin$list[["Big"]][[length(group_images_dapi$list[["Big"]])+1]] <- rgbImage(green=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin[,,2]))

    }

    if(nrow(classified_cells$df %>% filter(group == "Big"))==40){

      writeImage(rgbImage(blue = tile(EBImage::combine(group_images_dapi$list[["Big"]]), fg.col = "white"),
                          red = tile(EBImage::combine(group_images_actin$list[["Big"]]), fg.col = "white"),
                          green = tile(EBImage::combine(group_images_tubulin$list[["Big"]]), fg.col = "white")),
                 files = c("images_groups/Big/actin.tif",
                           "images_groups/Big/tubulin.tif",
                           "images_groups/Big/dapi.tif"))

    }

  })
  
  #################################################################################################################################################################
  
  # CONDENSED GROUP
  
  observeEvent(input$Condensed,{
    
    cell_row <- readRDS(paste0(images_folder, list_images_files[image_number$vector]))$features %>% as.data.frame() %>% mutate(field = as.character(field))
    
    if(colnames(cell_row)[1] == "barcode"){
      
      cell_row <- cell_row %>% select(-barcode)
      
    }
    
    cell_row["good_segm"] <- "Yes"
    
    cell_row["context"] <- "Isolated"
    
    cell_row["group"] <- "Condensed"
    
    cell_row <- cell_row %>% dplyr::select(screen, well, plate, field, number, group, everything())  
    
    classified_cells$df <- classified_cells$df %>% rbind(cell_row)
    
    classified_cells$used_numbers <- c(classified_cells$used_numbers, image_number$vector)
    
    
    if(length(group_images_dapi$list[["Condensed"]])<40){
      
      group_images_dapi$list[["Condensed"]][[length(group_images_dapi$list[["Condensed"]])+1]]<-  rgbImage(blue=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi[,,3]))
      
      group_images_actin$list[["Condensed"]][[length(group_images_dapi$list[["Condensed"]])+1]] <-  rgbImage(red=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin[,,1]))
      
      group_images_tubulin$list[["Condensed"]][[length(group_images_dapi$list[["Condensed"]])+1]] <- rgbImage(green=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin[,,2]))
      
    }
    
    if(nrow(classified_cells$df %>% filter(group == "Condensed"))==40){
      
      writeImage(rgbImage(blue = tile(EBImage::combine(group_images_dapi$list[["Condensed"]]), fg.col = "white"),
                          red = tile(EBImage::combine(group_images_actin$list[["Condensed"]]), fg.col = "white"),
                          green = tile(EBImage::combine(group_images_tubulin$list[["Condensed"]]), fg.col = "white")),
                 files = c("images_groups/Condensed/actin.tif",
                           "images_groups/Condensed/tubulin.tif",
                           "images_groups/Condensed/dapi.tif"))
      
    }
    
  })
  
  
#######################################################################################################################################################################  
  
  #SKIP
  
  observeEvent(input$Skip,{
    
    classified_cells$used_numbers <- c(classified_cells$used_numbers, image_number$vector)
    
  })
  
  
  #The dataframe with the classified cells is shown in the second tab. It is updated every time a user classifies a cell, this everytime it clicks one of the buttons. 

  output$table <- DT::renderDataTable({
    
  classified_cells$df 
    
  })
  
  #It is important that the dataframe with the manually classified cells is saved becuase later on it'll be used as the training set for the classifier. 
  #Since the table will not become very large, it can be saved every time a cell is classified.
  #observeEvent can be reactive to multiple inputs. In this case it will execute the code if one of the buttons is clicked. 
   
  observeEvent(c(input$Normal,input$Irregular_nucleus,input$Elongated,input$Big, input$Condensed, input$Skip),{

    saveRDS(classified_cells$df, file = "training_set_dataframe/training_set.rds")
    
    saveRDS(classified_cells$used_numbers, file = "objects/used_numbers.rds")


  })
   


    # It will be interesting to show the first 40 classified cells of each group. The images will be tiled and combined and shown
    # to the user in the three channels. By using the radioButtons a group can be chosen and only cells belonging to this cluster will be shown

  output$dapi_group <- renderDisplay({

      if(length(group_images_dapi$list[[input$group]])>0){

        display(tile(EBImage::combine(group_images_dapi$list[[input$group]]), fg.col = "white"))

      }

    })
 
  output$actin_group <- renderDisplay({

      if(length(group_images_dapi$list[[input$group]])>0){

        display(tile(EBImage::combine(group_images_actin$list[[input$group]]), fg.col = "white"))

      }

  })
  
  output$tubulin_group <- renderDisplay({

      if(length(group_images_dapi$list[[input$group]])>0){

        display(tile(EBImage::combine(group_images_tubulin$list[[input$group]]), fg.col = "white"))

      }

  })
  
  #In the tab Group numbers, a table with the total number of cells by group will be shown.
  
  observeEvent(c(input$Normal,input$Irregular_nucleus,input$Elongated,input$Big, input$Condensed, input$Unknown),{
    
    number_cells$table[1, "number_cells"] <- nrow(classified_cells$df %>% filter(group == "Normal"))
    
    number_cells$table[2, "number_cells"] <- nrow(classified_cells$df %>% filter(group == "Irregular_nucleus"))
    
    number_cells$table[3, "number_cells"] <- nrow(classified_cells$df %>% filter(group == "Elongated"))
    
    number_cells$table[4, "number_cells"] <- nrow(classified_cells$df %>% filter(group == "Big"))
    
    number_cells$table[5, "number_cells"] <- nrow(classified_cells$df %>% filter(group == "Condensed"))
    
    #number_cells$table[6, "number_cells"] <- nrow(classified_cells$df %>% filter(group == "Unknown"))
    
    number_cells$table[6, "number_cells"] <- nrow(classified_cells$df)
    
    number_cells$table["group"] <- c("Normal", "Irregular_nucleus", "Elongated", "Big", "Condensed", "Total")
    
    saveRDS(number_cells$table, file = "objects/numbers_table.rds")
    
  })
  
  output$number_cells <- DT::renderDataTable({
    
    number_cells$table
    
  })
  
  
  #################################################################################################################################################################
  
  #AVOID REPEATING CELLS AND SHOWING CLASSIFIER PROBABILITIES
  
  #The aim is that after clicking one of the buttons, a new cell will be displayed, therefore we will updated the image_number value every time the user clicks
  #one of the buttons. In principle several events can be added to the observeEvent function. 
  #I added action buttons to select whether we want random cells to be selected or cells which the classifier cannot assign to a specific group with high degree
  # of confidence (class probability between the first and second group is lower than 0.4). 
  
  observeEvent(c(input$Normal,input$Irregular_nucleus,input$Elongated,input$Big, input$Condensed, input$Skip),{
    
    if(input$cell_selection == "Random"){
    
      for (number in sample(1:length(list_images_files))){
      
        if(number %in% classified_cells$used_numbers == F){
        
          image_number$vector <- number
          
          break
          
        }else{next}
      }    
    }    
    
    
    
    else{
      
      for (number in sample(1:length(list_images_files))){
        
        if(number %in% classified_cells$used_numbers == F){
          
          #A random cell will be picked. Then the class will be predicted by the classifier. Each class gets a probability and then the cell is assigned to the class with
          # the higher probability. 
          #We need a for loop in case the classifier predicts the image with high accuracy and then is not shown to the user. 
          
            cell <- readRDS(paste0(images_folder, list_images_files[number]))$features 
          
            classifier$probabilities <- predict(classifier$model, cell, type = "prob") %>% as.data.frame()
          
          #the function order gives indexes to the columns based on the value. We will order the groups so that the groups with higher probability are in the first two places. 
          
            classifier$probabilities <- classifier$probabilities[,order(classifier$probabilities, decreasing = TRUE)]
          
          #Only if the difference between the probabilities of the classes with higher chances is lower than 0.4 the cell will be chosen and shown to the user. 
          #Otherwise another cell will be sampled and the filtering will be applied. 
          
            if(classifier$probabilities[1,1]-classifier$probabilities[1,2]<0.4){
            
              image_number$vector <- number
            
              #Only if the difference in class probabilities is surpassed the image will be shown. 
            
              break
            }else{next}
          
        }
      }
      
    }
    
})
  
  
  ######################################################################################################################################################################
  
#DISPLAYED IMAGES
  

  # #The single channel images will be tiled and combined in one single image. 
  # #It is updated every time the user clicks the button classify as the values of image number change. 
  #The first row corresponds to the segmented and cropped image while the second row show the raw image in the context of the cell of interest. 
  
  #To make use of the slider we will just multiply the image by the slider value. 
  
  output$images <- renderDisplay({

    display(tile(EBImage::combine(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi*input$dapi_adjust,
                                  readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin*input$actin_adjust,
                                  readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin*input$tubulin_adjust,
                                  readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi_raw*input$dapi_adjust,
                                  readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin_raw*input$actin_adjust,
                                  readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin_raw*input$tubulin_adjust), nx = 3, fg.col = "white"))

  })
  
  ####################################################################################################################################################################
  
  #SLIDERS back to default values
  
  #In order to make the cells more comparable, the sliders which adjust the brightness of the imgaes could go back to the default values (1) after the image change. 
  #If we do that at least the initial image will be comparable among cells. 
  
  observeEvent(c(input$Normal,input$Irregular_nucleus,input$Elongated,input$Big, input$Condensed, input$Skip),{
    
    updateSliderInput(session, "dapi_adjust", value = 1.75)
    
    updateSliderInput(session, "actin_adjust", value = 1.75)
    
    updateSliderInput(session, "tubulin_adjust", value = 1.75)
    
  })
  
  ######################################################################################################################################################################
  
  #DISPLAYED CELL INFORMATION
  
  #We will show the screen, plate, well, field, number and dist.10.nn of the cell that is currently shown to the user. this could be useful to associate phenotypes to 
  #specific conditions or time points. 
  #I will also show the query and target genes corresponding to the condition. 
  
  output$cell_info <- renderTable({
    
    info <- readRDS(paste0(images_folder, list_images_files[image_number$vector]))$features %>% select(screen, plate, well, field, number, dist.10.nn) 
    
    if(str_split(info$well, "")[[1]][2] == "0"){
      
      info["well"] <- paste0(str_split(info$well, "")[[1]][1], str_split(info$well, "")[[1]][3])
      
    }
    
    info <- info %>% mutate(screen = as.character(screen), plate = as.character(plate), well = as.character(well), field = as.character(field),
                            query = query_genes %>% filter(screen == info %>% pull(screen)) %>% pull(query_name) %>% as.character(),
                            target = target_genes %>% filter(plate == info %>% pull(plate), well == info %>% pull(well)) %>% pull(gene_symbol) %>% as.character())
    
    info 
    
  })
  
  ######################################################################################################################################################################
  
  #CLASS PROBABILITIES TABLE
  
  #It will only be displayed when difficult cells are shown to the user. 
  #Next to the images we will show for the displayed cell which are the probabilities that the model assigns to each of the classes. This could help improve the accuracy
  #of the model. For instance if the model assigns 0.5 % of probability to one class and 0.4 % to another the most rational decision will be to decide between the two most
  # likely groups.
  
  output$class_probabilities <- renderTable({
    
    if(input$cell_selection == "Difficult"){
      
      probablilities <- t(classifier$probabilities) %>% as.data.frame() %>% dplyr::rename(probability = 1) %>% rownames_to_column(var = "class")
      
      probablilities
    }
  })

  
  ######################################################################################################################################################################
  
  #CLASSIFIER
  
  
  #We would like to see how the error rate of the random forest changes as we add more cells to the training set. That is why everytime 50 cells are classified, the random
  #forest will be trained on the dataset. Then the classification error will be extracted and plot against the number of cells in the training set. 
  #We would expect the error of the model to be very high at the beginning, then drop very fast as the number of classified cells increases and eventually form 
  # a plateau. When this level is reach we would know that classifying more cells won't make the model more accurate. It will start to be trained when more than 100 cells
  # have been classified. 
  #To run the Random forest at least 20 cells need to be classified and the group column of the dataframe has to be transformed from character to factor.
  #we could only update the classifier when a certain number of new cells have been added to the training set. This is because when the dataframe becomes larger it takes
  #a little bit of time to calculate the random forest.  that is why perhaps setting it to every 50 cells will be optimal. 
  
  
  observeEvent(c(input$Normal,input$Irregular_nucleus,input$Elongated,input$Big, input$Condensed),{

    if(nrow(classified_cells$df)>10 & ((number_cells$table[6,"number_cells"]/50) %%1) == 0){

      training_data <- classified_cells$df %>% as.data.frame() %>% select(group, numeric_features)

      training_data[,"group"] <- as.factor(training_data[,"group"])

      classifier$model <- randomForest(group~., importance = TRUE, data = training_data)

      row <- as.data.frame(matrix(ncol = 0, nrow=0))

      row[1,"number_cells"] <- nrow(classified_cells$df)

      row[1,"error"] <- as.vector(1-classifier$model$err.rate[500,1])

      error$table <- error$table %>% rbind(as.data.frame(row))

      saveRDS(error$table, file = "objects/classifier_error.rds")

    }

    #What is more we would like to save the classifier every time it gets updated and
    #this process also requires a short but noticeable period of time which increases as the training set is bigger. That is why it might be optimal to save it only every
    #50 or 100 classified cells.

    if((number_cells$table[6,"number_cells"]/50) %%1 ==0){

    saveRDS(classifier$model, file = "objects/classifier.rds")


    }

  })

  ######################################################################################################################################################################
  
  #PLOT OF THE MODEL ACCURACY
  
  output$model_error <- renderPlot({
    
    if(nrow(classified_cells$df)>100){
      
      plot <- print(ggplot(error$table, aes(x = number_cells, y = error)) +
                      geom_line() +
                      ylab("Accuracy (1-out of bag average error)") +
                      xlab("Number of classified cells")+
                      theme_classic()+
                      theme(
                        axis.title.x = element_text(size = 16),
                        axis.title.y = element_text(size = 16)
                      ))
      
    }
    

  }) 
  
  ######################################################################################################################################################################
  
  #CONFUSION MATRIX
  
   output$heat_map <- renderPlotly({

     if(nrow(classified_cells$df)>50){

       data <- classifier$model$confusion %>% as.data.frame() %>% rownames_to_column(var = "group1") %>% select(-class.error) %>%
         gather(key="group2", value="number_cells",-group1) %>% arrange(group1, group2)

       proportion <- vector()

       for(row in 1:nrow(data)){

         group <- data[row,1]

         total_cells <- sum(data %>% filter(group1==group) %>% select(number_cells))

         row_proportion <- data[row, "number_cells"]/total_cells

         proportion <- c(proportion, row_proportion)

       }

       data <- data %>% mutate(proportion_cells=proportion)

       heat_map <- ggplot(data, aes(x=group1, y= group2, fill = proportion_cells))+
         geom_tile()+
         theme_classic()+
         scale_fill_gradient(low = "black", high = "white")+
         xlab("user-defined group")+
         ylab("predicted group")+
         ggtitle("Confusion matrix (proportion of cells)")

       ggplotly(heat_map)
     }

   })
  
  ###################################################################################################################################################################
  
  #UMAP
  
  #I will show the distribution of cells in a UMAP once the 100 first cells are classified. 
  #I will only update this every 50 cells. 
  
  
  output$umap <- renderPlotly({
    
    if(nrow(classified_cells$df)>15){
      
      set.seed(100)
      umap_data <- classified_cells$df %>% select(numeric_features) %>% uwot::umap(n_neighbors = 15) %>% as.data.frame() %>% 
        dplyr::rename(dimension1 = V1, dimension2 = V2) %>% mutate(group = classified_cells$df %>% pull(group))
      
      umap_plot <- ggplot(umap_data, aes(x = dimension1, y = dimension2, colour = group)) + 
        geom_point(size = 2) +
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5, size = 16),
              legend.text = element_text(size = 14),
              legend.position = "bottom",
              legend.title = element_blank(),
              legend.spacing.x = unit(1, "cm"),
              axis.title = element_text(size = 13)) +
        guides(colour = guide_legend(nrow = 1))+
        scale_color_manual(values = colours)+
        guides(colour=guide_legend(override.aes = list(size = 6)))+
        ylab("UMAP dimension 1")+
        xlab("UMAP dimension 2")
      
      ggplotly(umap_plot, height = 700, width = 900)
      
      
    }
    
    
  })


 #   #The importance of the features for the classifier will be shown in one of the tabs
 # 
 #   output$features <- renderPlot({
 # 
 #     varImpPlot(classifier$model)
 #   })
 # 
#   output$image_number <- renderText({
# 
#       image_number$vector
# 
#   })
# 
#   output$used_cells <- renderText({
# 
#       if(image_number$vector %in% classified_cells$used_numbers){
#         
#         "Repeated cell"
#       
#       }else{
#         
#         "New cell"
#         
#       }
# 
#   })
# #  # 
  }

shinyApp(ui, server)
