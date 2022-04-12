library(tidyverse)
library(dplyr)
library(corrr)
library(caret)
library(umap)
library(cluster)
library(factoextra)
library(randomForest)
library(EBImage)
library(e1071)
library(RColorBrewer)
library(plotly)
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(pryr)
library(DT)


images_folder <- "~/Desktop/Sergi_master_thesis/Aim3_large_screen_analysis/3.2_data_preparation/3.2.1_images_classifiers/single_cell_images/"

list_images_files <- list.files(images_folder)

load("objects/numeric_features.RData")

#Gene annotation

query_genes <- readRDS("objects/query_genes.rds") %>% separate(plate_barcode, into = c("barcode", "screen"), sep = "_")

target_genes <- readRDS("objects/target_genes.rds")

#I need the segmentation classifier to classifiy the single cells and only show those which are well segmented. 

filtering_classifier <- readRDS("objects/filtering_classifier.rds")

google_green <- "#3cba54"

google_red <- "#db3236"

colours <- c(google_green, google_red)


ui <- fluidPage(theme = "journal",
                fluidRow(
                  
                  #We will add an attribute to each cell called context There will be two groups: "isolated" and "crowded".
                  #We will create buttons. When the user clicks one of them, the displayed cells will be assigned to the selected group. 
                  #A skip button will discard the shown cell. This could be useful for occasions in which the user is not sure about the class to which the shown cell
                  #belongs. By doing this, we hope to decrease the error of the classifier. 
                  
                  
                  column(12, offset = 2, actionGroupButtons(inputIds = c("Crowded","Isolated", "skip"),
                                                            label =c("Crowded","Isolated", "Skip")))),
                  
                  fluidRow(
                  
                  #radio buttons with the group names will also be added. When a user selects one of the groups the images of the first 40 classified cells from the 
                  #group will be displayed in the group images tab.
                  
                  
                  column(12,offset = 1, radioButtons(inputId = "group",
                                                     label = "Local_context",
                                                     choices = list("Crowded","Isolated"),
                                                     inline = TRUE))),
                
                mainPanel(tabsetPanel(
                  
                  #In the main panel we would like to have in the main tab the images and the buttons. In the main tab we will show the segmented cropped cells on the first
                  #row and the cell in the raw image on the 2nd row. 
                  #Sliders to adjust the brightness of the image will also be placed. 
                  ##The screen, plate, number, well. will be displayed for the shown image.

                  
                  tabPanel("Images", fluidRow(column(10, displayOutput("image_field")),
                                              column(2,
                                                sliderInput("dapi_adjust", "Dapi:", min = 0.1, max = 10, value = 2, step = 0.01),
                                                sliderInput("actin_adjust", "Actin:", min = 0.1, max = 10, value = 1.5, step = 0.01),
                                                sliderInput("tubulin_adjust", "Tubulin:", min = 0.1, max = 10, value = 1, step = 0.01),
                                                tableOutput("features"))),
                           
                  #Below the single cell images I will show an image covering a bigger field of view. This will allow us to determine better whether the envirnoment
                  #is crowded or isolated. 
                  #After the first 100 cells the classifier will be trained. In order to make the algorithm more accurate 
                  #we could show the probabilities that the model assigns to each of the groups to facilitate the classification. 
                           
                                     fluidRow(column(10, displayOutput("images")),
                                              column(2, tableOutput("class_probabilities")))),
                  
                  
                  
                  
                  #The second tab contains the dataframe with the labelled cells so far. 
                  
                  tabPanel("Classified cells", DT::dataTableOutput("table")),
                  
                  #A third tab will show for the selected group 50 cells in the dapi, actin and tubulin channel. This could be useful to determine the quality of the 
                  #classification and how biologically meaningful is it. 
                  
                  tabPanel("Group images", fluidRow(column(12,displayOutput("dapi_group"))),
                                           fluidRow(column(12,displayOutput("actin_group"))),
                                           fluidRow(column(12,displayOutput("tubulin_group")))),
                  
                  #A foruth tab will show how many cells have been classified for each group.
                  
                  tabPanel("Group numbers",fluidRow(column(12,DT::dataTableOutput("number_cells")))),
                  
                  #A fifth tab will show the the efficiency of the Random forest by plotting the classification error as the number of classified cells increases. 
                  
                  tabPanel("Classification error", plotOutput("model_error")),
                  
                  #It is intersting to see which classes are classified more accurately by the model. This can be accomplished by showing the confusion matrix.
                  #It will also be good to show which features are important to build the model, therefore in one of the tabs the plot of feature importance will 
                  #be shown.
                  
                  tabPanel("Confusion matrix", plotlyOutput("heat_map")),
                  
                  #I will show the distribution of the cells in a UMAP and we will fill by the group.
                  
                  tabPanel("UMAP", plotlyOutput("umap")))
              
              ))

server <- function(input, output, session) {
  
  #A data table will be shown with the user classified cells. It has to be a reactive object in order to bind the rows corresponding to the classified cells. The table
  # will be referred throughout the app as classified_cells$df. 
  #Remaining cells will be an object which has a series of numbers which ranks from 1 to the number of total images. Every time an image is shown, the corresponding number
  #(which correspond to the position of the file in the list_images_file object) will be eliminated. By doing this we avoid showing the same image twice. 
  #As a double checking step we will put all used numbers into the object used_numbers. If the cell is repeated we will print a string in the app. 
  
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
  
  #CROWDED GROUP
  
  observeEvent(input$Crowded,{

    cell_row <- readRDS(paste0(images_folder, list_images_files[image_number$vector]))$features %>% as.data.frame()
    
    cell_row["good_segm"] <- "Yes"

    cell_row["context"] <- "Crowded"
    
    if(colnames(cell_row)[1] == "barcode"){
      
      cell_row <- cell_row %>% select(-barcode)
      
    }

    cell_row <- cell_row %>% dplyr::select(screen, well, plate, field, good_segm, context, everything()) 
    
    classified_cells$df <- classified_cells$df %>% rbind(cell_row)
    
    classified_cells$used_numbers <- c(classified_cells$used_numbers, image_number$vector)
    
    #The array corresponding to the images will be append to the group entry of the lists. Each image has to have a number and this is given
    #by the term [[length(group_images_dapi$list[[input$group]])+1]]. So that it will start with 1 and go all the way to 40. It plays a similar role as if we were doing a for
    #loop and include a count variable.
    
    #In this case it makes more sense to save the image of the bigger field, since it shows much better the local context of the cell. 


    if(length(group_images_dapi$list[["Crowded"]])<40){

      group_images_dapi$list[["Crowded"]][[length(group_images_dapi$list[["Crowded"]])+1]]<-  rgbImage(blue=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi_field[,,3]))

      group_images_actin$list[["Crowded"]][[length(group_images_dapi$list[["Crowded"]])+1]] <-  rgbImage(red=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin_field[,,1]))

      group_images_tubulin$list[["Crowded"]][[length(group_images_dapi$list[["Crowded"]])+1]] <- rgbImage(green=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin_field[,,2]))

    }
    
    #The next step will be to save the images. One way to do it is once a group has reached 40 members, the images are saved. They are saved as RGB images. Each channel image
    #is saved in a different tif file called dapi, actin, tubulin. The image is only saved one time when a group reaches 40 cells.

    if(nrow(classified_cells$df %>% filter(context == "Crowded"))==40){

      writeImage(rgbImage(blue = tile(EBImage::combine(group_images_dapi$list[["Crowded"]]), fg.col = "white"),
                          red = tile(EBImage::combine(group_images_actin$list[["Crowded"]]), fg.col = "white"),
                          green = tile(EBImage::combine(group_images_tubulin$list[["Crowded"]]), fg.col = "white")),
                 files = c("images_groups/crowded/actin.tif",
                           "images_groups/crowded/tubulin.tif",
                           "images_groups/crowded/dapi.tif"))

    }

   })
  
  
  #################################################################################################################################################################
  
  # ISOLATED GROUP
  
  observeEvent(input$Isolated,{
    
    cell_row <- readRDS(paste0(images_folder, list_images_files[image_number$vector]))$features %>% as.data.frame()
    
    cell_row["good_segm"] <- "Yes"
    
    cell_row["context"] <- "Isolated"
    
    if(colnames(cell_row)[1] == "barcode"){
      
      cell_row <- cell_row %>% select(-barcode)
      
    }
    
    cell_row <- cell_row %>% select(screen, well, plate, field, good_segm, context, everything()) 
    
    classified_cells$df <- classified_cells$df %>% rbind(cell_row)
    
    classified_cells$used_numbers <- c(classified_cells$used_numbers, image_number$vector)
    
    
    #The array corresponding to the images will be append to the group entry of the lists. Each image has to have a number and this is given
    #by the term [[length(group_images_dapi$list[[input$group]])+1]]. So that it will start with 1 and go all the way to 40. It plays a similar role as if we were doing a for
    #loop and include a count variable.

    
    if(length(group_images_dapi$list[["Isolated"]])<40){
      
      group_images_dapi$list[["Isolated"]][[length(group_images_dapi$list[["Isolated"]])+1]]<-  rgbImage(blue=(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi_field[,,3]))
      
      group_images_actin$list[["Isolated"]][[length(group_images_dapi$list[["Isolated"]])+1]] <-  rgbImage(red=(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin_field[,,1]))
      
      group_images_tubulin$list[["Isolated"]][[length(group_images_dapi$list[["Isolated"]])+1]] <- rgbImage(green=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin_field[,,2]))
      
    }
    
    #The next step will be to save the images. One way to do it is once a group has reached 40 members, the images are saved. They are saved as RGB images. Each channel image
    #is saved in a different tif file called dapi, actin, tubulin. The image is only saved one time when a group reaches 40 cells.
    
    if(nrow(classified_cells$df %>% filter(context == "Isolated"))==40){
      
      writeImage(rgbImage(blue = tile(EBImage::combine(group_images_dapi$list[["Isolated"]]), fg.col = "white"),
                          red = tile(EBImage::combine(group_images_actin$list[["Isolated"]]), fg.col = "white"),
                          green = tile(EBImage::combine(group_images_tubulin$list[["Isolated"]]), fg.col = "white")),
                 files = c("images_groups/isolated/actin.tif",
                           "images_groups/isolated/tubulin.tif",
                           "images_groups/isolated/dapi.tif"))
      
    }
    
  })
  
 
  #################################################################################################################################################################
  
  #DATAFRAME
  
  #The dataframe with the classified cells is shown in the second tab. It is updated every time a user classifies a cell, this everytime it clicks one of the buttons. 

  output$table <- DT::renderDataTable({
    
  classified_cells$df
    
  })
  
  #It is important that the dataframe with the manually classified cells is saved becuase later on it'll be used as the training set for the Random Forest. 
  #Since the table will not become very large, it can be saved every time a cell is classified.
  #observeEvent can be reactive to multiple inputs. In this case it will execute the code if one of the buttons is clicked. 
   
  observeEvent(c(input$Crowded,input$Isolated),{

    saveRDS(classified_cells$df, file = "training_set_dataframe/training_set.rds")


  })
  
  #################################################################################################################################################################
  
  #REPRESENTATIVE IMAGES
  
  # It will be interesting to show the first 50 classified cells of each group. The images will be tiled and combined and shown
  # to the user in the three channels. By using the radioButtons a group can be chosen and only cells belonging to this cluster will be shown

  output$dapi_group <- renderDisplay({

      if(length(group_images_dapi$list[[input$group]])>0){

        display(tile(EBImage::combine(group_images_dapi$list[[input$group]]), fg.col = "white"))

      }

    })
    # 
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
  
  
  #################################################################################################################################################################
  
  #GROUP NUMBERS
  
  #In the tab Group numbers, a table with the total number of cells by group will be shown.
  
  observeEvent(c(input$Crowded, input$Isolated),{
    
    number_cells$table[1, "number_cells"] <- nrow(classified_cells$df %>% filter(context == "Crowded"))
    
    number_cells$table[2, "number_cells"] <- nrow(classified_cells$df %>% filter(context == "Isolated"))
    
    number_cells$table[3, "number_cells"] <- nrow(classified_cells$df)
    
    number_cells$table["group"] <- c("Crowded", "Isolated", "Total")
    
    saveRDS(number_cells$table, file = "objects/numbers_table.rds")
    
  })
  
  output$number_cells <- DT::renderDataTable({
    
    number_cells$table
    
  })
  
  
  #################################################################################################################################################################
  
  #AVOID REPEATING CELLS AND SHOWING CLASSIFIER PROBABILITIES
  
  #The aim is that after clicking one of the buttons, a new cell will be displayed, therefore we will updated the image_number value every time the user clicks
  #one of the buttons. In principle several events can be added to the observeEvent function.  
  #When we click a labelling button, the image file number will be saved. Before selecting a new cell to show to the user the program will check that it has not been shown 
  #before (the index number of the file is not present in the used_numbers vector) and that the cell is well segmented (the filtering classifier assigns it to the 
  #well segmented group).
  
  observeEvent(c(input$Crowded,input$Isolated, input$skip),{
    
    if(nrow(classified_cells$df)<2000){
      
      for (number in sample(1:length(list_images_files))){
        
        
        cell_row <- readRDS(paste0(images_folder, list_images_files[number]))$features
        
        
        if(number %in% classified_cells$used_numbers == F & predict(filtering_classifier, cell_row) == "Yes"){
          
          image_number$vector <- number
          
          break
          
        }else{next}
      }    
    }    
    
    
    # #When the first 200 cells have been randomly sampled and classified, the filter for cells which are difficult for the model to classify will be applied. 
    # else{
    #   
    #   for (number in sample(1:length(list_images_files))){
    #     
    #     if(number %in% classified_cells$used_numbers == F){
    #       
    #       #A random cell will be picked. Then the class will be predicted by the classifier. Each class gets a probability and then the cell is assigned to the class with
    #       # the higher probability. 
    #       #We need a for loop in case the classifier predicts the image with high accuracy and then is not shown to the user. 
    #       
    #       cell <- readRDS(paste0(images_folder, list_images_files[number]))$features 
    #       
    #       classifier$probabilities <- predict(classifier$model, cell, type = "prob") %>% as.data.frame()
    #       
    #       #the function order gives indexes to the columns based on the value. We will order the groups so that the groups with higher probability are in the first two places. 
    #       
    #       classifier$probabilities <- classifier$probabilities[,order(classifier$probabilities, decreasing = TRUE)]
    #       
    #       #Only if the difference between the probabilities of the classes with higher chances is lower than 0.4 the cell will be chosen and shown to the user. 
    #       #Otherwise another cell will be sampled and the filtering will be applied. 
    #       
    #       if(classifier$probabilities[1,1]-classifier$probabilities[1,2]<0.4){
    #         
    #         image_number$vector <- number
    #         
    #         #Only if the difference in class probabilities is surpassed the image will be shown. 
    #         
    #         break
    #       }else{next}
    
  })
  
  
  ######################################################################################################################################################################
  
  # #CLASS PROBABILITIES TABLE
  # 
  # #Next to the images we will show for the displayed cell which are the probabilities that the model assigns to each of the classes. This could help improve the accuracy
  # #of the model. For instance if the model assigns 0.5 % of probability to one class and 0.4 % to another the most rational decision will be to decide between the two most
  # # likely groups.
  # 
  # output$class_probabilities <- renderTable({
  # 
  #   if(nrow(classified_cells$df)>200){
  # 
  #     probablilities <- t(classifier$probabilities) %>% as.data.frame() %>% dplyr::rename(probability = 1) %>% rownames_to_column(var = "class")
  # 
  #     probablilities
  #   }
  # })
  # 
  
#######################################################################################################################################################################
  
#DISPLAYED IMAGES
  

  # #The single channel images will be tiled and combined in one single image. 
  # #It is updated every time the user clicks the button classify as the values of image number change. 
  #The first row corresponds to the segmented and cropped image while the second row show the raw image in the context of the cell of interest. 
  
  #To make use of the slider we will just multiply the image by the slider value. 
  
  output$images <- renderDisplay({

    display(tile(EBImage::combine(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi_raw*input$dapi_adjust,
                                  readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin_raw*input$actin_adjust,
                                  readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin_raw*input$tubulin_adjust,
                                  readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi*input$dapi_adjust,
                                  readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin*input$actin_adjust,
                                  readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin*input$tubulin_adjust), nx = 3, fg.col = "white"))

  })
  
  #Below that I will show the field images to observe better the context. 
  
  output$image_field <- renderDisplay({
    
    display(tile(EBImage::combine(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi_field*input$dapi_adjust,
                                  readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin_field*input$actin_adjust,
                                  readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin_field*input$tubulin_adjust), nx = 3, fg.col = "white"))
    
  })
  
  
  
  
 ######################################################################################################################################################################
  
  #DISPLAYED CELL INFORMATION
  
  #We will show the screen, plate, well, field, number and of the cell that is currently shown to the user. this could be useful to associate phenotypes to 
  #specific conditions or time points. 
  #I will also show the query and target genes corresponding to the condition. 
  
  output$features <- renderTable({
    
    info <- readRDS(paste0(images_folder, list_images_files[image_number$vector]))$features %>% select(screen, plate, well, field, number, dist.10.nn) 
    
    if(str_split(info$well, "")[[1]][2] == "0"){
      
      info["well"] <- paste0(str_split(info$well, "")[[1]][1], str_split(info$well, "")[[1]][3])
      
    }
    
    # info <- info %>% mutate(query = query_genes %>% filter(screen == info %>% pull(screen)) %>% pull(query_name) %>% as.character(),
    #                         target = target_genes %>% filter(plate == info %>% pull(plate), well == info %>% pull(well)) %>% pull(gene_symbol) %>% as.character())
    
    info
    
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
  
  
  observeEvent(c(input$Crowded, input$Isolated),{

    if(nrow(classified_cells$df)>10 & (number_cells$table[3,"number_cells"]/50) %%1 ==0){

      training_data <- classified_cells$df %>% as.data.frame() %>% select(context, numeric_features)

      training_data[,"context"] <- as.factor(training_data[,"context"])

      classifier$model <- randomForest(context~., importance = TRUE, data = training_data, ntree = 500)

      row <- as.data.frame(matrix(ncol = 0,nrow=0))
      
      row[1,"number_cells"] <- nrow(classified_cells$df)
      
      row[1,"error"] <- as.vector(1-classifier$model$err.rate[500,1])

      error$table <- error$table %>% rbind(as.data.frame(row))
      
      saveRDS(error$table, file = "objects/classifier_error.rds")
      
    }  
    
    #What is more we would like to save the classifier every time it gets updated and 
    #this process also requires a short but noticeable period of time which increases as the training set is bigger. That is why it might be optimal to save it only every
    #50 or 100 classified cells.
    
    if((number_cells$table[3,"number_cells"]/50) %%1 ==0){
    
    saveRDS(classifier$model, file = "objects/classifier.rds")
      
    saveRDS(classified_cells$used_numbers, file = "objects/used_numbers.rds")
    
    
    }
    
  })
  
  ####################################################################################################################################################################
  
  #SLIDERS back to default values
  
  #In order to make the cells more comparable, the sliders which adjust the brightness of the imgaes could go back to the default values (1) after the image change. 
  #If we do that at least the initial image will be comparable among cells. 
  
  observeEvent(c(input$Crowded, input$Isolated, input$skip),{
    
    updateSliderInput(session, "dapi_adjust", value = 2)
    
    updateSliderInput(session, "actin_adjust", value = 1.5)
    
    updateSliderInput(session, "tubulin_adjust", value = 1)
    
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
         scale_fill_gradient(low = "black", high = "white")+
         xlab("user-defined group")+
         ylab("predicted group")+
         theme_classic()+
         theme(legend.position = "none")+
         ggtitle("Confusion matrix (proportion of cells)")

       ggplotly(heat_map)
     }

   })

  #################################################################################################################################################################
  
  
  #UMAP
  
  #I will show the distribution of cells in a UMAP once the 100 first cells are classified. 
  
  
  output$umap <- renderPlotly({
    
    important_features_rf <- classifier$model$importance %>% as.data.frame() %>% rownames_to_column(var = "feature") %>% arrange(desc(MeanDecreaseAccuracy))
    
    #We will pick the 24 most important features for the model. 
    
    important_features_rf <- important_features_rf[1:10,1]
    
    if(nrow(classified_cells$df)>60){
      
      set.seed(100)
      umap_data <- classified_cells$df %>% select(important_features_rf) %>% uwot::umap(n_neighbors = 15) %>% as.data.frame() %>% 
                dplyr::rename(dimension1 = V1, dimension2 = V2) %>% mutate(context = classified_cells$df %>% pull(context))
      
      umap_plot <- ggplot(umap_data, aes(x = dimension1, y = dimension2, colour = context)) + 
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
 }

shinyApp(ui, server)
