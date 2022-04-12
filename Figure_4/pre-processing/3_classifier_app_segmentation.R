library(tidyverse)
library(dplyr)
library(corrr)
library(caret)
library(umap)
library(cluster)
library(factoextra)
library(EBImage)
library(umap)
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

query_genes <- readRDS("objects/query_genes.rds")

target_genes <- readRDS("objects/target_genes.rds")

ui <- fluidPage(theme = "journal",
                fluidRow(
                  
                  #We will add an attribute to each cell called good_segm. There will be two groups: "yes" for well-segmented cells and "no" for badly segmented cells
                  #We will create buttons. When the user clicks one of them, the displayed cells will be assigned to the selected group. 
                  #A skip button will discard the shown cell. This could be useful for occasions in which the user is not sure about the class to which the shown cell
                  #belongs. By doing this, we hope to decrease the error of the classifier. 
                  
                  
                  column(12, offset = 2, actionGroupButtons(inputIds = c("Yes","No", "skip"),
                                                            label =c("Yes", "No", "Skip")))),
                  
                  fluidRow(
                  
                  #radio buttons with the group names will also be added. When a user selects one of the groups the images of the first 40 classified cells from the 
                  #group will be displayed in the group images tab.
                  
                  
                  column(12,offset = 1, radioButtons(inputId = "group",
                                                     label = "Well_segmented",
                                                     choices = list("Yes", "No"),
                                                     inline = TRUE))),
                
                mainPanel(tabsetPanel(
                  
                  #In the main panel we would like to have in the main tab the images and the buttons. In the main tab we will show the segmented cropped cells on the first
                  #row and the cell in the raw image on the 2nd row. 

                  
                  tabPanel("Images", fluidRow(column(10, displayOutput("images")),
                                              column(2,
                                                sliderInput("dapi_adjust", "Dapi:", min = 0.1, max = 10, value = 2, step = 0.01),
                                                sliderInput("actin_adjust", "Actin:", min = 0.1, max = 10, value = 1.5, step = 0.01),
                                                sliderInput("tubulin_adjust", "Tubulin:", min = 0.1, max = 10, value = 1, step = 0.01),
                                                tableOutput("features"))),
                           
                  #The screen, plate, number, well and dist.10. will be displayed for the shown image. 
                  #After the first 100 cells only cells which are difficult for the model to classify are shown. In order to make the algorithm more accurate 
                  #we could show the probabilities that the model assigns to each of the groups to facilitate the classification. 
                           
                                     fluidRow(column(6, tableOutput("class_probabilities")),
                                              (column(6, tableOutput("cell_info"))))),
                  
                  #The second tab contains the dataframe with the labelled cells so far. 
                  
                  tabPanel("Classified cells", dataTableOutput("table")),
                  
                  #A third tab will show for the selected group 50 cells in the dapi, actin and tubulin channel. This could be useful to determine the quality of the 
                  #classification and how biologically meaningful is it. 
                  
                  tabPanel("Group images", fluidRow(column(12,displayOutput("dapi_group"))),
                                           fluidRow(column(12,displayOutput("actin_group"))),
                                           fluidRow(column(12,displayOutput("tubulin_group")))),
                  
                  #A foruth tab will show how many cells have been classified for each group.
                  
                  tabPanel("Group numbers",fluidRow(column(12,dataTableOutput("number_cells")))),
                  
                  #A fifth tab will show the the efficiency of the Random forest by plotting the classification error as the number of classified cells increases. 
                  
                  tabPanel("Classification error", plotOutput("model_error")),
                  
                  #It is intersting to see which classes are classified more accurately by the model. This can be accomplished by showing the confusion matrix.
                  #It will also be good to show which features are important to build the model, therefore in one of the tabs the plot of feature importance will 
                  #be shown.
                  
                  tabPanel("Confusion matrix", plotlyOutput("heat_map")))
                  #tabPanel("Feature importance", plotOutput("feature")))
                              
              ))

server <- function(input, output, session) {
  
  #A data table will be shown with the user classified cells. It has to be a reactive object in order to bind the rows corresponding to the classified cells. The table
  # will be referred throughout the app as classified_cells$df. 
  #Remaining cells will be an object which has a series of numbers which ranks from 1 to the number of total images. Every time an image is shown, the corresponding number
  #(which correspond to the position of the file in the list_images_file object) will be eliminated. By doing this we avoid showing the same image twice. 
  #As a double checking step we will put all used numbers into the object used_numbers. If the cell is repeated we will print a string in the app. 
  
  classified_cells <- reactiveValues(df = as.data.frame(matrix(ncol = 0, nrow = 0)), remaining_cells = 1:length(list_images_files), used_numbers = vector())
  
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
  
  #WELL SEGMENTED GROUP
  
  observeEvent(input$Yes,{

    cell_row <- readRDS(paste0(images_folder, list_images_files[image_number$vector]))$features %>% as.data.frame() %>% 
                  mutate(plate = as.character(plate), well = as.character(well), screen = as.character(screen), field = as.character(field))

    cell_row["good_segm"] <- "Yes"
    
    if(colnames(cell_row)[1] == "barcode"){
      
      cell_row <- cell_row %>% select(-barcode)
      
    }

    cell_row <- cell_row %>% dplyr::select(screen, well, plate, field, good_segm, everything()) 
    
    classified_cells$df <- classified_cells$df %>% rbind(cell_row)

    #The array corresponding to the images will be append to the group entry of the lists. Each image has to have a number and this is given
    #by the term [[length(group_images_dapi$list[[input$group]])+1]]. So that it will start with 1 and go all the way to 40. It plays a similar role as if we were doing a for
    #loop and include a count variable.


    if(length(group_images_dapi$list[["Yes"]])<40){

      group_images_dapi$list[["Yes"]][[length(group_images_dapi$list[["Yes"]])+1]]<-  rgbImage(blue=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi[,,3]))

      group_images_actin$list[["Yes"]][[length(group_images_dapi$list[["Yes"]])+1]] <-  rgbImage(red=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin[,,1]))

      group_images_tubulin$list[["Yes"]][[length(group_images_dapi$list[["Yes"]])+1]] <- rgbImage(green=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin[,,2]))

    }
    
    #The next step will be to save the images. One way to do it is once a group has reached 40 members, the images are saved. They are saved as RGB images. Each channel image
    #is saved in a different tif file called dapi, actin, tubulin. The image is only saved one time when a group reaches 40 cells.

    if(nrow(classified_cells$df %>% filter(good_segm == "Yes"))==40){

      writeImage(rgbImage(blue = tile(EBImage::combine(group_images_dapi$list[["Yes"]]), fg.col = "white"),
                          red = tile(EBImage::combine(group_images_actin$list[["Yes"]]), fg.col = "white"),
                          green = tile(EBImage::combine(group_images_tubulin$list[["Yes"]]), fg.col = "white")),
                 files = c("images_groups/well_segm/actin.tif",
                           "images_groups/well_segm/tubulin.tif",
                           "images_groups/well_segm/dapi.tif"))

    }

   })
  
  
  #################################################################################################################################################################
  
  # WRONG-SEGMENTED GROUP
  
  observeEvent(input$No,{
    
    cell_row <- readRDS(paste0(images_folder, list_images_files[image_number$vector]))$features %>% as.data.frame() %>% 
                   mutate(plate = as.character(plate), well = as.character(well), screen = as.character(screen), field = as.character(field))
    
    cell_row["good_segm"] <- "No"
    
    if(colnames(cell_row)[1] == "barcode"){
      
      cell_row <- cell_row %>% select(-barcode)
      
    }
    
    cell_row <- cell_row %>% select(screen, well, plate, field, good_segm, everything()) 
    
    classified_cells$df <- classified_cells$df %>% rbind(cell_row)
    
    
    #The array corresponding to the images will be append to the group entry of the lists. Each image has to have a number and this is given
    #by the term [[length(group_images_dapi$list[[input$group]])+1]]. So that it will start with 1 and go all the way to 40. It plays a similar role as if we were doing a for
    #loop and include a count variable.

    
    if(length(group_images_dapi$list[["No"]])<40){
      
      group_images_dapi$list[["No"]][[length(group_images_dapi$list[["No"]])+1]]<-  rgbImage(blue=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$dapi[,,3]))
      
      group_images_actin$list[["No"]][[length(group_images_dapi$list[["No"]])+1]] <-  rgbImage(red=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$actin[,,1]))
      
      group_images_tubulin$list[["No"]][[length(group_images_dapi$list[["No"]])+1]] <- rgbImage(green=normalize(readRDS(paste0(images_folder, list_images_files[image_number$vector]))$tubulin[,,2]))
      
    }
    
    #The next step will be to save the images. One way to do it is once a group has reached 40 members, the images are saved. They are saved as RGB images. Each channel image
    #is saved in a different tif file called dapi, actin, tubulin. The image is only saved one time when a group reaches 40 cells.
    
    if(nrow(classified_cells$df %>% filter(good_segm == "No"))==40){
      
      writeImage(rgbImage(blue = tile(EBImage::combine(group_images_dapi$list[["No"]]), fg.col = "white"),
                          red = tile(EBImage::combine(group_images_actin$list[["No"]]), fg.col = "white"),
                          green = tile(EBImage::combine(group_images_tubulin$list[["No"]]), fg.col = "white")),
                 files = c("images_groups/bad_segm/actin.tif",
                           "images_groups/bad_segm/tubulin.tif",
                           "images_groups/bad_segm/dapi.tif"))
      
    }
    
  })
  
 
  #################################################################################################################################################################
  
  #DATAFRAME
  
  #The dataframe with the classified cells is shown in the second tab. It is updated every time a user classifies a cell, this everytime it clicks one of the buttons. 

  output$table <- renderDataTable({
    
  classified_cells$df
    
  })
  
  #It is important that the dataframe with the manually classified cells is saved becuase later on it'll be used as the training set for the Random Forest. 
  #Since the table will not become very large, it can be saved every time a cell is classified.
  #observeEvent can be reactive to multiple inputs. In this case it will execute the code if one of the buttons is clicked. 
   
  observeEvent(c(input$Yes,input$No),{

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
  
  observeEvent(c(input$Yes,input$No),{
    
    number_cells$table[1, "number_cells"] <- nrow(classified_cells$df %>% filter(good_segm == "Yes"))
    
    number_cells$table[2, "number_cells"] <- nrow(classified_cells$df %>% filter(good_segm == "No"))
    
    number_cells$table[3, "number_cells"] <- nrow(classified_cells$df)
    
    number_cells$table["group"] <- c("Well_segmented", "Bad_segmented", "Total")
    
    saveRDS(number_cells$table, file = "objects/numbers_table.rds")
    
  })
  
  output$number_cells <- renderDataTable({
    
    number_cells$table
    
  })
  
  
  #################################################################################################################################################################
  
  #AVOID REPEATING CELLS AND SHOWING CLASSIFIER PROBABILITIES
  
  #The aim is that after clicking one of the buttons, a new cell will be displayed, therefore we will updated the image_number value every time the user clicks
  #one of the buttons. In principle several events can be added to the observeEvent function. The random sampling of cells will only be done for the initial 100-200 cells.
  #Afterwards a filter towards cells which are difficult to classifiy by the model will be applied. 
  #When we click a labelling button, the image file will be eliminated so that the cell is not shown again.
  
  observeEvent(c(input$Yes,input$No, input$skip),{
    
    for (number in sample(classified_cells$remaining_cells, 1000)){
      
      if(nrow(classified_cells$df)<=100){
        
        image_number$vector <- number
        
        classified_cells$remaining_cells <- classified_cells$remaining_cells[!classified_cells$remaining_cells %in% image_number$vector]
        
        classified_cells$used_numbers <- c(classified_cells$used_numbers, image_number$vector)
        
        break
        
        #When the first 100 cells have been randomly sampled and classified, the filter for cells which are difficult for the model to classify will be applied. 
      }
      
      else{
        
        #A random cell will be picked. Then the class will be predicted by the classifier. Each class gets a probability and then the cell is assigned to the class with
        # the higher probability. 
        #We need a for loop in case the classifier predicts the image with high accuracy and then is not shown to the user. 
        
        cell <- readRDS(paste0(images_folder, list_images_files[number]))$features 

        classifier$probabilities <- predict(classifier$model, cell, type = "prob") %>% as.data.frame()
        
        #the function order gives indexes to the columns based on the value. We will order the groups so that the groups with higher probability are in the first two places. 
        
        classifier$probabilities <- classifier$probabilities[,order(classifier$probabilities, decreasing = TRUE)]
        
        #Only if the difference between the probabilities of the classes with hgiher chances is lower than 0.3 the cell will be chosen and shown to the user. 
        #Otherwise another cell will be sampled and the filtering will be applied. 
        
        if(classifier$probabilities[1,1]-classifier$probabilities[1,2]<0.4){
          
          image_number$vector <- number
          
          #Only if the difference in class probabilities is surpassed the image will be shown. 
          
          classified_cells$remaining_cells <- classified_cells$remaining_cells[!classified_cells$remaining_cells %in% image_number$vector]
          
          classified_cells$used_numbers <- c(classified_cells$used_numbers, image_number$vector)
          
          break
        }
        
      }
      
    }
    
 
  })
  
  
  ######################################################################################################################################################################
  
  #CLASS PROBABILITIES TABLE
  
  #Next to the images we will show for the displayed cell which are the probabilities that the model assigns to each of the classes. This could help improve the accuracy
  #of the model. For instance if the model assigns 0.5 % of probability to one class and 0.4 % to another the most rational decision will be to decide between the two most
  # likely groups.
  
  output$class_probabilities <- renderTable({

    if(nrow(classified_cells$df)>100){

      probablilities <- t(classifier$probabilities) %>% as.data.frame() %>% rename(probability = 1) %>% rownames_to_column(var = "class")

      probablilities
    }
  })
  
  
#######################################################################################################################################################################
  
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
  
 ######################################################################################################################################################################
  
  #DISPLAYED CELL INFORMATION
  
  #We will show the screen, plate, well, field, number and dist.10.nn of the cell that is currently shown to the user. this could be useful to associate phenotypes to 
  #specific conditions or time points. 
  
  output$cell_info <- renderTable({
    
    readRDS(paste0(images_folder, list_images_files[image_number$vector]))$features %>% select(screen, plate, well, field, number, dist.10.nn)
    
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
  
  
  observeEvent(c(input$Yes,input$No),{

    if(nrow(classified_cells$df)>10 & (number_cells$table[3,"number_cells"]/50) %%1 ==0){

      training_data <- classified_cells$df %>% as.data.frame() %>% select(good_segm, numeric_features) %>% mutate(good_segm = factor(good_segm))
      
      classifier$model <- randomForest(good_segm~., importance = TRUE, data = training_data, ntree = 500)

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
  
  observeEvent(c(input$Yes, input$No, input$skip),{
    
    updateSliderInput(session, "dapi_adjust", value = 2)
    
    updateSliderInput(session, "actin_adjust", value = 2)
    
    updateSliderInput(session, "tubulin_adjust", value = 2)
    
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
         ggtitle("Confusion matrix (proportion of cells)")

       ggplotly(heat_map)
     }

   })

  #################################################################################################################################################################
  
  

   # The importance of the features for the classifier will be shown in one of the tabs
   # 
   # output$features <- renderPlot({
   # 
   #   varImpPlot(classifier$model)
   # })
  

 }

shinyApp(ui, server)
