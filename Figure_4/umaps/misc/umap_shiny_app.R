library(tidyverse)
library(dplyr)
library(dbplyr)
library(umap)
library(cluster)
library(e1071)
library(RColorBrewer)
library(plotly)
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(uwot)
library(ggpubr)
library(ggrastr)


#In this app I am going to show a umap with the distribution of a subsample of cells colored by the state. 
#I have retrieved the single cell images for each of the classified cells. 
#I will include a few conditions which show interesting group proportions compared to the average. 

#Next to the UMAP I will show a stacked bar plot with the group proportions. 

#This table contains the normalized feature values as well as the embedded dimensions. 
#A handful of interesting conditions were selected, including RLUC controls. 

#It includes the option of saving the UMAP plot alone or to save and arranged plot of the UMAP and stacked barplot.


data_sc <- sample_data_complete

#In order to colour the heat maps I will round the values to 2 and -2 in order to see the differences more clearly.

round_values <- function(values){
  
  new_values <- if_else(values > 2, 2, values)
  
  new_values <- if_else(new_values < (-2), -2, new_values)
  
}

#theme for the UMAP

load("objects/umap_theme.RData")

#theme for the bar plot

load("objects/bar_plot_theme.RData")

#This table contains the proportion of cells in the groups by condition. It will be used to make the stacked
#bar plot. 

data_prop <- sample_grp_prop


#####################################################################################################################################################################



ui <- fluidPage(theme = "journal",
                fluidRow(
                  
                  
                  #There are four drop down bars. The first one is to select the query gene.
                  
                  column(2, selectInput(inputId = "query", label = "Query gene", 
                                        choices = sort(unique(data_sc$query_name)), 
                                        selected = "RLUC")),
                  
                  #The second one is to select the target gene. The options will be updated based on the query gene selected. 
                  
                  column(2, selectInput(inputId = "target", label = "Target gene", 
                                        choices = sort(unique(data_sc$target_name)),
                                        selected = "RLUC")),
                
                  #The last one allows to select the feature by which the points in the map will be coloured. 
                  
                  column(2, selectInput(inputId = "feature", label = "Feature",
                                        choices = c("group", sort(imp_features)))),
                  
                  #Since I can colour the UMAP by feature values, it will be interesting to have a button to save the UMAP that show relevant information. 
                  
                  column(1, actionButton(inputId = "save_umap", label = "Save UMAP")),
                  
                  column(1, actionButton(inputId = "save_arranged", label = "Save arranged"))
                  
                  
                  
                  ),
                
                mainPanel(tabsetPanel(
                  
                  #The first tab will show the UMAP for the condition selected colored by state as a plotly. On the side of the plot 
                  #I will show the image corresponding to the point at which we're hovering. 
                  #On top of the single cell image I'll indicate the query and target genes which were knockdown.
                  
                  tabPanel("UMAP", fluidRow(column(10, plotOutput("umap", height = 700, width = 800)),
                                            column(2, plotOutput("bar_plot", height = 610, width = 200))))
                           
                  
                  
                )))




server <- function(input, output, session) {
  
  #####################################################################################################################################################################
  
  
  #Slider target genes RLUC
  
  #For RLUC I only have the negative control conditions, therefore I will update the drop down menu of target genes.
  
  observeEvent(input$query,{
    
    if(input$query == "RLUC"){
      
      updateSelectizeInput(session, inputId = "target", label = "Target gene", choices = "RLUC")
      
      
    }else{
      
      updateSelectizeInput(session, inputId = "target", label = "Target gene", 
                           choices = sort(unique(data_sc$target_name)))
    }
    
    
  })
  
  #browser()
  
  #Reactive objects
  
  #Table with the single cell features
  
  data <- reactiveValues(umap = data_sc, barplot = data_prop)
  
  
  #####################################################################################################################################################################
  
  #Filter dataframe
  
  #I will filter the single cell dataframe for the query, target and time point selected by the user. That is why it has to be a reactive object.

  observeEvent(c(input$query, input$target),{

    if (nrow(data_sc %>% filter(query_name == input$query, target_name == input$target))>0){

      data$umap <- data_sc %>% filter(query_name == input$query, target_name == input$target) %>% ungroup() %>% sample_n(1500,replace = T)
      
      data$barplot <- data_prop %>% filter(query_name == input$query, target_name == input$target)
  
   }

})


  ######################################################################################################################################################################

  #UMAP

  output$umap <- renderPlot({

    if (input$feature == "group"){

      umap_kd_plot <-  ggplot(data$umap, aes(x=dimension1, y= dimension2, color = factor(group)))+
        geom_point(size = 1.5)+
        umap_theme()+
        theme(
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 19),
          plot.title = element_text(size = 24, hjust = 0.5))+
        ggtitle(paste0(input$query, " + ", input$target))+
        guides(color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 6)))+
        scale_color_manual(values = colours)+
        xlab("UMAP dimension 1")+
        ylab("UMAP dimension 2")
      
      umap_kd_plot


    }else{

      new_data <- data$umap %>% mutate_at(.vars = c(input$feature), round_values)

      umap_kd_plot <-  ggplot(new_data, aes_string(x="dimension1", y= "dimension2", color = input$feature))+
        geom_point(size = 1.5)+
        umap_theme()+
        theme(
          legend.position = "bottom",
          legend.title = element_text(size = 19),
          plot.title = element_text(size = 24, hjust = 0.5))+
        ggtitle(paste0(input$query, " + ", input$target))+
        guides(color = guide_colorbar(title.vjust = 0.5,
                                      barwidth = unit(7, "cm"), barheight = unit(0.75, "cm")
                                      ))+
        scale_color_gradient2(low = google_blue, mid = "white", high = google_red, 
                              midpoint = new_data %>% pull(input$feature) %>% median())+
        xlab("UMAP dimension 1")+
        ylab("UMAP dimension 2")
      
      umap_kd_plot



    }
  })
  
###################################################################################################################

#Stacked bar plot

  output$bar_plot <- renderPlot({

    bar_plot <- ggplot(data$barplot, aes(x = context, y = prop*100, fill = group))+
                                 geom_bar(stat = "identity", position = "stack")+
                                 bar_plot_theme()+
                                 scale_y_continuous(expand = c(0,0), limits = c(0,101))+
                                 scale_x_discrete(expand = c(0,0))+
                                 theme(legend.title = element_blank(),
                                       legend.position = "none",
                                       plot.title = element_blank(),
                                       axis.text.x = element_blank(),
                                       axis.ticks.x = element_blank())+
                                 guides(fill=guide_legend(nrow=3,byrow=TRUE))+
                                 scale_fill_manual(values = c(google_yellow,
                                                              google_red,lacoste_green,
                                                              google_green,
                                                              google_blue))+
                                 geom_text(data = data$barplot  %>% filter(prop >0.04),
                                           aes(label = round(prop*100, 1)),
                                           position = position_stack(vjust = .5),
                                           colour = "white", size = 7)+
                                 ylab("Percentage of cells (%)")
      bar_plot


  })

#####################################################################################################################################################################

  #Save the corresponding UMAP

  observeEvent(input$save_umap,{

    if (input$feature == "group"){

      umap_kd_plot <-  ggplot(data$umap, aes(x=dimension1, y= dimension2, 
                                           color = factor(group)))+
                                  geom_point(size = 1.5)+
                                  umap_theme()+
                                  theme(
                                    legend.title = element_blank())+
                                  scale_color_manual(values = c(google_green, lacoste_green,
                                                                google_yellow, 
                                                                google_red, google_blue), 
                                                     labels = c("Enlarged", "Condensed", "Elongated",
                                                                "Irregular nucleus", "Normal"))+
                                  guides(colour = guide_legend(nrow = 2, byrow = T, 
                                                               override.aes = list(size = 6)))+
                                  xlab("UMAP dimension 1")+
                                  ylab("UMAP dimension 2")
      

      }else{new_data <- data$umap %>% mutate_at(.vars = c(input$feature), round_values)

              umap_kd_plot <-  ggplot(new_data, aes_string(x="dimension1", y= "dimension2", 
                                                           color = input$feature))+
                geom_point(size = 1.5)+
                umap_theme()+
                theme(
                  legend.text = element_text(size = 16),
                  legend.position = "bottom",
                  legend.title = element_text(size = 18, vjust = 0.5),
                  legend.spacing.x = unit(1.5, "cm"),
                  axis.title = element_text(size = 13),
                  legend.key.size = unit(1, "cm"),
                  plot.title = element_blank())+
                scale_color_gradient2(low = google_blue, mid = "white", high = google_red, 
                                        midpoint = new_data %>% pull(input$feature) %>% median())+
                xlab("UMAP dimension 1")+
                ylab("UMAP dimension 2")+
                guides(color = guide_colourbar(title.position="top", title.hjust = 0.5))
      }


      ggsave(umap_kd_plot, filename = paste0("plots/app/", input$query, "_", 
                                             input$target, "_", input$feature, ".png"),
             height = 7, width = 7.5)

    
  })
  
  #####################################################################################################################################################################
  
  #Save the arranged plot
  
  observeEvent(input$save_arranged,{
    
    if (input$feature == "group"){
      
      umap_kd_plot <-  ggplot(data$umap, aes(x=dimension1, y= dimension2, 
                                             color = factor(group)))+
                                geom_point_rast(size = 1.5)+
                                umap_theme()+
                                theme(
                                  legend.title = element_blank())+
                                scale_color_manual(values = c(google_green, lacoste_green,
                                                              google_yellow, 
                                                              google_red, google_blue), 
                                                   labels = c("Enlarged", "Condensed", "Elongated",
                                                              "Irregular nucleus", "Normal"))+
                                guides(colour = guide_legend(nrow = 2, byrow = T, 
                                                             override.aes = list(size = 6)))+
                                xlab("UMAP dimension 1")+
                                ylab("UMAP dimension 2")
                              
      bar_plot <- ggplot(data$barplot, aes(x = context, y = prop*100, fill = group))+
                            geom_bar(stat = "identity", position = "stack")+
                            bar_plot_theme()+
                            scale_y_continuous(expand = c(0,0), limits = c(0,101))+
                            scale_x_discrete(expand = c(0,0))+
                            theme(legend.title = element_blank(),
                                  legend.position = "none",
                                  plot.title = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank())+
                            guides(fill=guide_legend(nrow=3,byrow=TRUE))+
                            scale_fill_manual(values = c(google_yellow,
                                                         google_red,lacoste_green,
                                                         google_green,
                                                         google_blue))+
                            geom_text(data = data$barplot  %>% filter(prop >0.04),
                                      aes(label = round(prop*100, 1)),
                                      position = position_stack(vjust = .5),
                                      colour = "white", size = 7)+
                            ylab("Percentage of cells (%)")
      
      
      
      arranged_plot <- ggarrange(bar_plot, umap_kd_plot, ncol = 2, widths = c(1,3.5))
      
      if(input$query == "RLUC"){
        
        arranged_plot_tit <- annotate_figure(arranged_plot, 
                                             top = text_grob(paste0("RLUC RNAi: Isolated cells"), 
                                                                    size = 26, hjust = 0.5))
        
      }else{arranged_plot_tit <- annotate_figure(arranged_plot, 
                                                   top = text_grob(paste0(input$query, " + ", input$target,
                                                                          " RNAi: Isolated cells"), 
                                                    size = 26, hjust = 0.5))}
      

      
      ggsave(arranged_plot_tit, filename = paste0("plots/app/", input$query, "_", 
                                                  input$target, "_arranged.pdf"), 
             height = 8, width = 10)
    }
    
    
  })


}

shinyApp(ui, server)

                  