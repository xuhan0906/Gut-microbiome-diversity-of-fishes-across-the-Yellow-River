library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(networkD3)
library(htmlwidgets)
matrix_input<-read.delim("assembly_deterministic_stochastic.txt")
a1<-read.delim('five_processes.txt')
#
matrix_input_t<-gather(matrix_input,key='DLMC',value='AREA',-name) %>% filter(!is.na(AREA))
a2<-gather(a1,key='DLMC',value='AREA',-name)
matrix_input_t<-rbind(matrix_input_t,a2)
#
#matrix_input_t$DLMC<-paste(" ",matrix_input_t$DLMC)
#
nodes_input<-c(matrix_input_t$name,matrix_input_t$DLMC) %>% unique() %>% as.data.frame()
names(nodes_input)<-'name'
#
nodes_input$id<- 0:(nrow(nodes_input)-1)
matrix_input_t$sourceid<-nodes_input$id[match(matrix_input_t$name,nodes_input$name)]
matrix_input_t$endid<-nodes_input$id[match(matrix_input_t$DLMC,nodes_input$name)]
#
#nodes_input$col[nodes_input$id<4]<-"col1"
#nodes_input$col[nodes_input$id>3]<-'col2'
nodes_input$col<-paste0("col",nodes_input$id)
sn<-sankeyNetwork(Links = matrix_input_t,Nodes = nodes_input,Source = "sourceid",
                  Target = "endid",Value = "AREA",nodeWidth = 30,
                  nodePadding = 14,
                  NodeGroup = "col",LinkGroup = "name",fontSize = 13.8,
                  margin = list("right"=100)
                )
#move

onRender(sn,
         '
         function(el, x) {
         var sankey = this.sankey;
         var path = sankey.link();
         var nodes = d3.selectAll(".node");
         var link = d3.selectAll(".link")
         var width = el.getBoundingClientRect().width - 40;
         var height = el.getBoundingClientRect().height - 40;
         
         window.dragmove = function(d) {
         d3.select(this).attr("transform", 
         "translate(" + (
         d.x = Math.max(0, Math.min(width - d.dx, d3.event.x))
         ) + "," + (
         d.y = Math.max(0, Math.min(height - d.dy, d3.event.y))
         ) + ")");
         sankey.relayout();
         link.attr("d", path);
         };
         
         nodes.call(d3.drag()
         .subject(function(d) { return d; })
         .on("start", function() { this.parentNode.appendChild(this); })
         .on("drag", dragmove));
         }
         '
)

require(htmlwidgets)
saveWidget(sn, file="name_of_your_file.html")

require(webshot)
webshot("file:///C:/Users/hp/Desktop/name_of_your_file.html", "name_of_your_pdf.pdf")
