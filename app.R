## An shiny app for using gene symbol and exon numbers 
## to search matching structure variants, deletion and duplication only, 
## in gnomAD v4 structure database (SVs and CNVs). Default to h38 genome.

library(shiny)
library(bslib)
library(Gviz)
library(rtracklayer)

## MANE and MANE Plus Clinical exon table, v1.4, hg38.
mane_exon = read.table("data/mane1.4_exon", header = TRUE, sep = "\t")

## gnomAD SV (del and dup only) and CNV BED files with selected columns and 
## 1 added to start to be in GenBank style. SV read in below. 
cnv = read.table("data/gnomADcnv", header = TRUE, sep = "\t")

## chain file for hg38 to hg19 liftover.
hg38_19 = import.chain("data/hg38ToHg19.over.chain")


ui <- page_sidebar(
    theme = bs_theme(version = 5, bootswatch = "journal"),
    title = markdown("**gnomADsvSeek**"),
    
    ## input
    sidebar = sidebar(
        width = 300,
        HTML("Search gnomAD structure variant database for matching entries.<br>
                 <br>Deletion and duplication only.<br>
                 <br>MANE and MANE Plus Clinical transcripts only.<br>"),
        HTML('<p>Visit app <a href="https://github.com/wshuchen/gnomADsvSeek"
             target="_blank">GitHub repository</a></p>'),
        textInput("gene", markdown("**Gene name**"),
              value = "",
              placeholder = "TP53 OR tp53"),
        numericInput("exon_from", markdown("**Exon from**"),
                 value = ""),
        numericInput("exon_to", markdown("**Exon to**"),
                 value = ""),
        radioButtons("variant", markdown("**Variant**"), 
                     choices = c("SV", "CNV"),
                     selected = "SV"),
        radioButtons("type", markdown("**Type**"), 
                     choices = c("del", "dup"),
                     selected = "del"), 
        radioButtons("boundary", HTML("<b>Boundary</b>
                                      <br>by exon not first or last"), 
                     choices = c("bound", "unbound"),
                     selected = "unbound"),
        radioButtons("genome", HTML("<b>Genome</b>
                                    <br>h38 to hg19 liftover for the result"), 
                     choices = c("hg19", "hg38"),
                     selected = "hg38"),
        actionButton("search", h5("Search", style = "color: green;")),
        sliderInput("zoom", HTML("<b>Zoom</b>"),
                    min = 0, max = 3, step = 0.5,
                    value = 0),
        actionButton("clear", h5("Clear result", style = "color: red;"))
    ),
    
    ## output
     card(card_header(HTML("<b>Query gene</b>")),
          height = "350px",
          tableOutput("gene_exon")),
     textOutput("gene_info"),
     uiOutput("error"),
     card(height = "300px",
          plotOutput("gene_viz")),
     card(card_header(HTML("<b>Overlapping intervals</b>")),
          plotOutput("matching_sv_viz"),
          height = "500px"),
     card(tableOutput("matching_sv"),
          height = "300px")
)

server <- function(input, output, session) {

    ## query gene data and view
    gene_df = reactive({
        req(input$gene)
        gene_df = mane_exon[mane_exon$symbol == toupper(input$gene), ]
        gene_info = paste(gene_df$gene[1],
                          gene_df$chrom[1],
                          gene_df$start[1], "-", gene_df$end[length(gene_df$end)],
                          gene_df$transcript[1],
                          gene_df$exon[length(gene_df$exon)], "exons")
        gene_track = GeneRegionTrack(gene_df, 
                                     genome = "hg38",
                                     name = gene_df$gene[1])
        list(gene_df, gene_info, gene_track)
    })
    
    ## query CNV data and view
    cnv_df = reactive({
        exon_from = input$exon_from
        exon_to = input$exon_to
        gene_df = gene_df()[[1]]
        gene_track = gene_df()[[3]]
        
        cnv_df = gene_df[exon_from:exon_to, ]
        cnv_track = GeneRegionTrack(cnv_df, name = "CNV",
                                    background.panel = "#FFFEDB",
                                    background.title = "red")
        cnv_track = HighlightTrack(list(gene_track, cnv_track),
                                    start = min(cnv_df$start),
                                    end = max(cnv_df$end))
        list(cnv_df, cnv_track)
    })
    
    ## Search the database for matching SVs.
    match_cnv = reactive({
        exon_from = input$exon_from
        exon_to = input$exon_to
        
        gene_df = gene_df()[[1]]
        cnv_df = cnv_df()[[1]]
        cnv_gr = GRanges(seqnames = cnv_df$chrom[1],
                               ranges = IRanges(min(cnv_df$start),
                                                max(cnv_df$end)),
                               strand = cnv_df$strand[1])

        # Get the interval containing the query CNV.
        chr = unique(cnv_df$chrom)
        sv_file = paste0("data/SV/gnomADsv_", chr)
        sv = read.table(sv_file, header = TRUE, sep = "\t")
        if (input$variant == "SV") {
            sv_select = sv[sv$chrom == cnv_df$chrom[1] &
                               sv$type == toupper(input$type), ]    
        } else {
            sv_select = cnv[cnv$chrom == cnv_df$chrom[1] &
                   cnv$type == toupper(input$type), ] 
        }
        sv_select_gr = GRanges(seqnames = sv_select$chrom,
                                ranges = IRanges(sv_select$start, sv_select$end))
        if (input$variant == "SV") {
            mcols(sv_select_gr) = sv_select[4:9]    
        } else {
            mcols(sv_select_gr) = sv_select[4:10]
        }
        
        sv_hits = findOverlaps(cnv_gr, sv_select_gr, type = "within")
        hits_index = subjectHits(sv_hits)
        sv_found = sv_select_gr[hits_index]

        # Restrict the range to within the exon not first or last,
        # open in both ends (exon 1 and last exon), with strand consideration.
        if (input$boundary == "bound") {
            if (exon_from > 1) {
                if (gene_df$strand[1] == "+") {
                    lower_bound = gene_df$end[exon_from-1]
                    sv_found = sv_found[start(sv_found) >= lower_bound[1]]
                } else {
                    lower_bound = gene_df$start[exon_from-1]
                    sv_found = sv_found[end(sv_found) <= lower_bound[1]]
                }
            }
            if (exon_to < max(gene_df$exon)) {
                if (gene_df$strand[1] == "+") {
                    upper_bound = gene_df$start[exon_to+1]
                    sv_found = sv_found[end(sv_found) <= upper_bound[1]]
                } else {
                    upper_bound = gene_df$end[exon_to+1]
                    sv_found = sv_found[start(sv_found) >= upper_bound[1]]
                }           
            }
        }
        
        # Add back the frequency data removed for smaller file size.
        if (input$variant == "SV") {
            mcols(sv_found)$frequency = format(sv_found$ac/sv_found$an,
                                               scientific = FALSE)
        } else {
            mcols(sv_found)$frequency = format(sv_found$sc/sv_found$sn, 
                                               scientific = FALSE)
        }
        # Name in uppercase in gnomAD browser
        sv_found$name = toupper(sv_found$name) 
        
        # Add gene symbol(s) to the data, so we know the gene(s) an interval overlaps.
        # This helps to see the actual match because we have both ends open.
        # Note that original SV BED file does not include gene symbol, while
        # CNV file has a "genes" column with only one symbol in a row, 
        # which may not include all genes in the interval. For CNV result, 
        # "gene" column is the original "genes" column; the new "genes" column
        # gathers all the genes (with symbols up to 3) overlapping an interval.
        get_symbol = function(sv_name) {
            sv_gr = sv_found[sv_found$name == sv_name]
            mane = mane_exon[mane_exon$chrom == seqlevels(sv_gr), ]
            mane_gr = GRanges(seqnames = mane$chrom, 
                              ranges = IRanges(start = mane$start,
                                               end = mane$end))
            mcols(mane_gr)$symbol = mane$symbol
            ol = findOverlaps(sv_gr, mane_gr)
            mane_ol = mane_gr[subjectHits(ol)]
            symbols = unique(mane_ol$symbol)
            if (length(symbols) <= 3) {
                return(paste(symbols, collapse = ", "))
            } else {
                n = length(symbols)
                symbols = paste(symbols[1:3], collapse = ", ")
                return(paste(symbols, "and other", n, "MANE genes"))
            }
        }

        if (length(sv_found) >= 1) {
             symbols = lapply(sv_found$name, function(x) (get_symbol(x)))
            if (input$variant == "SV") {
                sv_found$genes = unlist(symbols) 
            } else {
                sv_found$genes = unlist(symbols)
            }
        }
        sv_found
    })
           
    ## Outputs
    # Display gene info and view when gene name is provided.
    observe({
        req(toupper(input$gene) %in% mane_exon$symbol)
        gene_df = gene_df()[[1]]
        gene_info = gene_df()[[2]]
        gene_track = gene_df()[[3]]
        ax <- GenomeAxisTrack()
        output$gene_viz = renderPlot({plotTracks(list(ax, gene_track))})
        output$gene_exon = renderTable({gene_df},
                                       striped = TRUE, hover = TRUE, bordered = TRUE,
                                       align = "c")
        output$gene_info = renderText(gene_info)
    })
    
    # Display gene and CNV info; validate the input exon range.
    observe({
        gene_df = gene_df()[[1]]
        gene_track = gene_df()[[3]]
        req(input$exon_to >= input$exon_from & input$exon_to <= max(gene_df$exon))

        output$error = renderUI({
            validate(need(input$exon_to >= input$exon_from,
                          "Error: second exon number must >= first exon number."),
                     need(input$exon_to <= max(gene_df$exon),
                          "Error: last exon number out of range."),
                     paste("You enter", "gene:", toupper(input$gene), 
                           "exon", input$exon_from,
                           "to", input$exon_to))
        })
        cnv_df = cnv_df()[[1]]
        cnv_track = cnv_df()[[2]]
        ax <- GenomeAxisTrack()
        output$gene_viz = renderPlot({plotTracks(list(ax, cnv_track))})
    })
    
    # Display matching SVs if found.
    observeEvent(c(input$search, input$zoom), {
        gene_df = gene_df()[[1]]
        req(input$exon_to >= input$exon_from & input$exon_to <= max(gene_df$exon))

        cnv_track = cnv_df()[[2]]
        ax <- GenomeAxisTrack()
        
        sv_found = match_cnv()
        if (length(sv_found) == 0) {
            output$error = renderUI({p(paste("No overlapping SV found."), style = "color: red;")})
        }  
        if (length(sv_found) >= 1) {
            SV = paste0(as.character(sv_found$name), collapse = ", ")
            output$error = renderUI({p(paste("Overlapping SVs:", SV), style = "color: green;")})
            
            # Display in the scale of longest interval,
            # which help to judge the match.
            sv_found = as.data.frame(rev(sv_found))
            colnames(sv_found)[1] = "chrom"
            sv_pos = sv_found[, 1:3]
            colnames(sv_pos) = c("chrom", "start", "end")
            pos = rbind(gene_df[, 1:3], sv_pos)

            found_track = lapply(1:nrow(sv_found),
                                 function(x) {GeneRegionTrack(sv_found[x, ],
                                              name = as.character(sv_found[x, ]$name))
                            })
            
            # Zoom: extend the range at half of the gene length both sides at 1
            if (input$zoom == 0) {
                output$matching_sv_viz = renderPlot({
                                plotTracks(c(list(ax, cnv_track), found_track),
                                           from = min(pos$start), to = max(pos$end))
                })
            } else {
                fold = input$zoom
                if (gene_df$strand[1] == "+") {
                    gene_length = max(gene_df$end) - min(gene_df$start) + 1
                    output$matching_sv_viz = renderPlot({
                        plotTracks(c(list(ax, cnv_track), found_track),
                            from = min(gene_df$start) - round(gene_length * (3-fold)),
                            to = max(gene_df$end) + round(gene_length * (3-fold)))
                    })
                } else {
                    gene_length = max(gene_df$start) - min(gene_df$end) - 1
                    output$matching_sv_viz = renderPlot({
                        plotTracks(c(list(ax, cnv_track), found_track),
                            from = min(gene_df$end) - round(gene_length * (3-fold)),
                            to = max(gene_df$start) + round (gene_length * (3-fold)))                    
                    })
                }
            }
            output$matching_sv = renderTable({sv_found},
                                             striped = TRUE, hover = TRUE, 
                                             bordered = TRUE, align = "c")
        }
    })
    
    # Liftover may break a large interval into pieces with different gaps. 
    # We would merge all pieces from a sv into one (hopefully) when that happens.
    # We would also get rid of any pieces mapped to other chromosomes.
    observeEvent(input$genome == "hg19", {
        sv_found = match_cnv()
        req(length(sv_found) >= 1)
        chrom = seqlevels(sv_found)
        sv_found = unlist(liftOver(sv_found, hg38_19))
        sv_found = reduce(split(sv_found, sv_found$name), 
                          min.gapwidth = 10000000)
        sv_found = as.data.frame(sv_found)
        colnames(sv_found) = c("", "name", "chrom", "start", "end", "width", "strand")
        sv_found = sv_found[sv_found$chrom == chrom, ]
        output$matching_sv = renderTable({sv_found},
                                         striped = TRUE, hover = TRUE, 
                                         bordered = TRUE, align = "c")
    })
        
    # Clear search result when needed.
    observeEvent(input$clear, {
        output$matching_sv_viz = renderPlot({})
        output$matching_sv = renderTable({})
        output$error = renderUI({})
    })
}

shinyApp(ui, server)