## A shiny app for using gene symbol and exon number 
## to search matching exon deletion and duplication in gnomADv4  
## structure variant data (SVs and CNVs). Default to h38 genome.
## Use processed gnomAD SV (del and dup only) and CNV BED files 
## with selected columns and 1 added to start to be in GenBank style.

library(shiny)
library(bslib)
library(GenomicRanges)
library(Gviz)
library(DT)

ui <- page_sidebar(
    theme = bs_theme(version = 5, bootswatch = "journal"),
    title = tagList(
        shiny::tags$span(
            shiny::tags$i("gnomADsvSeek"), 
            style = "margin-left: 30px; font-weight: bold; font-size: 1.5rem;"),
        shiny::tags$span(
            shiny::tags$i("Searching exon del/dup in gnomAD structure variant data"), 
            style = "margin-top: 5px; font-weight: bold; font-size: 1rem;"),
        shiny::tags$a(
            href = "https://github.com/wshuchen/gnomADsvSeek.git", 
            "GitHub repository",
            target = "_blank",
            style = "margin-top: 5px; margin-right: 30px; 
                    font-size: 1rem; color: white; 
                    text-decoration: underline;"
        )
    ),
    
    ## Input
    sidebar = sidebar(
        width = 300,
        p("hg38 MANE transcripts", 
          style = "color: green; font-weight: bold; font-size: 18px;"),
        p("(e.g. PKD1 deletion of exons 22-30)", style = "font-size: 1rem;"),
        p("Please verify the result.", style = "color: green;"),
        
        textInput("gene", HTML("<b>Gene name</b>"),
                  value = "",
                  placeholder = "PKD1 OR pkd1"),
        numericInput("exon_from", HTML("<b>Exon from</b>"),
                     value = ""),
        numericInput("exon_to", HTML("<b>Exon to</b>"),
                     value = ""),
        radioButtons("variant", HTML("<b>Variant</b>"), 
                     choices = c("SV", "CNV"),
                     selected = "SV",
                     inline = TRUE),
        radioButtons("type", HTML("<b>Type</b>"), 
                     choices = c("del", "dup"),
                     selected = "del",
                     inline = TRUE),  
        radioButtons("limit", HTML("<b>Limit</b>
                                      <br>by exon not first or last"), 
                     choices = c("yes", "no"),
                     selected = "no", 
                     inline = TRUE),
        radioButtons("genome", HTML("<b>Genome</b>
                                    <br>h38 to hg19 liftover for the result"), 
                     choices = c("hg19", "hg38"),
                     selected = "hg38",
                     inline = TRUE),
        div(style = "text-align: center;",
            actionButton("search", "Search", 
                         style = "background-color: #01774e;
                                 font-size: 18px; color: white;
                                 height: 40px; width: 120px;")
        ),
        sliderInput("zoom", HTML("<b>Zoom in</b>"),
                    min = 0, max = 3, step = 0.5, value = 0),
        radioButtons("get_met", HTML("<b>getMet</b>
                                      <br>retrieve methionine position"),
                     choices = c("yes", "no"),
                     selected = "no",
                     inline = TRUE),
        div(style = "text-align: center;",
            actionButton("clear", "Clear", 
                         style = "background-color: #01774e;
                                 font-size: 18px; color: white;
                                 height: 40px; width: 120px;")
        )
    ),
    
    ## Output
    uiOutput("gene_info"),
    card(card_header(HTML("<b>Query gene</b>")),
         height = "400px",
         DTOutput("gene_exon")),
    card(card_header(HTML("<b>Query CNV</b>")),
         height = "200px",
         tableOutput("cnv_info")),
    uiOutput("error"),
    card(width = "100%",
         height = "auto",
         plotOutput("gene_viz")),
    card(tableOutput("matching_sv"),
         height = "300px")
)

server <- function(input, output, session) {
    
    ## MANE and MANE Plus Clinical exon table, v1.4, hg38.
    mane_tx = readRDS("rdsData/mane1.4_transcript.rds")

    ## Query gene data and view
    gene_df = reactive({
        req(input$gene)
        gene_df = mane_tx[mane_tx$symbol == toupper(input$gene), ]
        rownames(gene_df) = NULL
        exon_number = as.numeric(gene_df$exon[grep("[0-9.*]", gene_df$exon)])
        exon_number = exon_number[length(exon_number)]
        gene_info = paste(toupper(input$gene), 
                          paste0(unique(gene_df$chrom), ":",
                                 gene_df$start[1], "-", 
                                 gene_df$end[length(gene_df$end)]),
                                 unique(gene_df$transcript),
                                 exon_number, "exons")
        
        options(ucscChromosomeNames = FALSE)
        gene_track = GeneRegionTrack(gene_df, 
                                     genome = "hg38",
                                     name = toupper(input$gene))
        list(gene_df, gene_info, gene_track, exon_number)
    })
    
    ## Query CNV data and view
    cnv_df = reactive({
        exon_from = input$exon_from
        exon_to = input$exon_to
        gene_df = gene_df()[[1]]
        gene_track = gene_df()[[3]]
        
        cnv_df = gene_df[gene_df$exon %in% exon_from:exon_to, ]
        exon_range = paste0(exon_from, "-", exon_to)
        cnv_cds = cnv_df$feature =="cds"
        cds_size = sum(cnv_df[cnv_cds, ]$width)
        cds_start = min(as.numeric(gsub("\\-.*", "", cnv_df[cnv_cds, ]$CDS)))
        cds_end = max(as.numeric(gsub(".*\\-", "", cnv_df[cnv_cds, ]$CDS)))
        cds_range = paste0(cds_start, "_", cds_end)
        gene_cds = gene_df$feature =="cds"
        percentage = round(cds_size/sum(gene_df[gene_cds, ]$width)*100, 1)
        frame = ifelse(cds_size %% 3 == 0, "in frame", "out of frame")
        cnv_info = data.frame(chrom = unique(cnv_df$chrom),
                              start = min(cnv_df$start),
                              end = max(cnv_df$end),
                              width = sum(cnv_df$width),
                              "exon range" = exon_range,
                              "CDS range" = cds_range,
                              "CDS size" = cds_size,
                              "CDS %" = percentage,
                              frame = frame, 
                              check.names = FALSE)
        
        options(ucscChromosomeNames = FALSE)
        cnv_track = GeneRegionTrack(cnv_df, 
                                    genome = "hg38",
                                    name = "CNV",
                                    background.panel = "#F0FFFF",
                                    background.title = "#d65b00")
        cnv_track = HighlightTrack(list(gene_track, cnv_track),
                                   start = min(cnv_df$start),
                                   end = max(cnv_df$end))
        list(cnv_df, cnv_info, cnv_track)
    })
    
    ## Search the database for matching SVs.
    match_cnv = reactive({
        exon_from = input$exon_from
        exon_to = input$exon_to
        
        gene_df = gene_df()[[1]]
        exon_number = gene_df()[[4]]
        
        cnv_df = cnv_df()[[1]]
        cnv_gr = GRanges(seqnames = cnv_df$chrom[1],
                         ranges = IRanges(min(cnv_df$start),
                                          max(cnv_df$end)),
                         strand = cnv_df$strand[1])
        
        # Get the interval containing the query CNV.
        cnv_chrom = unique(cnv_df$chrom)
        if (input$variant == "SV") {
            sv_file = paste0("rdsData/SV/gnomADsv_", cnv_chrom, ".rds")
            sv = readRDS(sv_file)
            sv_select = sv[sv$chrom == cnv_chrom &
                               sv$type == toupper(input$type), ]    
        } else {
            cnv = readRDS("rdsData/gnomADcnv.rds")
            sv_select = cnv[cnv$chrom == cnv_chrom &
                                cnv$type == toupper(input$type), ] 
        }
        sv_select_gr = GRanges(seqnames = sv_select$chrom,
                               ranges = IRanges(sv_select$start, sv_select$end))
        mcols(sv_select_gr) = sv_select[4:length(sv_select)]    
        
        sv_hits = findOverlaps(cnv_gr, sv_select_gr, type = "within")
        hits_index = subjectHits(sv_hits)
        sv_found = sv_select_gr[hits_index]
        
        # Restrict the range to within the exon not first or last,
        # open in both ends (exon 1 and last exon), with strand consideration.
        # Prepare a exon range table from the gene data, merge utr in exon
        # if they have the same exon number. 
        exons = gene_df[, c("chrom", "start", "end", "strand", "exon")]
        starts = aggregate(start ~ exon, data = exons, min)
        ends = aggregate(end ~ exon, data = exons, max)
        start_end = merge(starts, ends)
        start_end$chrom = c(unique(exons$chrom))
        start_end$strand = c(unique(exons$strand))
        exon_range = start_end[, c("chrom", "start", "end", "strand", "exon")]

        if (input$limit == "yes") {
            if (exon_from > 1) {
                if (exon_range$strand[1] == "+") {
                    start_limit = exon_range$end[exon_from-1]
                    sv_found = sv_found[start(sv_found) >= start_limit[1]]
                } else {
                    start_limit = exon_range$start[exon_from-1]
                    sv_found = sv_found[end(sv_found) <= start_limit[1]]
                }
            }
            
            if (exon_to < exon_number) {
                if (exon_range$strand[1] == "+") {
                    end_limit = exon_range$start[exon_to+1]
                    sv_found = sv_found[end(sv_found) <= end_limit[1]]
                } else {
                    end_limit = exon_range$end[exon_to+1]
                    sv_found = sv_found[start(sv_found) >= end_limit[1]]
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
            mane = mane_tx[mane_tx$chrom == seqlevels(sv_gr), ]
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
    
    ## Get methionine position
    ## plot protein sequence as a matrix. Function courtesy of Google Gemini.
    plot_protein <- function(seq, target = "M", width = 60, col = "#FF4D4D") {
        aa <- unlist(strsplit(seq, ""))
        rows <- ceiling(length(aa) / width)
        mat <- matrix(c(aa, rep("", rows * width - length(aa))), 
                      ncol = width, byrow = TRUE)
        ticks <- unique(c(seq(1, width, by = 10), width))
        
        par(mar = c(0.3, 2.2, 1.2, 0.3)) # Tight margins
        plot(1, type = "n", xlim = c(0.5, width + 0.5), 
             ylim = c(0.5, rows + 0.8), 
             xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
             bty = "n", asp = 1, yaxs = "i", xaxs = "i")
        
        text(ticks, rows + 1.0, labels = ticks, font = 2, cex = 1.0, xpd = TRUE)
        
        for (r in 1:rows) {
            y <- rows - r + 1
            valid <- mat[r, ] != ""
            rect((1:width)[valid] - 0.5, y - 0.5, (1:width)[valid] + 0.5, y + 0.5, 
                 col = ifelse(mat[r, valid] %in% target, col, "#EAEAEA"), 
                 border = "white")
            text((1:width)[valid], y, mat[r, valid], font = 2, cex = 0.75)
            text(-0.2, y, labels = (r - 1) * width + 1, font = 2, cex = 1.0, 
                 adj = c(1, 0.5), xpd = TRUE)
        }
    }
    
    get_met = reactive({
        mane_sum = readRDS("rdsData/mane1.4_summary.rds")
        gene_df = gene_df()[[1]]
        
        transcript = unique(gene_df$transcript)
        protein = mane_sum[mane_sum$transcript == transcript, ]$protein
        seq <- refseq_AAseq(protein)[[1]]
        met_pos = matchPattern("M", seq)
        met_pos = as.data.frame(met_pos)$start
        fasta = as.character(seq)
        fasta_length = nchar(fasta)
        fasta_url = paste0("https://ncbi.nlm.nih.gov/protein/", protein, "?report=fasta")
        fasta_link = paste0('<a href=', fasta_url, ' target="_blank"', '>(fasta)</a>')
        size_before = round((met_pos-1)/fasta_length*100, 1)
        met_pos = paste0(met_pos, "(", size_before, ")", collapse = ", ")
        list(met_pos, fasta, protein, fasta_link)
    })
    
    ## Outputs
    # Display gene info when gene name is provided.
    observe({
        req(toupper(input$gene) %in% mane_tx$symbol)
        gene_df = gene_df()[[1]]
        gene_info = gene_df()[[2]]
        gene_track = gene_df()[[3]]
        ax <- GenomeAxisTrack()
        output$gene_viz = renderPlot({plotTracks(list(ax, gene_track))})
        output$gene_exon = renderDT({datatable(gene_df, 
                                options = list(
                                    columnDefs = list(
                                        list(className = 'dt-center', 
                                             targets = "_all")),
                                    searchHighlight = TRUE), 
                                    rownames = FALSE
                                )
                            })
        output$gene_info = renderUI(p(gene_info, style = "color: green;"))
    })
    
    # Display gene and CNV info; validate the input exon range.
    observe({
        gene_df = gene_df()[[1]]
        gene_info = gene_df()[[2]]
        gene_track = gene_df()[[3]]
        exon_number = gene_df()[[4]]
        
        req(input$exon_to >= input$exon_from & input$exon_to <= exon_number)
        output$error = renderUI({
            validate(need(input$exon_to >= input$exon_from,
                          "Error: second exon number must >= first exon number."),
                     need(input$exon_to <= exon_number,
                          "Error: last exon number out of range."),
                     paste("You enter", "gene:", toupper(input$gene), 
                           "exon", input$exon_from,
                           "to", input$exon_to))
        })
        
        cnv_df = cnv_df()[[1]]
        cnv_info = cnv_df()[[2]]
        cnv_track = cnv_df()[[3]]
        ax <- GenomeAxisTrack()
        output$gene_viz = renderPlot({plotTracks(list(ax, cnv_track))})
        output$cnv_info = renderTable({cnv_info},
                                       striped = TRUE, hover = TRUE, 
                                       bordered = TRUE, align = "c", digit = 1)
    })
    
    # Display matching SVs if found.
    observeEvent(c(input$search, input$zoom), {
        gene_df = gene_df()[[1]]
        exon_number = gene_df()[[4]]
        req(input$exon_to >= input$exon_from & 
                input$exon_to <= exon_number)
        
        cnv_df = cnv_df()[[1]]
        cnv_info = paste0(cnv_df$chrom[1], ":", 
                          min(cnv_df$start), "-", max(cnv_df$end),
                          input$type)
        
        cnv_track = cnv_df()[[3]]
        ax <- GenomeAxisTrack()
        
        sv_found = match_cnv()
        if (length(sv_found) == 0) {
            output$error = renderUI({p(paste("No overlapping SV found."), 
                                       style = "color: red;")})
        }  
        if (length(sv_found) >= 1) {            
            sv_found = as.data.frame(rev(sv_found))
            sv_found = sv_found[order(sv_found$width), ]
            colnames(sv_found)[1] = "chrom"
            
            # Add a link to gnomAD SVs and CNVs database for the matches.
            # https://gnomad.broadinstitute.org/variant/
            # DEL_CHR16_7A0D76E1?dataset=gnomad_sv_r4
            # 211102__DEL?dataset=gnomad_cnv_r4
            
            if (input$variant == "SV") {
                sv_url = paste0("https://gnomad.broadinstitute.org/variant/", 
                            sv_found$name, "?dataset=gnomad_sv_r4")
            }
            if (input$variant == "CNV") {
                sv_url = paste0("https://gnomad.broadinstitute.org/variant/", 
                                sv_found$name, "?dataset=gnomad_cnv_r4")
            }
            
            sv_link = paste0('<a href=', sv_url, ' target="_blank"', '>open</a>')
            sv_found$gnomAD = sv_link

            # Write out the list of SVs. 
            SV = paste0(as.character(sv_found$name), collapse = ", ")
            output$error = renderUI({
                HTML(paste0(
                    "<span style = 'color: green;'>Query variant - ", cnv_info, "</span>",
                    "<span style = 'color: red;'>Wholely overlapping SVs:</span>",
                    "<span style = 'color: blue;'>", SV, "</span>"))
                                            
            })
            
            # Display in the scale of longest interval,
            # which helps to judge the match.
            sv_pos = sv_found[c("chrom", "start", "end")]
            pos = rbind(gene_df[, 1:3], sv_pos)
            
            found_track = lapply(1:nrow(sv_found),
                                 function(x) {GeneRegionTrack(sv_found[x, ],
                                        name = as.character(sv_found[x, ]$name))
                                 })
            
            # Zoom 1: extend the range at half of the gene length 
            # on both sides
            if (input$zoom == 0) {
                output$gene_viz = renderPlot({
                    plotTracks(c(list(ax, cnv_track), found_track),
                               from = min(pos$start), to = max(pos$end))
                })
            } else {
                fold = input$zoom
                if (gene_df$strand[1] == "+") {
                    gene_length = max(gene_df$end) - min(gene_df$start) + 1
                    output$gene_viz = renderPlot({
                        plotTracks(c(list(ax, cnv_track), found_track),
                               from = min(gene_df$start) - round(gene_length * (3-fold)),
                               to = max(gene_df$end) + round(gene_length * (3-fold)))
                    })
                } else {
                    gene_length = max(gene_df$start) - min(gene_df$end) - 1
                    output$gene_viz = renderPlot({
                        plotTracks(c(list(ax, cnv_track), found_track),
                               from = min(gene_df$end) - round(gene_length * (3-fold)),
                               to = max(gene_df$start) + round (gene_length * (3-fold)))                    
                    })
                }
            }
            output$matching_sv = renderTable({sv_found},
                                         striped = TRUE, hover = TRUE, 
                                         bordered = TRUE, align = "c",
                                         sanitize.text.function = function(x) paste(x) )
        }
    })
    
	# h38 to hg19 coordinate liftover for the result if requested.
    observeEvent(input$genome == "hg19", {
    	library(rtracklayer)
        # Chain file for hg38 to hg19 coordinate liftover.
        hg38_19 = readRDS("rdsData/hg38ToHg19.over.chain.rds")
        
        sv_found = match_cnv()
        req(length(sv_found) >= 1)
        chrom = seqlevels(sv_found)
        mini_width = 0.8 * min(width(sv_found)) 
        sv_found = unlist(liftOver(sv_found, hg38_19))
        
        # Liftover may break a large interval into pieces with different gaps. 
        # We would hopefully merge all pieces from a sv into one.
        sv_found = reduce(split(sv_found, sv_found$name), min.gapwidth = 1e+8)
                          
        # We would also remove any pieces mapped to other chromosomes, 
        # and anything that is less than 80% of the size of the original entries.
        sv_found = as.data.frame(sv_found)
        sv_found = sv_found[, 2:length(sv_found)]
        colnames(sv_found) = c("name", "chrom", "start", "end", "width", "strand")
        sv_found = sv_found[sv_found$chrom == chrom & sv_found$width >= mini_width, ]
        sv_found = sv_found[order(sv_found$width), ]
        sv_list = paste(paste0(sv_found$name, " (", sv_found$chrom, ":", 
                         sv_found$start, "-", sv_found$end, ")"), collapse = ", ")

        output$matching_sv = renderTable({sv_found},
                                         striped = TRUE, hover = TRUE, 
                                         bordered = TRUE, align = "c")
        output$error = renderUI({
                    HTML(paste0(
                        "<span style = 'color: red;'>hg19 lifeover:</span>",
                        "<span style = 'color: blue;'>", sv_list, "</span>"))
        })
    })
    
    # Show methionine position for the corresponding protein sequence.
    observeEvent(input$get_met == "yes", {
        library(Biostrings)
        library(refseqR)
        req(input$gene)
        
        met_pos = get_met()[[1]]
        fasta = get_met()[[2]]
        protein = get_met()[[3]]
        fasta_link = get_met()[[4]]
        # Only display up to 900aa. 
        if (nchar(fasta) > 900) fasta = substr(fasta, 1, 900)
        output$error = renderUI({
            HTML(paste0(
                "<span style='color:red;'>The methionine position (% of size before):</span>",
                "<span style='color:blue;'>", met_pos, "</span>",
                "<span style='color:green;'>Display protein sequence (up to 900aa): ", 
                        protein, " ", fasta_link,"</span>"))
        })
        output$gene_viz = renderPlot({plot_protein(fasta, target = "M", width = 100)})
    })

    # Clear search result when needed.
    observeEvent(input$clear, {
        output$gene_viz = renderPlot({})
        output$matching_sv = renderTable({})
        output$error = renderUI({})
    })
}

shinyApp(ui, server)
